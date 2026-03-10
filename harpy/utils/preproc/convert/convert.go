// preproc_barcodes - renames and records demultiplexed linked-read barcodes
// coming out of Pheniqs. Reads a Pheniqs JSON config and a SAM/BAM file,
// rewrites BX and VX tags on every record that has an RX tag, and writes
// SAM to stdout (no header).
//
// Usage:
//
//	preproc_barcodes <pheniqs.json> <input.bam>
//
// Build:
//
//	go mod init preproc_barcodes
//	go get github.com/biogo/hts@latest
//	go build -ldflags="-s -w" -o preproc_barcodes preproc_barcodes.go
package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ─── JSON schema ─────────────────────────────────────────────────────────────
// We only need the two codec maps from the Pheniqs JSON output.

type codec struct {
	Barcode []string `json:"barcode"`
}

type decoder struct {
	Stagger struct {
		Codec map[string]codec `json:"codec"`
	} `json:"stagger"`
	Segment struct {
		Codec map[string]codec `json:"codec"`
	} `json:"segment"`
}

type pheniqsJSON struct {
	Decoder decoder `json:"decoder"`
}

// ─── barcode reconstruction ───────────────────────────────────────────────────

// reconstructBarcode mirrors the Python reconstruct_barcode function.
// nuc_bc is the value of the RX SAM tag, e.g. "001-001-002-003-004".
// Returns the BX string and a VX validity flag (1=valid, 0=invalid).
func reconstructBarcode(nucBC string, stagger, bc map[string]string) (string, int) {
	seg := strings.SplitN(nucBC, "-", 4)
	if len(seg) != 4 {
		return "A00C00B00D00", 0
	}

	A := "A" + lookup(bc, seg[1])
	B := "B" + lookup(bc, seg[2])
	C := "C" + lookup(bc, seg[3])
	D := "D" + lookup(stagger, seg[0])

	valid := 1
	if A == "A00" || B == "B00" || C == "C00" || D == "D00" {
		valid = 0
	}
	return A + C + B + D, valid
}

func lookup(m map[string]string, key string) string {
	if v, ok := m[key]; ok {
		return v
	}
	return "00"
}

// ─── main ─────────────────────────────────────────────────────────────────────

func main() {
	if len(os.Args) != 3 {
		fmt.Fprintf(os.Stderr, "Usage: preproc_barcodes <pheniqs.json> <input.bam>\n")
		os.Exit(1)
	}
	jsonPath := os.Args[1]
	bamPath := os.Args[2]

	// ── parse JSON ──────────────────────────────────────────────────────────
	jf, err := os.Open(jsonPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening JSON: %v\n", err)
		os.Exit(1)
	}
	var pj pheniqsJSON
	if err := json.NewDecoder(jf).Decode(&pj); err != nil {
		fmt.Fprintf(os.Stderr, "Error parsing JSON: %v\n", err)
		os.Exit(1)
	}
	jf.Close()

	// Build lookup maps: barcode-sequence -> numeric code (prefix stripped).
	stagger := make(map[string]string, len(pj.Decoder.Stagger.Codec))
	for k, v := range pj.Decoder.Stagger.Codec {
		if len(v.Barcode) > 0 {
			stagger[v.Barcode[0]] = strings.TrimPrefix(k, "@")
		}
	}
	bc := make(map[string]string, len(pj.Decoder.Segment.Codec))
	for k, v := range pj.Decoder.Segment.Codec {
		if len(v.Barcode) > 0 {
			bc[v.Barcode[0]] = strings.TrimPrefix(k, "@")
		}
	}

	// ── open BAM ────────────────────────────────────────────────────────────
	bf, err := os.Open(bamPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening BAM: %v\n", err)
		os.Exit(1)
	}
	defer bf.Close()

	br, err := bam.NewReader(bufio.NewReaderSize(bf, 1<<20), 1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM reader: %v\n", err)
		os.Exit(1)
	}
	defer br.Close()

	// ── open SAM writer (stdout, no header, matching Python behaviour) ───────
	outBuf := bufio.NewWriterSize(os.Stdout, 1<<20)
	defer outBuf.Flush()

	sw, err := sam.NewWriter(outBuf, br.Header(), sam.FlagDecimal)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM writer: %v\n", err)
		os.Exit(1)
	}

	// ── process records ──────────────────────────────────────────────────────
	for {
		rec, err := br.Read()
		if err != nil {
			break // EOF or error — either way stop
		}

		// Find RX tag.
		if rx, ok := getStringTag(rec, "RX"); ok {
			bxVal, vxVal := reconstructBarcode(rx, stagger, bc)
			setStringTag(rec, "BX", bxVal)
			setIntTag(rec, "VX", vxVal)
		}

		if err := sw.Write(rec); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing record: %v\n", err)
			os.Exit(1)
		}
	}
}

// ─── SAM tag helpers ──────────────────────────────────────────────────────────
// biogo/hts stores tags as []sam.Aux. Each Aux is a raw byte slice:
// [tag0, tag1, type, ...value bytes]. We manipulate them directly to avoid
// allocating new sam.Record objects.

// getStringTag returns the string value of a two-letter tag, and whether it
// was found.
func getStringTag(r *sam.Record, tag string) (string, bool) {
	t := sam.Tag{tag[0], tag[1]}
	for _, aux := range r.AuxFields {
		if aux.Tag() == t {
			if s, ok := aux.Value().(string); ok {
				return s, true
			}
		}
	}
	return "", false
}

// setStringTag overwrites an existing tag or appends a new Z-type tag.
func setStringTag(r *sam.Record, tag, value string) {
	t := sam.Tag{tag[0], tag[1]}
	newAux, err := sam.NewAux(t, value)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating string tag %s: %v\n", tag, err)
		os.Exit(1)
	}
	for i, aux := range r.AuxFields {
		if aux.Tag() == t {
			r.AuxFields[i] = newAux
			return
		}
	}
	r.AuxFields = append(r.AuxFields, newAux)
}

// setIntTag overwrites an existing tag or appends a new i-type tag.
func setIntTag(r *sam.Record, tag string, value int) {
	t := sam.Tag{tag[0], tag[1]}
	newAux, err := sam.NewAux(t, int32(value))
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating int tag %s: %v\n", tag, err)
		os.Exit(1)
	}
	for i, aux := range r.AuxFields {
		if aux.Tag() == t {
			r.AuxFields[i] = newAux
			return
		}
	}
	r.AuxFields = append(r.AuxFields, newAux)
}


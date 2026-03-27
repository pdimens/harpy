package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"runtime"
	"strconv"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

var invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")
var stlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)
var vxTag = sam.Tag{'V', 'X'}
var bxTag = sam.Tag{'B', 'X'}

// SetBX sets a string aux tag on a record
func SetBX(rec *sam.Record, val string) {
	aux := make(sam.Aux, 3+len(val)+1)
	aux[0] = 'B'
	aux[1] = 'X'
	aux[2] = 'Z'
	copy(aux[3:], val)
	aux[len(aux)-1] = 0 // null terminator

	for i, a := range rec.AuxFields {
		if a.Tag() == bxTag {
			rec.AuxFields[i] = aux
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, aux)
}

// SetVX sets an integer (0/1) auxiliary tag on a record.
func SetVX(rec *sam.Record, val bool) {
	b := byte(1)
	if val {
		b = 0
	}
	for i, a := range rec.AuxFields {
		if a.Tag() == vxTag {
			a[3] = b
			rec.AuxFields[i] = a
			return
		}
	}
	rec.AuxFields = append(rec.AuxFields, sam.Aux{'V', 'X', 'c', b})
}
func GetStringTag(r *sam.Record, tag string) (string, bool) {
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

func copyHeaderWithPG(src *sam.Header, infile string) (*sam.Header, error) {
	// Clone the header via marshal/unmarshal
	b, err := src.MarshalText()
	if err != nil {
		return nil, err
	}
	dst := &sam.Header{}
	if err := dst.UnmarshalText(b); err != nil {
		return nil, err
	}

	// Build the @PG record
	pg := sam.NewProgram(
		"djinn",                     // ID
		"djinn",                     // name (PN)
		"djinn standardize "+infile, // command line (CL)
		lastPGID(src),               // previous PG ID (PP), or "" if none
		"1.0",                       // version (VN) — set as appropriate
	)

	if err := dst.AddProgram(pg); err != nil {
		return nil, err
	}

	return dst, nil
}

// lastPGID returns the ID of the last @PG record in the header,
// which becomes the PP (previous program) of the new entry.
// Returns "" if there are no existing @PG records.
func lastPGID(h *sam.Header) string {
	progs := h.Progs()
	if len(progs) == 0 {
		return ""
	}
	return strconv.Itoa(progs[len(progs)-1].ID())
}

func main() {
	nThreads := flag.Int("threads", 1, "Number of threads for BAM io. No real value beyond ~4-5.")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: standardize <--threads> input.bam > output.bam\n")
	}
	flag.Parse()
	args := flag.Args()
	if len(args) != 1 {
		flag.Usage()
		os.Exit(1)
	}
	infile := args[0]
	threads := min(max(*nThreads, 1), runtime.NumCPU())
	readThread := 1
	writeThread := 1

	if threads > 2 {
		readThread = 2
		writeThread = threads - readThread
	}
	// ── BAM reader ──────────────────────────────────────────────────────────────
	bf, err := os.Open(infile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening BAM: %v\n", err)
		os.Exit(1)
	}
	defer bf.Close()

	br, err := bam.NewReader(bufio.NewReaderSize(bf, 2<<20), readThread)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM reader: %v\n", err)
		os.Exit(1)
	}
	defer br.Close()

	// ── BAM writer ──────────────────────────────────────────────────────────────
	newHeader, err := copyHeaderWithPG(br.Header(), infile)
	if err != nil {
		log.Fatal(err)
	}

	// WRITES SAM OUTPUT
	//w, err := sam.NewWriter(os.Stdout, newHeader, sam.FlagDecimal)
	w, err := bam.NewWriterLevel(os.Stdout, newHeader, 4, writeThread)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	// ── loop through BAM records ───────────────────────────────────────────────
	for {
		rec, err := br.Read()
		if err != nil {
			break // EOF
		}

		bxVal, hasBX := GetStringTag(rec, "BX")
		if hasBX {
			SetVX(rec, invalid.MatchString(bxVal))
		} else {
			matches := stlfTell.FindStringSubmatch(rec.Name)
			if matches != nil {
				switch {
				case len(matches) > 1 && matches[1] != "":
					// matches[1] is the tellseq barcode e.g. "ATCGN"
					bxVal = matches[1]
				case len(matches) > 2 && matches[2] != "":
					// matches[2] is the stlfr barcode e.g. "1_2_3"
					bxVal = matches[2]
				}
			}
			SetVX(rec, invalid.MatchString(bxVal))
			SetBX(rec, bxVal)
		}

		w.Write(rec)
	}
}

// preproc_barcodes - renames and records demultiplexed linked-read barcodes
// coming out of Pheniqs. Reads a Pheniqs JSON config and a BAM file, rewrites
// BX and VX tags on every record that has an RX tag, and writes R1/R2 FASTQ
// output files compressed with pgzip.
//
// Usage:
//
//	convert [--threads N] <pheniqs.json> <input.bam> <out.R1.fq.gz> <out.R2.fq.gz>
//
// Build:
//
//	go mod init convert
//	go get github.com/biogo/hts@latest
//	go get github.com/klauspost/pgzip
//	go build -ldflags="-s -w" -o gih-convert convert.go
package main

import (
	"bufio"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/klauspost/pgzip"
)

// ─── JSON schema ──────────────────────────────────────────────────────────────

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

func reconstructBarcode(nucBC, origBC string, stagger, bc *map[string]string) (string, int, int) {
	var corrected int
	seg := strings.SplitN(nucBC, "-", 4)
	if nucBC != origBC {
		corrected = 1
	}
	if len(seg) != 4 {
		return "A00C00B00D00", 0, corrected
	}

	A := lookup(bc, &seg[1])
	B := lookup(bc, &seg[2])
	C := lookup(bc, &seg[3])
	D := lookup(stagger, &seg[0])

	valid := 1
	if *A == "00" || *B == "00" || *C == "00" || *D == "00" {
		valid = 0
	}

	return "A" + *A + "C" + *B + "B" + *C + "D" + *D, valid, corrected
}

func lookup(m *map[string]string, key *string) *string {
	if v, ok := (*m)[*key]; ok {
		return &v
	}
	empty := "00"
	return &empty
}

// ─── FASTQ writer ─────────────────────────────────────────────────────────────

// fastqWriter wraps a pgzip.Writer with a bufio layer.
// Writes go directly into the bufio buffer, which flushes to pgzip only
// when its 1MB capacity is reached — no intermediate per-record copy.
type fastqWriter struct {
	f   *os.File
	gz  *pgzip.Writer
	buf *bufio.Writer
	dir string
}

func newFastqWriter(path string, readnum string, level, threads int) (*fastqWriter, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	gz, err := pgzip.NewWriterLevel(f, level)
	if err != nil {
		f.Close()
		return nil, err
	}
	gz.SetConcurrency(2<<20, threads)
	return &fastqWriter{
		f:   f,
		gz:  gz,
		buf: bufio.NewWriterSize(gz, 2<<20),
		dir: readnum,
	}, nil
}

func (fw *fastqWriter) close() error {
	var errs []error

	if err := fw.buf.Flush(); err != nil {
		errs = append(errs, fmt.Errorf("flushing buffer: %w", err))
	}
	if err := fw.gz.Close(); err != nil {
		errs = append(errs, fmt.Errorf("closing gzip writer: %w", err))
	}
	if err := fw.f.Close(); err != nil {
		errs = append(errs, fmt.Errorf("closing file: %w", err))
	}

	return errors.Join(errs...)
}

// writeRecord writes a single FASTQ record directly into the bufio buffer.
// qual is raw Phred scores (0-40) as stored in sam.Record — +33 applied inline.
// bufio only flushes to pgzip when its 1MB buffer is full, so no intermediate
// per-record copy is needed.
func (fw *fastqWriter) writeRecord(name string, bxTag string, vxTag int, seq []byte, qual []uint8) error {
	w := errWriter{buf: fw.buf}
	w.writeString("@")
	w.writeString(name)
	w.writeString(fw.dir)
	w.writeString("\tVX:i:")
	w.writeInt(vxTag)
	w.writeString("\tBX:Z:")
	w.writeString(bxTag)
	w.writeByte('\n')
	w.write(seq)
	w.writeString("\n+\n")
	if w.err != nil {
		return w.err
	}
	for _, q := range qual {
		if err := fw.buf.WriteByte(q + 33); err != nil {
			return err
		}
	}
	return fw.buf.WriteByte('\n')
}

// errWriter wraps bufio.Writer and stops writing after the first error,
// so call sites can batch writes and check once at the end.
type errWriter struct {
	buf *bufio.Writer
	err error
}

func (w *errWriter) write(b []byte) {
	if w.err == nil {
		_, w.err = w.buf.Write(b)
	}
}

func (w *errWriter) writeString(s string) {
	if w.err == nil {
		_, w.err = w.buf.WriteString(s)
	}
}

func (w *errWriter) writeByte(b byte) {
	if w.err == nil {
		w.err = w.buf.WriteByte(b)
	}
}

func (w *errWriter) writeInt(n int) {
	if w.err == nil {
		w.err = writeInt(w.buf, n)
	}
}

// writeInt writes a non-negative integer to b without allocating.
func writeInt(b *bufio.Writer, n int) error {
	if n == 0 {
		return b.WriteByte('0')
	}
	var tmp [10]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	_, err := b.Write(tmp[i:])
	return err
}

// ─── SAM tag helpers ──────────────────────────────────────────────────────────

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

// ─── main ─────────────────────────────────────────────────────────────────────

func main() {
	nThreads := flag.Int("threads", 4, "Number of compression threads for gzip output")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: gih-convert [options] <pheniqs.json> <input.bam> <out_R1.fq.gz> <out_R2.fq.gz>\n\nOptions:\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	args := flag.Args()
	if len(args) != 4 {
		flag.Usage()
		os.Exit(1)
	}

	threads := min(max(*nThreads, 1), runtime.NumCPU())

	jsonPath := args[0]
	bamPath := args[1]
	r1Path := args[2]
	r2Path := args[3]

	// ── parse JSON ────────────────────────────────────────────────────────────
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

	// ── open BAM ──────────────────────────────────────────────────────────────
	bf, err := os.Open(bamPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening BAM: %v\n", err)
		os.Exit(1)
	}
	defer bf.Close()

	br, err := bam.NewReader(bufio.NewReaderSize(bf, 2<<20), 1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM reader: %v\n", err)
		os.Exit(1)
	}
	defer br.Close()

	// ── open FASTQ writers ────────────────────────────────────────────────────
	// Split thread budget evenly between the two writers so total
	// compression goroutines == --threads.
	threadsPerWriter := max(threads/2, 1)
	fw1, err := newFastqWriter(r1Path, "/1", 4, threadsPerWriter)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R1 output: %v\n", err)
		os.Exit(1)
	}
	defer fw1.close()

	fw2, err := newFastqWriter(r2Path, "/2", 4, threadsPerWriter)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R2 output: %v\n", err)
		os.Exit(1)
	}
	defer fw2.close()

	// ── process records ───────────────────────────────────────────────────────
	var valids int
	var corrected int
	set := make(map[string]struct{}, 1_000_000)
	for {
		rec, err := br.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintf(os.Stderr, "Error reading BAM: %v\n", err)
			os.Exit(1)
		}

		rx, hasRX := getStringTag(rec, "RX")
		ox, hasOX := getStringTag(rec, "OX")
		if !hasRX {
			continue // no RX tag — skip
		}
		if !hasOX {
			ox = rx
		}

		bxVal, vxVal, corrVal := reconstructBarcode(rx, ox, &stagger, &bc)
		valids += vxVal
		corrected += corrVal
		set[bxVal] = struct{}{}

		// Route to R1 or R2 by flag bits 0x40 / 0x80.
		var fw *fastqWriter
		if rec.Flags&sam.Read1 != 0 {
			fw = fw1
		} else {
			fw = fw2
		}

		if err := fw.writeRecord(rec.Name, bxVal, vxVal, rec.Seq.Expand(), rec.Qual); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing record: %v\n", err)
			os.Exit(1)
		}
	}

	fmt.Fprintf(os.Stdout, "Total unique barcodes:               %d\n", len(set))
	fmt.Fprintf(os.Stdout, "Total reads with valid barcodes:     %d\n", valids)
	fmt.Fprintf(os.Stdout, "Total reads with corrected barcodes: %d\n", corrected)
}

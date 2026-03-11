// stagger_gih - finds the ME sequence in R1 reads, removes it, adds the
// appropriate stagger prefix, and writes interleaved unaligned BAM to stdout.
//
// Replaces the cutadapt --info-file pipeline entirely. ME detection uses a
// bounded Hamming scan (≤2 mismatches, N in read matches anything) over the
// expected ME position window, matching cutadapt -e 0.11 --overlap 19
// --match-read-wildcards behaviour.
//
// Usage:
//
//	stagger_gih [--me <seq>] [--max-mismatch <n>] [--min-len <n>] <R1.fastq[.gz]> <R2.fastq[.gz]>
//
// Build:
//
//	go mod init stagger_gih
//	go get github.com/biogo/hts@latest
//	go build -ldflags="-s -w" -o stagger_gih stagger_gih.go
package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ── stagger pads ──────────────────────────────────────────────────────────────
// Index 0 = 7bp pad (ME found outside expected window, or not found but kept),
// index 7 = no pad (ME found exactly at position 58).

var padSeq = [8]string{
	"TTTTTTT", // plen 0 — default / out-of-window
	"CCCCCC",  // plen 1
	"GGGGG",   // plen 2
	"AAAA",    // plen 3
	"TTT",     // plen 4
	"CC",      // plen 5
	"GG",      // plen 6
	"",        // plen 7 — no pad needed
}

var padQual [8][]byte // ASCII Phred+33, filled in init()

func init() {
	for i, p := range padSeq {
		q := make([]byte, len(p))
		for j := range q {
			q[j] = 'I' // Phred 40
		}
		padQual[i] = q
	}
}

// ── FASTQ reader ──────────────────────────────────────────────────────────────

type fastqRecord struct {
	name    string
	comment string
	seq     []byte
	qual    []byte // ASCII Phred+33
}

type fastqReader struct{ r *bufio.Reader }

func newFastqReader(r io.Reader) *fastqReader {
	return &fastqReader{r: bufio.NewReaderSize(r, 1<<19)}
}

func (f *fastqReader) next() (fastqRecord, bool) {
	header, err := f.r.ReadString('\n')
	if err != nil {
		return fastqRecord{}, false
	}
	header = strings.TrimRight(header, "\r\n")
	if len(header) == 0 || header[0] != '@' {
		return fastqRecord{}, false
	}
	header = header[1:]

	var name, comment string
	if idx := strings.IndexByte(header, ' '); idx >= 0 {
		name, comment = header[:idx], header[idx+1:]
	} else {
		name = header
	}

	seqLine, _ := f.r.ReadString('\n')
	seqLine = strings.TrimRight(seqLine, "\r\n")
	_, _ = f.r.ReadString('\n') // skip '+'
	qualLine, _ := f.r.ReadString('\n')
	qualLine = strings.TrimRight(qualLine, "\r\n")

	return fastqRecord{name: name, comment: comment, seq: []byte(seqLine), qual: []byte(qualLine)}, true
}

// ── file helpers ──────────────────────────────────────────────────────────────

func openFastq(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		return struct {
			io.Reader
			io.Closer
		}{gz, multiCloser{gz, f}}, nil
	}
	return f, nil
}

type multiCloser []io.Closer

func (mc multiCloser) Close() error {
	var last error
	for _, c := range mc {
		if err := c.Close(); err != nil {
			last = err
		}
	}
	return last
}

// ── ME sequence search ────────────────────────────────────────────────────────

// findME searches for the ME sequence in read seq using a bounded Hamming scan
// over positions 44-65, covering the expected ME start range of 51-58 (after
// the barcode region + up to 7bp stagger) with a small margin either side.
// 'N' in the read matches any base and does not count as a mismatch.
// Returns the 0-based start position of the best match, or -1 if no match
// within maxMismatch is found.
func findME(seq, me []byte, maxMismatch int) int {
	meLen := len(me)
	readLen := len(seq)

	if readLen < meLen {
		return -1
	}

	// ME can start at positions 51 (max 7bp stagger) through 58 (no stagger).
	// Search a small margin either side to tolerate upstream variation.
	winStart := 44 // a few bp before the earliest possible ME start
	winEnd := 65   // a few bp past the latest possible ME start
	if winStart < 0 {
		winStart = 0
	}
	if winEnd+meLen > readLen {
		winEnd = readLen - meLen
	}
	if winStart > winEnd {
		return -1
	}

	bestPos := -1
	bestMM := maxMismatch + 1

	for pos := winStart; pos <= winEnd; pos++ {
		mm := 0
		for i := 0; i < meLen; i++ {
			r := seq[pos+i]
			m := me[i]
			if r == 'N' || m == 'N' {
				continue // wildcard — not a mismatch
			}
			if r != m {
				mm++
				if mm >= bestMM {
					break
				}
			}
		}
		if mm < bestMM {
			bestMM = mm
			bestPos = pos
			if mm == 0 {
				break // perfect match — no need to keep scanning
			}
		}
	}

	if bestMM <= maxMismatch {
		return bestPos
	}
	return -1
}

// ── BAM helpers ───────────────────────────────────────────────────────────────

func asciiQualToPhred(q []byte) []byte {
	out := make([]byte, len(q))
	for i, b := range q {
		out[i] = b - 33
	}
	return out
}

func pairedFlags(isRead1 bool) sam.Flags {
	flags := sam.Paired | sam.Unmapped | sam.MateUnmapped
	if isRead1 {
		flags |= sam.Read1
	} else {
		flags |= sam.Read2
	}
	return flags
}

func makeRecord(name string, seq, asciiQual []byte, isRead1 bool) (*sam.Record, error) {
	r, err := sam.NewRecord(
		name,
		nil, nil, // ref, mRef — unmapped
		-1, -1,   // pos, mPos
		0,        // tLen
		255,      // mapQ — unavailable
		nil,      // CIGAR
		seq,
		asciiQualToPhred(asciiQual),
		nil, // AUX
	)
	if err != nil {
		return nil, err
	}
	r.Flags = pairedFlags(isRead1)
	return r, nil
}

// ── main ──────────────────────────────────────────────────────────────────────

func main() {
	meSeq := flag.String("me", "CTGTCTCTTATACACATCT", "ME sequence to search for")
	maxMM := flag.Int("max-mismatch", 2, "Maximum mismatches allowed in ME match (default matches cutadapt -e 0.11 with 19bp ME)")
	minLen := flag.Int("min-len", 30, "Minimum biological sequence length after ME excision (excludes pad+barcode region); shorter reads are discarded")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: stagger_gih [options] <R1.fastq[.gz]> <R2.fastq[.gz]>\n\nOptions:\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	args := flag.Args()
	if len(args) != 2 {
		flag.Usage()
		os.Exit(1)
	}

	me := []byte(strings.ToUpper(*meSeq))
	meLen := len(me)

	fh1, err := openFastq(args[0])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R1: %v\n", err)
		os.Exit(1)
	}
	defer fh1.Close()

	fh2, err := openFastq(args[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R2: %v\n", err)
		os.Exit(1)
	}
	defer fh2.Close()

	rdr1 := newFastqReader(fh1)
	rdr2 := newFastqReader(fh2)

	// BAM header
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM header: %v\n", err)
		os.Exit(1)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted

	outBuf := bufio.NewWriterSize(os.Stdout, 1<<19)
	w, err := bam.NewWriter(outBuf, header, 1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM writer: %v\n", err)
		os.Exit(1)
	}

	var total, discarded, tooShort int

	for {
		fq1, ok1 := rdr1.next()
		fq2, ok2 := rdr2.next()
		if !ok1 || !ok2 {
			break
		}
		total++

		// Search for ME in R1.
		mePos := findME(fq1.seq, me, *maxMM)

		if mePos == -1 {
			// ME not found — discard both reads (equivalent to col2 == -1).
			discarded++
			continue
		}

		// Determine stagger pad length from ME start position.
		// Canonical window: positions 51-58 (for 19bp ME after 16bp barcode).
		// meStart = meLen + 32 maps to plen 0; meStart = meLen + 39 to plen 7.
		// More precisely: canonical first ME position = meLen + (meLen - meLen) = meLen.
		// The stagger occupies positions (meLen - 8) to (meLen - 1) before ME,
		// so plen = mePos - (meLen - 8)... but we match the original Python logic:
		//   col3 = mePos (0-based start of ME in original read)
		//   plen = col3 - 51  if 51 <= col3 <= 58, else 7 (max stagger / default)
		var plen int
		if mePos >= 51 && mePos <= 58 {
			plen = mePos - 51
		} else {
			plen = 7
		}

		// Build staggered R1: pad + barcode region (before ME) + sample sequence (after ME).
		// Structure: [pad][seq[:mePos]][seq[mePos+meLen:]]
		// plen==6 drops the first base of the barcode region (matching original logic).
		p := []byte(padSeq[plen])
		q := append([]byte{}, padQual[plen]...) // copy — never mutate the global

		after := fq1.seq[mePos+meLen:]
		afterQual := fq1.qual[mePos+meLen:]

		var seq, qual []byte
		if plen == 6 {
			seq = append(p, fq1.seq[1:mePos]...)
			seq = append(seq, after...)
			qual = append(q, fq1.qual[1:mePos]...)
			qual = append(qual, afterQual...)
		} else {
			seq = append(p, fq1.seq[:mePos]...)
			seq = append(seq, after...)
			qual = append(q, fq1.qual[:mePos]...)
			qual = append(qual, afterQual...)
		}

		// Discard if the biological sequence after ME is too short.
		if len(after) < *minLen {
			tooShort++
			continue
		}

		rec1, err := makeRecord(fq1.name, seq, qual, true)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error building R1 record %d: %v\n", total, err)
			os.Exit(1)
		}
		if err := w.Write(rec1); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing R1 record %d: %v\n", total, err)
			os.Exit(1)
		}

		rec2, err := makeRecord(fq2.name, fq2.seq, fq2.qual, false)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error building R2 record %d: %v\n", total, err)
			os.Exit(1)
		}
		if err := w.Write(rec2); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing R2 record %d: %v\n", total, err)
			os.Exit(1)
		}
	}

	if err := w.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "Error closing BAM writer: %v\n", err)
		os.Exit(1)
	}
	outBuf.Flush()

	fmt.Fprintf(os.Stderr, "Total read pairs processed:            %d\n", total)
	fmt.Fprintf(os.Stderr, "Discarded (ME sequence not found):     %d\n", discarded)
	fmt.Fprintf(os.Stderr, "Discarded (post-trim read too short):  %d\n", tooShort)
}

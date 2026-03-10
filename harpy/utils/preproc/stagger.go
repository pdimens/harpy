// stagger_gih - reads two FASTQ files and a cutadapt info file, adds staggers
// to barcodes, and writes interleaved unaligned BAM records to stdout.
//
// Replaces the previous FASTQ-to-stdout pipeline that required samtools import.
//
// Usage:
//
//	gih-stagger <R1.fastq[.gz]> <R2.fastq[.gz]> [info_file]
//	cat info | gih-stagger <R1.fastq[.gz]> <R2.fastq[.gz]>
//
// Build:
//
//	go mod init stagger
//	go get github.com/biogo/hts@latest
//	go build -ldflags="-s -w" -o gih-stagger stagger.go
package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// Stagger pads, mirroring the Python pad/qpad lists (index 0-7).
var pad = [8]string{
	"TTTTTTT",
	"CCCCCC",
	"GGGGG",
	"AAAA",
	"TTT",
	"CC",
	"GG",
	"",
}

var qpad [8][]byte // filled in init(); stored as ASCII Phred+33

func init() {
	for i, p := range pad {
		q := make([]byte, len(p))
		for j := range q {
			q[j] = 'I' // ASCII 73 = Phred 40 + 33
		}
		qpad[i] = q
	}
}

// fastqRecord holds one parsed FASTQ record.
type fastqRecord struct {
	name    string
	comment string
	seq     []byte
	qual    []byte // raw ASCII Phred+33
}

// fastqReader wraps a bufio.Reader and yields fastqRecords.
type fastqReader struct {
	r *bufio.Reader
}

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
		name = header[:idx]
		comment = header[idx+1:]
	} else {
		name = header
	}

	seqLine, _ := f.r.ReadString('\n')
	seqLine = strings.TrimRight(seqLine, "\r\n")

	_, _ = f.r.ReadString('\n') // skip '+'

	qualLine, _ := f.r.ReadString('\n')
	qualLine = strings.TrimRight(qualLine, "\r\n")

	return fastqRecord{
		name:    name,
		comment: comment,
		seq:     []byte(seqLine),
		qual:    []byte(qualLine),
	}, true
}

// openFastq opens a plain or gzip-compressed FASTQ file.
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

// asciiQualToPhred converts ASCII Phred+33 bytes to raw Phred scores (0-40),
// which is what biogo/hts sam.Record expects.
func asciiQualToPhred(q []byte) []byte {
	out := make([]byte, len(q))
	for i, b := range q {
		out[i] = b - 33
	}
	return out
}

// pairedFlags returns the SAM flags for a paired unmapped read:
//   - 0x1  = paired
//   - 0x4  = read unmapped
//   - 0x8  = mate unmapped
//   - 0x40 = read1 / 0x80 = read2
func pairedFlags(isRead1 bool) sam.Flags {
	flags := sam.Paired | sam.Unmapped | sam.MateUnmapped
	if isRead1 {
		flags |= sam.Read1
	} else {
		flags |= sam.Read2
	}
	return flags
}

// makeRecord builds a sam.Record for a paired unaligned read.
// sam.NewRecord does not accept Flags directly; they are set on the struct
// after construction. ref/mRef are nil for unmapped reads (*).
// mapQ is a byte (255 = unavailable).
func makeRecord(name string, seq, asciiQual []byte, isRead1 bool) (*sam.Record, error) {
	r, err := sam.NewRecord(
		name,
		nil, // ref  (*Reference) = nil -> RNAME "*"
		nil, // mRef (*Reference) = nil -> RNEXT "*"
		-1,  // p    (POS, 0-based; -1 -> "*")
		-1,  // mPos (PNEXT)
		0,   // tLen (TLEN)
		255, // mapQ (byte; 255 = unavailable)
		nil, // co   ([]CigarOp)
		seq,
		asciiQualToPhred(asciiQual),
		nil, // aux  ([]Aux)
	)
	if err != nil {
		return nil, err
	}
	r.Flags = pairedFlags(isRead1)
	return r, nil
}

func main() {
	args := os.Args[1:]
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "Usage: gih-stagger <R1.fastq[.gz]> <R2.fastq[.gz]> [info_file]\n")
		os.Exit(1)
	}

	fq1Path := args[0]
	fq2Path := args[1]

	var infoSrc io.Reader = os.Stdin
	if len(args) >= 3 {
		f, err := os.Open(args[2])
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening info file: %v\n", err)
			os.Exit(1)
		}
		defer f.Close()
		infoSrc = f
	}

	fh1, err := openFastq(fq1Path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R1: %v\n", err)
		os.Exit(1)
	}
	defer fh1.Close()

	fh2, err := openFastq(fq2Path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R2: %v\n", err)
		os.Exit(1)
	}
	defer fh2.Close()

	rdr1 := newFastqReader(fh1)
	rdr2 := newFastqReader(fh2)
	infoScanner := bufio.NewScanner(bufio.NewReaderSize(infoSrc, 1<<19))

	// Build a minimal uBAM header.
	// sam.Header exposes Version, SortOrder, GroupOrder directly — no HD sub-struct.
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM header: %v\n", err)
		os.Exit(1)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted

	// Stream BAM to stdout. gzip level 1 (fast) is appropriate for uBAM —
	// the downstream aligner will decompress and re-sort anyway.
	outBuf := bufio.NewWriterSize(os.Stdout, 1<<19)
	w, err := bam.NewWriter(outBuf, header, 1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM writer: %v\n", err)
		os.Exit(1)
	}

	var total, discarded int

	for infoScanner.Scan() {
		line := infoScanner.Text()

		fq1, ok1 := rdr1.next()
		fq2, ok2 := rdr2.next()
		if !ok1 || !ok2 {
			break
		}
		total++

		fields := strings.SplitN(line, "\t", 4)
		if len(fields) < 3 {
			fmt.Fprintf(os.Stderr, "Malformed info line %d: %q\n", total, line)
			os.Exit(1)
		}

		caName, _, _ := strings.Cut(fields[0], " ")

		col2, err := strconv.Atoi(fields[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "Bad col2 at record %d: %v\n", total, err)
			os.Exit(1)
		}

		if col2 == -1 {
			discarded++
			continue
		}

		if caName != fq1.name {
			fmt.Fprintf(os.Stderr,
				"Error: Read name mismatch at record %d. Expected:\n%s\nbut found:\n%s\n",
				total, caName, fq1.name)
			os.Exit(1)
		}

		col3, err := strconv.Atoi(fields[2])
		if err != nil {
			fmt.Fprintf(os.Stderr, "Bad col3 at record %d: %v\n", total, err)
			os.Exit(1)
		}

		var plen int
		if col3 < 51 || col3 > 58 {
			plen = 7
		} else {
			plen = col3 - 51
		}

		p := []byte(pad[plen])
		q := append([]byte{}, qpad[plen]...)

		var seq, qual []byte
		if plen == 6 {
			seq = append(p, fq1.seq[1:]...)
			qual = append(q, fq1.qual[1:]...)
		} else {
			seq = append(p, fq1.seq...)
			qual = append(q, fq1.qual...)
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

	if err := infoScanner.Err(); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading info: %v\n", err)
		os.Exit(1)
	}

	if err := w.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "Error closing BAM writer: %v\n", err)
		os.Exit(1)
	}
	outBuf.Flush()

	fmt.Fprintf(os.Stderr, "Reads discarded because of missing ME sequence: %d\n", discarded)
}

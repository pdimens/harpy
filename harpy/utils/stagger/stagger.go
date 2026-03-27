// gih-stagger - finds the ME sequence in R1 reads, removes it, adds the
// appropriate stagger prefix, and writes interleaved unaligned SAM to stdout.
//
// Replaces the cutadapt --info-file pipeline entirely. ME detection uses a
// bounded Hamming scan (≤2 mismatches, N in read matches anything) over the
// expected ME position window, matching cutadapt -e 0.11 --overlap 19
// --match-read-wildcards behavior.
//
// Usage:
//
//	gih-stagger [options] <R1.fastq[.gz]> <R2.fastq[.gz]>
//
// Build:
//
//	go mod init stagger
//	go get github.com/biogo/hts@latest
//	go get github.com/klauspost/pgzip
//	go build -ldflags="-s -w" -o gih-stagger stagger.go
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/sam"
	"github.com/klauspost/pgzip"
)

const R1flag int = 77
const R2flag int = 141

// ── stagger pads ──────────────────────────────────────────────────────────────

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

// padEntry stores each pad as a fixed [7]byte array with Phred scores
// pre-converted, alongside its true length n. Workers copy with no allocation.
type padEntry struct {
	seq  [7]byte
	qual [7]byte // raw Phred scores (40) — not ASCII
	n    int
}

var pads [8]padEntry

func init() {
	for i, s := range padSeq {
		var e padEntry
		e.n = len(s)
		for j, c := range []byte(s) {
			e.seq[j] = c
			e.qual[j] = 40 // Phred 40
		}
		pads[i] = e
	}
}

const batchSize = 2000 // read pairs per batch

// ── FASTQ reader ──────────────────────────────────────────────────────────────

type fastqRecord struct {
	name    string
	comment string
	seq     []byte
	qual    []byte // ASCII Phred+33
}

type readPair struct {
	fq1, fq2 fastqRecord
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
	if before, after, ok := strings.Cut(header, " "); ok {
		name, comment = before, after
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

func openFastq(path string, gzBlocks int) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(path, ".gz") {
		gz, err := pgzip.NewReaderN(f, 3<<20, gzBlocks)
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

// findME searches for the ME sequence using a bounded Hamming scan over
// positions 44-65, covering the expected ME start range of 51-58 with margin.
// 'N' in the read matches any base and does not count as a mismatch.
// Returns the 0-based start position of the best match, or -1 if not found.
func findME(seq, me *[]byte, maxMismatch *int) int {
	meLen := len(*me)
	readLen := len(*seq)
	if readLen < meLen {
		return -1
	}
	winStart := 44
	winEnd := 65
	if winEnd+meLen > readLen {
		winEnd = readLen - meLen
	}
	if winStart > winEnd {
		return -1
	}

	bestPos := -1
	bestMM := *maxMismatch + 1

	for pos := winStart; pos <= winEnd; pos++ {
		mm := 0
		for i := range meLen {
			r := (*seq)[pos+i]
			if r == 'N' {
				continue
			}
			if r != (*me)[i] {
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
				break
			}
		}
	}

	if bestMM <= *maxMismatch {
		return bestPos
	}
	return -1
}

// ── SAM formatting ────────────────────────────────────────────────────────────

// pairedFlag returns the SAM FLAG integer for a paired unmapped read.
// R1: 0x1|0x4|0x8|0x40 = 77   R2: 0x1|0x4|0x8|0x80 = 141
func pairedFlag(isRead1 bool) int {
	if isRead1 {
		return R1flag
	}
	return R2flag
}

// workerState holds reusable per-goroutine scratch buffers.
type workerState struct {
	seqBuf  []byte       // assembled seq (pad + barcode + biological)
	qualBuf []byte       // assembled qual in ASCII Phred+33
	outBuf  bytes.Buffer // SAM text accumulator for the whole batch
}

// writeSAMRecord appends a single SAM record as a text line to ws.outBuf.
// qual is expected in ASCII Phred+33 (as read from FASTQ).
// This avoids constructing a sam.Record object entirely for the write path.
func (ws *workerState) writeSAMRecord(name string, flag int, seq, qual []byte) {
	ws.outBuf.WriteString(name)
	ws.outBuf.WriteByte('\t')
	// write flag as decimal
	ws.writeInt(flag)
	// RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN — all * or 0 for uSAM
	ws.outBuf.WriteString("\t*\t0\t255\t*\t*\t0\t0\t")
	ws.outBuf.Write(seq)
	ws.outBuf.WriteByte('\t')
	ws.outBuf.Write(qual)
	ws.outBuf.WriteByte('\n')
}

// writeInt writes a non-negative integer to outBuf without allocating.
func (ws *workerState) writeInt(n int) {
	if n == 0 {
		ws.outBuf.WriteByte('0')
		return
	}
	var tmp [10]byte
	i := len(tmp)
	for n > 0 {
		i--
		tmp[i] = byte('0' + n%10)
		n /= 10
	}
	ws.outBuf.Write(tmp[i:])
}

// ── worker ────────────────────────────────────────────────────────────────────

func processBatch(
	ws *workerState,
	batch []readPair,
	me []byte,
	maxMismatch, minLen int,
	out io.Writer,
	mu *sync.Mutex,
	discarded, tooShort *atomic.Int64,
) error {
	meLen := len(me)
	ws.outBuf.Reset()

	for i := range batch {
		fq1 := &batch[i].fq1
		fq2 := &batch[i].fq2

		mePos := findME(&fq1.seq, &me, &maxMismatch)
		if mePos == -1 {
			discarded.Add(1)
			continue
		}

		var plen int
		if mePos >= 51 && mePos <= 58 {
			plen = mePos - 51
		} else {
			plen = 7
		}

		pad := &pads[plen]
		after := fq1.seq[mePos+meLen:]
		afterQual := fq1.qual[mePos+meLen:]

		if len(after) < minLen {
			tooShort.Add(1)
			continue
		}

		var bcSeq, bcQual []byte
		if plen == 6 {
			bcSeq = fq1.seq[1:mePos]
			bcQual = fq1.qual[1:mePos]
		} else {
			bcSeq = fq1.seq[:mePos]
			bcQual = fq1.qual[:mePos]
		}

		// Assemble seq: pad bases + barcode + biological sequence
		totalLen := pad.n + len(bcSeq) + len(after)
		if cap(ws.seqBuf) < totalLen {
			ws.seqBuf = make([]byte, totalLen)
			ws.qualBuf = make([]byte, totalLen)
		}
		ws.seqBuf = ws.seqBuf[:totalLen]
		ws.qualBuf = ws.qualBuf[:totalLen]

		n := copy(ws.seqBuf, pad.seq[:pad.n])
		n += copy(ws.seqBuf[n:], bcSeq)
		copy(ws.seqBuf[n:], after)

		// Assemble qual in ASCII Phred+33:
		// pad qual is stored as raw Phred (40), so add 33 back for SAM text output.
		for j := range pad.n {
			ws.qualBuf[j] = pad.qual[j] + 33 // 40+33 = 73 = 'I'
		}
		n = copy(ws.qualBuf[pad.n:], bcQual)
		copy(ws.qualBuf[pad.n+n:], afterQual)

		ws.writeSAMRecord(fq1.name, pairedFlag(true), ws.seqBuf, ws.qualBuf)
		ws.writeSAMRecord(fq2.name, pairedFlag(false), fq2.seq, fq2.qual)
	}

	// Single locked write for the entire batch.
	mu.Lock()
	_, err := out.Write(ws.outBuf.Bytes())
	mu.Unlock()
	return err
}

// ── main ──────────────────────────────────────────────────────────────────────

func main() {
	meSeq := flag.String("me", "CTGTCTCTTATACACATCT", "ME sequence to search for")
	statsfile := flag.String("stats", "stats.txt", "File name for stats output")
	maxMM := flag.Int("max-mismatch", 2, "Maximum mismatches allowed in ME match (matches cutadapt -e 0.11 with 19bp ME)")
	minLen := flag.Int("min-len", 30, "Minimum biological sequence length after ME excision; shorter reads are discarded")
	nThreads := flag.Int("threads", 4, "Number of worker threads")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: gih-stagger [options] <R1.fastq[.gz]> <R2.fastq[.gz]>\n\nOptions:\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	args := flag.Args()
	if len(args) != 2 {
		flag.Usage()
		os.Exit(1)
	}

	threads := min(max(*nThreads, 1), runtime.NumCPU())
	runtime.GOMAXPROCS(threads + 1) // +1 for the reader goroutine

	// Split decompression block budget evenly between the two input readers.
	gzBlocks := max(threads/2, 1)

	me := []byte(strings.ToUpper(*meSeq))

	fh1, err := openFastq(args[0], gzBlocks)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R1: %v\n", err)
		os.Exit(1)
	}
	defer fh1.Close()

	fh2, err := openFastq(args[1], gzBlocks)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening R2: %v\n", err)
		os.Exit(1)
	}
	defer fh2.Close()

	rdr1 := newFastqReader(fh1)
	rdr2 := newFastqReader(fh2)

	// Write SAM header then stream records.
	// sam.NewHeader + sam.NewWriter handles the @HD line; we write no @SQ lines
	// since all reads are unmapped.
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM header: %v\n", err)
		os.Exit(1)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted

	outBuf := bufio.NewWriterSize(os.Stdout, 2<<20)

	// Write header via sam.Writer, then use outBuf directly for record lines.
	headerWriter, err := sam.NewWriter(outBuf, header, sam.FlagDecimal)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM writer: %v\n", err)
		os.Exit(1)
	}
	// sam.NewWriter writes the header on construction; we only needed it for that.
	_ = headerWriter

	var (
		mu        sync.Mutex
		wg        sync.WaitGroup
		discarded atomic.Int64
		tooShort  atomic.Int64
		total     atomic.Int64
		writeErr  atomic.Value
	)

	jobs := make(chan []readPair, threads*2)

	for range threads {
		wg.Go(func() {
			ws := &workerState{
				seqBuf:  make([]byte, 512),
				qualBuf: make([]byte, 512),
			}
			ws.outBuf.Grow(batchSize * 300) // pre-size for ~300 bytes/record pair
			for batch := range jobs {
				if writeErr.Load() != nil {
					continue
				}
				if err := processBatch(ws, batch, me, *maxMM, *minLen, outBuf, &mu, &discarded, &tooShort); err != nil {
					writeErr.Store(err)
				}
			}
		})
	}

	batch := make([]readPair, 0, batchSize)
	for {
		fq1, ok1 := rdr1.next()
		fq2, ok2 := rdr2.next()
		if !ok1 || !ok2 {
			break
		}
		total.Add(1)
		batch = append(batch, readPair{fq1, fq2})
		if len(batch) == batchSize {
			jobs <- batch
			batch = make([]readPair, 0, batchSize)
		}
	}
	if len(batch) > 0 {
		jobs <- batch
	}
	close(jobs)
	wg.Wait()

	if err, ok := writeErr.Load().(error); ok && err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	outBuf.Flush()

	file, err := os.Create(*statsfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating stats file: %v\n", err)
		os.Exit(1)
	}
	defer file.Close()

	fmt.Fprintf(file, "Total read pairs processed:           %d\n", total.Load())
	fmt.Fprintf(file, "Discarded (ME sequence not found):    %d\n", discarded.Load())
	fmt.Fprintf(file, "Discarded (post-trim read too short): %d\n", tooShort.Load())
}

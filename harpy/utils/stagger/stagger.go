package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/klauspost/pgzip"
)

const R1flag = sam.Paired | sam.Unmapped | sam.MateUnmapped | sam.Read1
const R2flag = sam.Paired | sam.Unmapped | sam.MateUnmapped | sam.Read2
const winStart = 49
const batchSize = 2000

// ── stagger pads ──────────────────────────────────────────────────────────────

var padSeq = [8]string{
	"TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", "",
}

// padEntry stores raw Phred qual (not ASCII) since BAM stores raw scores.
type padEntry struct {
	n    int
	seq  [7]byte
	qual [7]byte
}

var pads [8]padEntry

func init() {
	for i, s := range padSeq {
		var e padEntry
		e.n = len(s)
		for j, c := range []byte(s) {
			e.seq[j] = c
			e.qual[j] = 40 // raw Phred 40, no +33 offset — BAM wants raw scores
		}
		pads[i] = e
	}
}

// ── FASTQ reader ──────────────────────────────────────────────────────────────
//
// Each record owns one combined seq+qual allocation (instead of two), and
// header/plus lines are read via ReadSlice (no throwaway string alloc).
// Qual is converted to raw Phred (ASCII - 33) in place while copying, since
// BAM needs raw scores and nothing downstream needs the ASCII form again.

type fastqRecord struct {
	name string
	seq  []byte
	qual []byte // raw Phred, NOT ASCII+33
}

type readPair struct {
	fq1, fq2 fastqRecord
}

type fastqReader struct{ r *bufio.Reader }

func newFastqReader(r io.Reader) *fastqReader {
	return &fastqReader{r: bufio.NewReaderSize(r, 1<<20)}
}

// next reads one record. Returns io.EOF on clean end-of-stream (empty read
// at a record boundary). Any other error means a truncated/malformed record.
func (f *fastqReader) next() (fastqRecord, error) {
	header, err := f.r.ReadSlice('\n')
	if err != nil {
		if err == io.EOF && len(header) == 0 {
			return fastqRecord{}, io.EOF
		}
		return fastqRecord{}, fmt.Errorf("truncated header line: %w", err)
	}
	if len(header) == 0 || header[0] != '@' {
		return fastqRecord{}, fmt.Errorf("malformed header, expected '@', got %q", header)
	}
	header = bytes.TrimRight(header[1:], "\r\n")
	name := header
	if i := bytes.IndexByte(header, ' '); i >= 0 {
		name = header[:i]
	}
	nameCopy := string(name)

	seqLine, err := f.r.ReadSlice('\n')
	if err != nil {
		return fastqRecord{}, fmt.Errorf("truncated seq line for %q: %w", nameCopy, err)
	}
	seqLine = bytes.TrimRight(seqLine, "\r\n")

	if _, err := f.r.ReadSlice('\n'); err != nil {
		return fastqRecord{}, fmt.Errorf("truncated '+' line for %q: %w", nameCopy, err)
	}

	qualLine, err := f.r.ReadSlice('\n')
	if err != nil {
		return fastqRecord{}, fmt.Errorf("truncated qual line for %q: %w", nameCopy, err)
	}
	qualLine = bytes.TrimRight(qualLine, "\r\n")

	if len(seqLine) != len(qualLine) {
		return fastqRecord{}, fmt.Errorf("seq/qual length mismatch for %q: %d vs %d", nameCopy, len(seqLine), len(qualLine))
	}

	buf := make([]byte, len(seqLine)+len(qualLine))
	seq := buf[:len(seqLine)]
	qual := buf[len(seqLine):]
	copy(seq, seqLine)
	for i, q := range qualLine {
		qual[i] = q - 33
	}

	return fastqRecord{name: nameCopy, seq: seq, qual: qual}, nil
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

func findME(seq, me *[]byte, maxMismatch *float64) int {
	meLen := len(*me)
	readLen := len(*seq)
	winEnd := min(65, readLen-meLen)
	if readLen < meLen || winStart > winEnd {
		return -1
	}
seqscan:
	for pos := winStart; pos <= winEnd; pos++ {
		misMatch := 0.0
		for i := range meLen {
			if r := (*seq)[pos+i]; r == 'N' {
				misMatch += 0.3
			} else if r != (*me)[i] {
				misMatch++
			}
			if misMatch > *maxMismatch {
				continue seqscan
			}
		}
		return pos
	}
	return -1
}

// ── worker ────────────────────────────────────────────────────────────────────
//
// Builds sam.Record structs directly (no SAM text) and hands them to a
// shared bam.Writer under one lock per batch. Records are reused across
// iterations via workerState — only Seq re-packing allocates (sam.NewSeq
// packs into 4-bit nybbles; unavoidable, it's the BAM wire format).

type workerState struct {
	seqBuf, qualBuf []byte
	rec1, rec2      sam.Record
}

func newWorkerState() *workerState {
	return &workerState{
		seqBuf:  make([]byte, 512),
		qualBuf: make([]byte, 512),
		rec1:    sam.Record{Ref: nil, Pos: -1, MapQ: 0, MateRef: nil, MatePos: -1},
		rec2:    sam.Record{Ref: nil, Pos: -1, MapQ: 0, MateRef: nil, MatePos: -1},
	}
}

func processBatch(
	ws *workerState,
	batch []readPair,
	me []byte,
	maxMismatch float64,
	bw *bam.Writer,
	mu *sync.Mutex,
	discarded *atomic.Int64,
) error {
	meLen := len(me)

	mu.Lock()
	defer mu.Unlock()

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
		afterStart := mePos + meLen
		afterSeq := fq1.seq[afterStart:]
		afterQual := fq1.qual[afterStart:] // already raw Phred

		var bcSeq, bcQual []byte
		if plen == 6 {
			bcSeq = fq1.seq[1:mePos]
			bcQual = fq1.qual[1:mePos]
		} else {
			bcSeq = fq1.seq[:mePos]
			bcQual = fq1.qual[:mePos]
		}

		totalLen := pad.n + len(bcSeq) + len(afterSeq)
		if cap(ws.seqBuf) < totalLen {
			ws.seqBuf = make([]byte, totalLen)
			ws.qualBuf = make([]byte, totalLen)
		}
		ws.seqBuf = ws.seqBuf[:totalLen]
		ws.qualBuf = ws.qualBuf[:totalLen]

		n := copy(ws.seqBuf, pad.seq[:pad.n])
		n += copy(ws.seqBuf[n:], bcSeq)
		copy(ws.seqBuf[n:], afterSeq)

		copy(ws.qualBuf, pad.qual[:pad.n]) // already raw Phred
		n = copy(ws.qualBuf[pad.n:], bcQual)
		copy(ws.qualBuf[pad.n+n:], afterQual)

		seq1 := sam.NewSeq(ws.seqBuf)
		ws.rec1.Name = fq1.name
		ws.rec1.Flags = R1flag
		ws.rec1.Seq = seq1
		ws.rec1.Qual = append(ws.rec1.Qual[:0], ws.qualBuf...)
		if err := bw.Write(&ws.rec1); err != nil {
			return err
		}

		seq2 := sam.NewSeq(fq2.seq)
		ws.rec2.Name = fq2.name
		ws.rec2.Flags = R2flag
		ws.rec2.Seq = seq2
		ws.rec2.Qual = append(ws.rec2.Qual[:0], fq2.qual...) // already raw Phred
		if err := bw.Write(&ws.rec2); err != nil {
			return err
		}
	}
	return nil
}

// ── main ──────────────────────────────────────────────────────────────────────

func main() {
	meSeq := flag.String("me", "AGATGTGTATAAGAGACAG", "ME sequence to search for")
	statsfile := flag.String("stats", "stats.txt", "File name for stats output")
	maxMM := flag.Float64("max-mismatch", 2, "Maximum allowable ME mismatch score")
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
	runtime.GOMAXPROCS(threads + 1)

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

	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating SAM header: %v\n", err)
		os.Exit(1)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted

	stdoutBuf := bufio.NewWriterSize(os.Stdout, 1<<20)
	// gzip.NoCompression -> BGZF container still valid BAM, but store blocks
	// instead of deflating. wc = compress worker count (irrelevant at level 0
	// but bam.NewWriterLevel still wants it).
	bw, err := bam.NewWriterLevel(stdoutBuf, header, gzip.NoCompression, threads)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating BAM writer: %v\n", err)
		os.Exit(1)
	}

	var (
		mu        sync.Mutex
		wg        sync.WaitGroup
		discarded atomic.Int64
		total     atomic.Int64
		writeErr  atomic.Value
	)
	jobs := make(chan []readPair, threads*2)

	for range threads {
		wg.Go(func() {
			ws := newWorkerState()
			for batch := range jobs {
				if writeErr.Load() != nil {
					continue
				}
				if err := processBatch(ws, batch, me, *maxMM, bw, &mu, &discarded); err != nil {
					writeErr.Store(err)
				}
			}
		})
	}

	batch := make([]readPair, 0, batchSize)
	for {
		fq1, err1 := rdr1.next()
		fq2, err2 := rdr2.next()

		eof1 := errors.Is(err1, io.EOF)
		eof2 := errors.Is(err2, io.EOF)

		if eof1 && eof2 {
			break // clean end on both — expected
		}
		if err1 != nil && !eof1 {
			fmt.Fprintf(os.Stderr, "Error: R1 malformed near pair %d: %v\n", total.Load()+1, err1)
			os.Exit(1)
		}
		if err2 != nil && !eof2 {
			fmt.Fprintf(os.Stderr, "Error: R2 malformed near pair %d: %v\n", total.Load()+1, err2)
			os.Exit(1)
		}
		if eof1 != eof2 {
			fmt.Fprintf(os.Stderr, "Error: R1/R2 read count mismatch near pair %d (R1 EOF=%v, R2 EOF=%v)\n", total.Load()+1, eof1, eof2)
			os.Exit(1)
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

	if err := bw.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "Error closing BAM writer: %v\n", err)
		os.Exit(1)
	}
	stdoutBuf.Flush()

	file, err := os.Create(*statsfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating stats file: %v\n", err)
		os.Exit(1)
	}
	defer file.Close()

	fmt.Fprintf(file, "Total read pairs processed:           %d\n", total.Load())
	fmt.Fprintf(file, "Discarded (ME sequence not found):    %d\n", discarded.Load())
}

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
//	stagger_gih [options] <R1.fastq[.gz]> <R2.fastq[.gz]>
//
// Build:
//
//	go mod init stagger_gih
//	go get github.com/biogo/hts@latest
//	go get github.com/klauspost/pgzip
//	go build -ldflags="-s -w" -o stagger_gih stagger_gih.go
package main

import (
	"bufio"
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

// padSeqBytes and padQualBytes are pre-converted byte slices of padSeq/padQual.
// Each entry is stored as a fixed-size array of the pad's max possible length (7)
// alongside its true length, so workers can copy with no allocation.
type padEntry struct {
	seq  [7]byte
	qual [7]byte // Phred scores (not ASCII) — subtracted 33 from 'I'=73 → 40
	n    int     // actual pad length
}

var pads [8]padEntry

func init() {
	for i, s := range padSeq {
		var e padEntry
		e.n = len(s)
		for j, c := range []byte(s) {
			e.seq[j] = c
			e.qual[j] = 40 // Phred 40 directly — no ASCII subtract needed at write time
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

func openFastq(path string) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(path, ".gz") {
		gz, err := pgzip.NewReader(f) // parallel gzip decompression
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

	winStart := 44
	winEnd := 65
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
		for i := range meLen {
			r := seq[pos+i]
			m := me[i]
			if r == 'N' {
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
				break // perfect match — stop scanning
			}
		}
	}

	if bestMM <= maxMismatch {
		return bestPos
	}
	return -1
}

// ── BAM helpers ───────────────────────────────────────────────────────────────

func pairedFlags(isRead1 bool) sam.Flags {
	flags := sam.Paired | sam.Unmapped | sam.MateUnmapped
	if isRead1 {
		flags |= sam.Read1
	} else {
		flags |= sam.Read2
	}
	return flags
}

// workerState holds reusable per-goroutine scratch buffers, eliminating
// per-record heap allocations for qual conversion and sequence assembly.
type workerState struct {
	phredBuf []byte // reusable Phred conversion scratch buffer
	seqBuf   []byte // reusable assembled seq buffer
	qualBuf  []byte // reusable assembled qual buffer
}

// asciiQualToPhredInto converts ASCII Phred+33 bytes into raw Phred scores
// (0-40) using the worker's reusable scratch buffer.
func (ws *workerState) asciiQualToPhredInto(q []byte) []byte {
	if cap(ws.phredBuf) < len(q) {
		ws.phredBuf = make([]byte, len(q))
	}
	ws.phredBuf = ws.phredBuf[:len(q)]
	for i, b := range q {
		ws.phredBuf[i] = b - 33
	}
	return ws.phredBuf
}

// recordPool recycles sam.Record objects to reduce GC pressure.
var recordPool = sync.Pool{
	New: func() any { return &sam.Record{} },
}

// makeRecord builds a sam.Record using a pooled object where possible.
// seq and asciiQual are the assembled (already-padded) sequences.
// The qual conversion uses the worker's scratch buffer.
func makeRecord(ws *workerState, name string, seq, asciiQual []byte, isRead1 bool) (*sam.Record, error) {
	phred := ws.asciiQualToPhredInto(asciiQual)
	r, err := sam.NewRecord(
		name,
		nil, nil,
		-1, -1,
		0,
		255,
		nil,
		seq,
		phred,
		nil,
	)
	if err != nil {
		return nil, err
	}
	r.Flags = pairedFlags(isRead1)
	return r, nil
}

// ── worker ────────────────────────────────────────────────────────────────────

// processBatch processes a batch of read pairs: finds ME, builds staggered
// sequences using worker-local buffers, and writes R1+R2 BAM records
// atomically under mu.
func processBatch(
	ws *workerState,
	batch []readPair,
	me []byte,
	maxMismatch, minLen int,
	w *bam.Writer,
	mu *sync.Mutex,
	discarded, tooShort *atomic.Int64,
) error {
	meLen := len(me)

	for i := range batch {
		fq1 := &batch[i].fq1
		fq2 := &batch[i].fq2

		mePos := findME(fq1.seq, me, maxMismatch)
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

		// Assemble seq into worker-local buffer:
		// [pad.seq[:pad.n]] + [barcode region] + [biological seq after ME]
		var bcSeq, bcQual []byte
		if plen == 6 {
			bcSeq = fq1.seq[1:mePos]
			bcQual = fq1.qual[1:mePos]
		} else {
			bcSeq = fq1.seq[:mePos]
			bcQual = fq1.qual[:mePos]
		}

		totalLen := pad.n + len(bcSeq) + len(after)
		if cap(ws.seqBuf) < totalLen {
			ws.seqBuf = make([]byte, totalLen)
			ws.qualBuf = make([]byte, totalLen)
		}
		ws.seqBuf = ws.seqBuf[:totalLen]
		ws.qualBuf = ws.qualBuf[:totalLen]

		// seq: pad bytes (raw bases) + barcode + after
		n := copy(ws.seqBuf, pad.seq[:pad.n])
		n += copy(ws.seqBuf[n:], bcSeq)
		copy(ws.seqBuf[n:], after)

		// qual: pad bytes (Phred scores, pre-converted) + barcode qual (ASCII→Phred) + after qual (ASCII→Phred)
		n = copy(ws.qualBuf, pad.qual[:pad.n])
		for j, b := range bcQual {
			ws.qualBuf[n+j] = b - 33
		}
		n += len(bcQual)
		for j, b := range afterQual {
			ws.qualBuf[n+j] = b - 33
		}

		// For R1 we pass already-Phred qual directly — makeRecord must not
		// subtract 33 again. We temporarily reuse ws.phredBuf to pass the
		// already-converted qual through the sam.NewRecord path by setting
		// phredBuf to point at qualBuf (no copy needed).
		ws.phredBuf = ws.qualBuf

		rec1, err := sam.NewRecord(fq1.name, nil, nil, -1, -1, 0, 255, nil, ws.seqBuf, ws.phredBuf, nil)
		if err != nil {
			return fmt.Errorf("building R1 record: %w", err)
		}
		rec1.Flags = pairedFlags(true)

		// R2: convert qual with worker scratch buffer
		rec2, err := makeRecord(ws, fq2.name, fq2.seq, fq2.qual, false)
		if err != nil {
			return fmt.Errorf("building R2 record: %w", err)
		}

		mu.Lock()
		err1 := w.Write(rec1)
		err2 := w.Write(rec2)
		mu.Unlock()

		if err1 != nil {
			return fmt.Errorf("writing R1: %w", err1)
		}
		if err2 != nil {
			return fmt.Errorf("writing R2: %w", err2)
		}
	}
	return nil
}

// ── main ──────────────────────────────────────────────────────────────────────

func main() {
	meSeq := flag.String("me", "CTGTCTCTTATACACATCT", "ME sequence to search for")
	maxMM := flag.Int("max-mismatch", 2, "Maximum mismatches allowed in ME match (matches cutadapt -e 0.11 with 19bp ME)")
	minLen := flag.Int("min-len", 30, "Minimum biological sequence length after ME excision; shorter reads are discarded")
	nThreads := flag.Int("threads", 4, "Number of worker threads")
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

	threads := min(max(*nThreads, 1), runtime.NumCPU())
	runtime.GOMAXPROCS(threads + 1) // +1 for the reader goroutine

	me := []byte(strings.ToUpper(*meSeq))

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

	var (
		mu        sync.Mutex
		wg        sync.WaitGroup
		discarded atomic.Int64
		tooShort  atomic.Int64
		total     atomic.Int64
		writeErr  atomic.Value
	)

	// Buffered jobs channel — enough depth to keep all workers fed.
	jobs := make(chan []readPair, threads*2)

	// Launch worker pool — each worker owns its own workerState.
	for range threads {
		wg.Go(func() {
			ws := &workerState{
				phredBuf: make([]byte, 256),
				seqBuf:   make([]byte, 512),
				qualBuf:  make([]byte, 512),
			}
			for batch := range jobs {
				if writeErr.Load() != nil {
					continue // drain channel after an error
				}
				if err := processBatch(ws, batch, me, *maxMM, *minLen, w, &mu, &discarded, &tooShort); err != nil {
					writeErr.Store(err)
				}
			}
		})
	}

	// Read loop — runs on main goroutine, fills batches and sends to workers.
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

	if err := w.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "Error closing BAM writer: %v\n", err)
		os.Exit(1)
	}
	outBuf.Flush()

	fmt.Fprintf(os.Stderr, "Total read pairs processed:            %d\n", total.Load())
	fmt.Fprintf(os.Stderr, "Discarded (ME sequence not found):     %d\n", discarded.Load())
	fmt.Fprintf(os.Stderr, "Discarded (post-trim read too short):  %d\n", tooShort.Load())
}

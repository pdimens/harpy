// gih-demux - FASTQ in, FASTQ out replacement for the gih-stagger | pheniqs mux | gih-convert chain.
//
// Per read pair it (1) finds the ME in R1 and applies the stagger padding (as gih-stagger),
// (2) decodes the D/A/B/C barcode segments with a port of the Pheniqs PAMLD decoder
// (pamld.cpp, barcode.h, phred.cpp), driven by the same pheniqs.config.json, and
// (3) writes R1/R2 FASTQ with BX/VX tags (as gih-convert). No BAM is ever made.
//
// Usage:
//
//	demux [options] <pheniqs.json> <R1.fq[.gz]> <R2.fq[.gz]> <out.R1.fq.gz> <out.R2.fq.gz>
//
// Build:
//
//	go build -ldflags="-s -w" -o gih-demux demux.go
package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/klauspost/pgzip"
)

const (
	winStart  = 49
	batchSize = 2000
	empty     = "00"
	asciiQ40  = 'I' // Phred 40 + 33, quality of stagger pad bases
)

// ── stagger pads ──────────────────────────────────────────────────────────────

var padSeq = [8]string{"TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", ""}

// ── PAMLD decoder (port of pheniqs pamld.cpp) ─────────────────────────────────

// score[q][e][o]: -10log10 P(observe o | expect e, Phred q); o==4 is N. Mirrors
// PhredScale::assemble_substitution_lookup, including q==0 scoring 0 (table never filled).
var score [128][4][5]float64

var baseCode [256]uint8 // A,C,G,T -> 0..3, everything else -> 4 (N)

func init() {
	for i := range baseCode {
		baseCode[i] = 4
	}
	for i, c := range "ACGT" {
		baseCode[c] = uint8(i)
		baseCode[c+32] = uint8(i)
	}
	for q := 1; q < 128; q++ {
		perr := pow10(-0.1 * float64(q))
		for e := 0; e < 4; e++ {
			for o := 0; o < 5; o++ {
				switch {
				case o == 4:
					score[q][e][o] = 10 * math.Log10(4)
				case o == e:
					score[q][e][o] = -10 * math.Log10(1-perr)
				default:
					score[q][e][o] = float64(q)
				}
			}
		}
	}
}

var ln10 = math.Log(10)

// pow10 returns math.Pow(10, y) bit for bit for finite y <= 0, minus the
// per-call Log(10) and the special-case ladder; same Modf/Exp/Frexp/Ldexp steps
// as math.pow. demux_test.go checks it against math.Pow.
func pow10(y float64) float64 {
	switch y {
	case 0:
		return 1
	case -0.5:
		return 1 / math.Sqrt(10)
	}
	yi, yf := math.Modf(-y)
	if yf > 0.5 {
		yf--
		yi++
	}
	a1 := 1.0
	if yf != 0 {
		a1 = math.Exp(yf * ln10)
	}
	x1, xe := 0.625, 4 // Frexp(10)
	ae := 0
	for i := int64(yi); i != 0; i >>= 1 {
		if i&1 == 1 {
			a1 *= x1
			ae += xe
		}
		x1 *= x1
		xe <<= 1
		if x1 < .5 {
			x1 += x1
			xe--
		}
	}
	return math.Ldexp(1/a1, -ae)
}

// kahan is a compensated running sum, as in pheniqs: c holds the low-order
// bits lost by each add and feeds them back into the next one.
type kahan struct{ sum, c float64 }

func (k *kahan) add(x float64) {
	y := x - k.c
	t := k.sum + y
	k.c = (t - k.sum) - y
	k.sum = t
}

type decoder struct {
	start, end int               // segment location in the padded R1
	codes      []uint8           // barcodes flattened, each end-start long
	ids        map[string]string // barcode sequence -> codec id without '@'
	seqs       []string
	conc       float64 // prior, (1-noise)/n
	noiseP     float64 // noise * random barcode probability
	randomP    float64 // 4^-length
	confidence float64
}

// decode returns the best barcode index and whether the call passed the PAMLD
// noise and confidence filters (false == QC fail).
//
// Every float operation matches the naive per-barcode loop (same operands, same
// order), so results are bit-identical; the speedup is only from hoisting the
// score lookups out of the barcode loop and scoring four barcodes at a time so
// their independent Kahan chains overlap in the CPU pipeline.
func (d *decoder) decode(obs, qual []byte) (best int, pass bool) {
	n := len(obs)
	var sc [64][4]float64 // sc[i][e]: score of expecting base e at position i
	for i := 0; i < n; i++ {
		o, q := baseCode[obs[i]], min(max(qual[i], 33)-33, 127)
		for e := 0; e < 4; e++ {
			sc[i][e] = score[q][e][o]
		}
	}
	nb := len(d.seqs)
	var stack [128]float64
	cond := stack[:min(nb, len(stack))]
	if nb > len(stack) {
		cond = make([]float64, nb)
	}
	b := 0
	for ; b+4 <= nb; b += 4 {
		c0, c1, c2, c3 := d.codes[b*n:(b+1)*n], d.codes[(b+1)*n:(b+2)*n], d.codes[(b+2)*n:(b+3)*n], d.codes[(b+3)*n:(b+4)*n]
		var k0, k1, k2, k3 kahan
		for i := 0; i < n; i++ {
			s := &sc[i]
			k0.add(s[c0[i]])
			k1.add(s[c1[i]])
			k2.add(s[c2[i]])
			k3.add(s[c3[i]])
		}
		cond[b] = pow10(-0.1 * k0.sum)
		cond[b+1] = pow10(-0.1 * k1.sum)
		cond[b+2] = pow10(-0.1 * k2.sum)
		cond[b+3] = pow10(-0.1 * k3.sum)
	}
	for ; b < nb; b++ {
		var k kahan
		code := d.codes[b*n : (b+1)*n]
		for i := 0; i < n; i++ {
			k.add(sc[i][code[i]])
		}
		cond[b] = pow10(-0.1 * k.sum)
	}
	best = -1
	var bestP, bestCond float64
	var sigma kahan
	for b, c := range cond {
		p := c * d.conc
		sigma.add(p)
		if p > bestP {
			best, bestP, bestCond = b, p, c
		}
	}
	if best < 0 || bestCond <= d.randomP {
		return -1, false
	}
	sigma.add(d.noiseP)
	return best, bestP/sigma.sum > d.confidence
}

// ── config ────────────────────────────────────────────────────────────────────

type codec struct {
	Barcode []string `json:"barcode"`
}

type config struct {
	Template struct {
		Transform struct{ Token []string } `json:"transform"`
	} `json:"template"`
	Molecular []struct {
		Base       string                   `json:"base"`
		Noise      float64                  `json:"noise"`
		Confidence float64                  `json:"confidence_threshold"`
		Transform  struct{ Token []string } `json:"transform"`
	} `json:"molecular"`
	Decoder map[string]struct {
		Codec map[string]codec `json:"codec"`
	} `json:"decoder"`
}

// parseToken reads a pheniqs "read:start:end" token ("" end means to the end of the read).
func parseToken(t string) (start, end int, err error) {
	f := strings.Split(t, ":")
	if len(f) != 3 || f[0] != "0" {
		return 0, 0, fmt.Errorf("unsupported token %q, expected 0:start:end", t)
	}
	if start, err = strconv.Atoi(f[1]); err != nil {
		return
	}
	if f[2] != "" {
		end, err = strconv.Atoi(f[2])
	}
	return
}

func loadConfig(path string) (decs []*decoder, r1Start int, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, 0, err
	}
	defer f.Close()
	var cfg config
	if err = json.NewDecoder(f).Decode(&cfg); err != nil {
		return
	}
	if len(cfg.Template.Transform.Token) == 0 {
		return nil, 0, errors.New("config has no template transform token")
	}
	if r1Start, _, err = parseToken(cfg.Template.Transform.Token[0]); err != nil {
		return
	}
	if len(cfg.Molecular) != 4 {
		return nil, 0, fmt.Errorf("expected 4 molecular decoders (D,A,B,C), got %d", len(cfg.Molecular))
	}
	for _, m := range cfg.Molecular {
		d := &decoder{ids: map[string]string{}}
		if len(m.Transform.Token) != 1 {
			return nil, 0, errors.New("each molecular decoder needs exactly one token")
		}
		if d.start, d.end, err = parseToken(m.Transform.Token[0]); err != nil {
			return
		}
		if d.end-d.start > 64 {
			return nil, 0, errors.New("barcode segment longer than 64 bases")
		}
		cdc := cfg.Decoder[m.Base].Codec
		keys := make([]string, 0, len(cdc))
		for k := range cdc {
			keys = append(keys, k)
		}
		// ponytail: pheniqs keeps file order for tie-breaking, sorted ids match it for the shipped config
		sort.Strings(keys)
		for _, k := range keys {
			if len(cdc[k].Barcode) == 0 {
				continue
			}
			bc := strings.ToUpper(cdc[k].Barcode[0])
			if len(bc) != d.end-d.start {
				return nil, 0, fmt.Errorf("barcode %s (%s) is not %d bases", k, bc, d.end-d.start)
			}
			for i := range bc {
				d.codes = append(d.codes, baseCode[bc[i]])
			}
			d.seqs = append(d.seqs, bc)
			d.ids[bc] = strings.TrimPrefix(k, "@")
		}
		if len(d.seqs) == 0 {
			return nil, 0, fmt.Errorf("no barcodes for decoder base %q", m.Base)
		}
		d.conc = (1 - m.Noise) / float64(len(d.seqs))
		d.randomP = math.Pow(4, -float64(d.end-d.start))
		d.noiseP = m.Noise * d.randomP
		d.confidence = m.Confidence
		decs = append(decs, d)
	}
	return decs, r1Start, nil
}

// ── ME search (as gih-stagger) ────────────────────────────────────────────────

func findME(seq, me []byte, maxMismatch float64) int {
	meLen := len(me)
	winEnd := min(65, len(seq)-meLen)
	if len(seq) < meLen || winStart > winEnd {
		return -1
	}
seqscan:
	for pos := winStart; pos <= winEnd; pos++ {
		mm := 0.0
		for i := range meLen {
			if r := seq[pos+i]; r == 'N' {
				mm += 0.3
			} else if r != me[i] {
				mm++
			}
			if mm > maxMismatch {
				continue seqscan
			}
		}
		return pos
	}
	return -1
}

// ── FASTQ reader (as gih-stagger, but qualities stay ASCII) ───────────────────

type fastqRecord struct{ name, seq, qual []byte }

type readPair struct{ fq1, fq2 fastqRecord }

type fastqReader struct{ r *bufio.Reader }

func (f *fastqReader) next() (fastqRecord, error) {
	header, err := f.r.ReadSlice('\n')
	if err != nil {
		if err == io.EOF && len(header) == 0 {
			return fastqRecord{}, io.EOF
		}
		return fastqRecord{}, fmt.Errorf("truncated header line: %w", err)
	}
	if header[0] != '@' {
		return fastqRecord{}, fmt.Errorf("malformed header, expected '@', got %q", header)
	}
	header = bytes.TrimRight(header[1:], "\r\n")
	if i := bytes.IndexByte(header, ' '); i >= 0 {
		header = header[:i]
	}
	name := string(header)
	// ReadSlice data is only valid until the next read, so copy each line as it arrives
	buf := append(make([]byte, 0, 512), name...)
	var ends [3]int
	for i := range ends {
		l, err := f.r.ReadSlice('\n')
		if err != nil {
			return fastqRecord{}, fmt.Errorf("truncated record %q: %w", name, err)
		}
		if i != 1 { // skip the '+' line
			buf = append(buf, bytes.TrimRight(l, "\r\n")...)
		}
		ends[i] = len(buf)
	}
	n, sEnd := len(name), ends[0]
	if sEnd-n != ends[2]-sEnd {
		return fastqRecord{}, fmt.Errorf("seq/qual length mismatch for %q: %d vs %d", name, sEnd-n, ends[2]-sEnd)
	}
	return fastqRecord{buf[:n:n], buf[n:sEnd:sEnd], buf[sEnd:ends[2]:ends[2]]}, nil
}

func openFastq(path string, gzBlocks int) (io.ReadCloser, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if !strings.HasSuffix(path, ".gz") {
		return f, nil
	}
	gz, err := pgzip.NewReaderN(f, 3<<20, gzBlocks)
	if err != nil {
		f.Close()
		return nil, err
	}
	return struct {
		io.Reader
		io.Closer
	}{gz, closers{gz, f}}, nil
}

type closers []io.Closer

func (c closers) Close() (err error) {
	for _, x := range c {
		if e := x.Close(); e != nil {
			err = e
		}
	}
	return
}

// ── per-pair processing ───────────────────────────────────────────────────────

type stats struct {
	total, discarded, qcfail, valid, corrected int64
}

func (s *stats) add(o stats) {
	s.total += o.total
	s.discarded += o.discarded
	s.qcfail += o.qcfail
	s.valid += o.valid
	s.corrected += o.corrected
}

// worker owns the scratch space for one goroutine.
type worker struct {
	decs        []*decoder
	r1Start     int
	me          []byte
	maxMismatch float64
	seq, qual   []byte
	bxSet       map[string]struct{}
	st          stats
}

// process handles one pair. keep is false when the ME was not found.
// It appends the two FASTQ records to o1/o2.
func (w *worker) process(p *readPair, o1, o2 *bytes.Buffer) {
	w.st.total++
	fq1, fq2 := &p.fq1, &p.fq2
	mePos := findME(fq1.seq, w.me, w.maxMismatch)
	if mePos == -1 {
		w.st.discarded++
		return
	}
	plen := 7
	if mePos >= 51 && mePos <= 58 {
		plen = mePos - 51
	}
	pad := padSeq[plen]
	bcFrom := 0
	if plen == 6 {
		bcFrom = 1
	}
	// padded R1 = pad + barcode region + sequence after the ME
	w.seq, w.qual = w.seq[:0], w.qual[:0]
	for range pad {
		w.qual = append(w.qual, asciiQ40)
	}
	w.seq = append(w.seq, pad...)
	w.seq = append(w.seq, fq1.seq[bcFrom:mePos]...)
	w.seq = append(w.seq, fq1.seq[mePos+len(w.me):]...)
	w.qual = append(w.qual, fq1.qual[bcFrom:mePos]...)
	w.qual = append(w.qual, fq1.qual[mePos+len(w.me):]...)

	// decode D, A, B, C
	var orig [4][]byte
	var corr [4]string
	qcFail, changed := false, false
	for i, d := range w.decs {
		orig[i] = w.seq[d.start:d.end]
		idx, pass := d.decode(orig[i], w.qual[d.start:d.end])
		if !pass {
			qcFail = true
		}
		if idx >= 0 {
			corr[i] = d.seqs[idx]
			changed = changed || corr[i] != string(orig[i])
		}
	}
	// failed reads keep the observed barcode, passing reads the corrected one (as gih-convert)
	var id [4]string
	valid := 1
	for i, d := range w.decs {
		s := corr[i]
		if qcFail {
			s = string(orig[i])
		}
		var ok bool
		if id[i], ok = d.ids[s]; !ok {
			id[i], valid = empty, 0
		}
	}
	bx := "A" + id[1] + "C" + id[2] + "B" + id[3] + "D" + id[0]
	w.bxSet[bx] = struct{}{}
	// counts are per read, and both mates carry the same tags
	w.st.valid += int64(2 * valid)
	if qcFail {
		w.st.qcfail += 2
	} else if changed {
		w.st.corrected += 2
	}

	writeFastq(o1, fq1.name, "/1", bx, valid, w.seq[w.r1Start:], w.qual[w.r1Start:])
	writeFastq(o2, fq2.name, "/2", bx, valid, fq2.seq, fq2.qual)
}

func writeFastq(b *bytes.Buffer, name []byte, dir, bx string, vx int, seq, qual []byte) {
	b.WriteByte('@')
	b.Write(name)
	b.WriteString(dir)
	b.WriteString("\tVX:i:")
	b.WriteByte(byte('0' + vx))
	b.WriteString("\tBX:Z:")
	b.WriteString(bx)
	b.WriteByte('\n')
	b.Write(seq)
	b.WriteString("\n+\n")
	b.Write(qual)
	b.WriteByte('\n')
}

// ── output ────────────────────────────────────────────────────────────────────

type gzOut struct {
	f   *os.File
	gz  *pgzip.Writer
	buf *bufio.Writer
}

func newGzOut(path string, threads int) (*gzOut, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	gz, err := pgzip.NewWriterLevel(f, 4)
	if err != nil {
		f.Close()
		return nil, err
	}
	gz.SetConcurrency(2<<20, threads)
	return &gzOut{f, gz, bufio.NewWriterSize(gz, 2<<20)}, nil
}

func (o *gzOut) close() error {
	return errors.Join(o.buf.Flush(), o.gz.Close(), o.f.Close())
}

// ── main ──────────────────────────────────────────────────────────────────────

type result struct{ r1, r2 *bytes.Buffer }

func fail(format string, a ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", a...)
	os.Exit(1)
}

func main() {
	meSeq := flag.String("me", "AGATGTGTATAAGAGACAG", "ME sequence to search for")
	statsfile := flag.String("stats", "stats.txt", "File name for ME stats output")
	maxMM := flag.Float64("max-mismatch", 2, "Maximum allowable ME mismatch score")
	nThreads := flag.Int("threads", 4, "Number of worker threads")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: gih-demux [options] <pheniqs.json> <R1.fq[.gz]> <R2.fq[.gz]> <out.R1.fq.gz> <out.R2.fq.gz>\n\nOptions:\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	args := flag.Args()
	if len(args) != 5 {
		flag.Usage()
		os.Exit(1)
	}
	threads := min(max(*nThreads, 1), runtime.NumCPU())

	decs, r1Start, err := loadConfig(args[0])
	if err != nil {
		fail("Error reading config: %v", err)
	}
	gzBlocks := max(threads/2, 1)
	fh1, err := openFastq(args[1], gzBlocks)
	if err != nil {
		fail("Error opening R1: %v", err)
	}
	defer fh1.Close()
	fh2, err := openFastq(args[2], gzBlocks)
	if err != nil {
		fail("Error opening R2: %v", err)
	}
	defer fh2.Close()
	rdr1 := &fastqReader{bufio.NewReaderSize(fh1, 1<<19)}
	rdr2 := &fastqReader{bufio.NewReaderSize(fh2, 1<<19)}

	out1, err := newGzOut(args[3], max(threads/2, 1))
	if err != nil {
		fail("Error opening R1 output: %v", err)
	}
	out2, err := newGzOut(args[4], max(threads/2, 1))
	if err != nil {
		fail("Error opening R2 output: %v", err)
	}

	// batches are decoded in parallel and written in input order
	type job struct {
		batch []readPair
		res   chan result
	}
	jobs := make(chan job, threads*2)
	ordered := make(chan chan result, threads*2)
	var wg, wwg sync.WaitGroup
	var mu sync.Mutex
	var tot stats
	bxSet := map[string]struct{}{}

	for range threads {
		wg.Go(func() {
			w := &worker{decs: decs, r1Start: r1Start, me: []byte(strings.ToUpper(*meSeq)), maxMismatch: *maxMM, bxSet: map[string]struct{}{}}
			for j := range jobs {
				r := result{new(bytes.Buffer), new(bytes.Buffer)}
				for i := range j.batch {
					w.process(&j.batch[i], r.r1, r.r2)
				}
				j.res <- r
			}
			mu.Lock()
			tot.add(w.st)
			for k := range w.bxSet {
				bxSet[k] = struct{}{}
			}
			mu.Unlock()
		})
	}
	var writeErr error
	wwg.Go(func() {
		for res := range ordered {
			r := <-res
			if writeErr == nil {
				_, e1 := out1.buf.Write(r.r1.Bytes())
				_, e2 := out2.buf.Write(r.r2.Bytes())
				writeErr = errors.Join(e1, e2)
			}
		}
	})

	submit := func(b []readPair) {
		res := make(chan result, 1)
		ordered <- res
		jobs <- job{b, res}
	}
	batch := make([]readPair, 0, batchSize)
	var npairs int64
	for {
		fq1, err1 := rdr1.next()
		fq2, err2 := rdr2.next()
		eof1, eof2 := errors.Is(err1, io.EOF), errors.Is(err2, io.EOF)
		if eof1 && eof2 {
			break
		}
		if err1 != nil && !eof1 {
			fail("Error: R1 malformed near pair %d: %v", npairs+1, err1)
		}
		if err2 != nil && !eof2 {
			fail("Error: R2 malformed near pair %d: %v", npairs+1, err2)
		}
		if eof1 != eof2 {
			fail("Error: R1/R2 read count mismatch near pair %d (R1 EOF=%v, R2 EOF=%v)", npairs+1, eof1, eof2)
		}
		npairs++
		batch = append(batch, readPair{fq1, fq2})
		if len(batch) == batchSize {
			submit(batch)
			batch = make([]readPair, 0, batchSize)
		}
	}
	if len(batch) > 0 {
		submit(batch)
	}
	close(jobs)
	wg.Wait()
	close(ordered)
	wwg.Wait()

	if err := errors.Join(writeErr, out1.close(), out2.close()); err != nil {
		fail("Error writing output: %v", err)
	}

	sf, err := os.Create(*statsfile)
	if err != nil {
		fail("Error creating stats file: %v", err)
	}
	fmt.Fprintf(sf, "Total read pairs processed:           %d\n", tot.total)
	fmt.Fprintf(sf, "Discarded (ME sequence not found):    %d\n", tot.discarded)
	sf.Close()

	fmt.Printf("Total unique barcodes:               %d\n", len(bxSet))
	fmt.Printf("Total reads with barcode failure:    %d\n", tot.qcfail)
	fmt.Printf("Total reads with valid barcodes:     %d\n", tot.valid)
	fmt.Printf("Total reads with corrected barcodes: %d\n", tot.corrected)
}

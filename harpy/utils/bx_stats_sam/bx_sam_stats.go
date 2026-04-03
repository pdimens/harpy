// bxstats – linked-read molecule metrics from a coordinate-sorted BAM/SAM file.
//
// Mirrors the behaviour of the Python bx_stats_sam script.
// Build:
//   go mod init bxstats
//   go get github.com/biogo/hts@latest
//   go build -o bxstats .
//
// Usage:
//   bxstats [-d <distance>] <input.bam>

package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ---------------------------------------------------------------------------
// ReadCloud
// ---------------------------------------------------------------------------

type position struct{ start, end int }

// ReadCloud accumulates alignment records sharing the same BX barcode.
// All slices are grown in-place to avoid repeated small allocations.
type ReadCloud struct {
	chromosome string
	positions  []position
	bp         []int
	inserts    []int
	count      []bool // true  → count this read (read1, or unpaired)
	barcode    string
	suffix     int
	valid      bool
}

func newReadCloud(valid bool) *ReadCloud {
	return &ReadCloud{
		positions: make([]position, 0, 16),
		bp:        make([]int, 0, 16),
		inserts:   make([]int, 0, 16),
		count:     make([]bool, 0, 16),
		valid:     valid,
	}
}

// add appends the relevant fields of a SAM record.
func (rc *ReadCloud) add(rec *sam.Record) {
	rc.chromosome = rec.Ref.Name()
	rc.bp = append(rc.bp, queryAlignmentLength(rec))
	isRead1 := rec.Flags&sam.Read1 != 0
	isPaired := rec.Flags&sam.Paired != 0
	rc.count = append(rc.count, isRead1 || !isPaired)

	if rc.valid {
		rc.positions = append(rc.positions, position{rec.Pos, rec.End()})
		rc.inserts = append(rc.inserts, insertSize(rec))
		if bx, err := bxTag(rec); err == nil {
			rc.barcode = bx
		}
	}
}

// reset clears per-molecule data while keeping barcode, suffix and valid.
func (rc *ReadCloud) reset() {
	rc.chromosome = ""
	rc.positions = rc.positions[:0]
	rc.bp = rc.bp[:0]
	rc.inserts = rc.inserts[:0]
	rc.count = rc.count[:0]
}

// ---------------------------------------------------------------------------
// sort helper (avoids reflect-based sort.Slice on hot path)
// ---------------------------------------------------------------------------

type byStart struct {
	pos     []position
	bp      []int
	inserts []int
	count   []bool
}

func (s byStart) Len() int           { return len(s.pos) }
func (s byStart) Less(i, j int) bool { return s.pos[i].start < s.pos[j].start }
func (s byStart) Swap(i, j int) {
	s.pos[i], s.pos[j] = s.pos[j], s.pos[i]
	s.bp[i], s.bp[j] = s.bp[j], s.bp[i]
	s.inserts[i], s.inserts[j] = s.inserts[j], s.inserts[i]
	s.count[i], s.count[j] = s.count[j], s.count[i]
}

// ---------------------------------------------------------------------------
// deconvolve / write
// ---------------------------------------------------------------------------

// deconvolve deconvolutes the cloud using cutoff, writes TSV rows to w, then
// resets the cloud.  w is expected to be a *bufio.Writer.
func (rc *ReadCloud) deconvolve(cutoff int, w *bufio.Writer) {
	if len(rc.positions) == 0 && len(rc.bp) == 0 {
		return
	}

	if !rc.valid {
		totalBP := sumInts(rc.bp)
		totalCount := countTrue(rc.count)
		w.WriteString(rc.chromosome)
		w.WriteByte('\t')
		w.WriteString("invalid\t")
		w.WriteString(statsLine(0, 0, 0, totalBP, totalCount))
		rc.reset()
		return
	}

	// sort by start position
	sort.Sort(byStart{rc.positions, rc.bp, rc.inserts, rc.count})

	start := rc.positions[0].start
	end := rc.positions[0].end
	insert := rc.inserts[0]
	bp := rc.bp[0]
	count := boolToInt(rc.count[0])

	writeRow := func(s, e, ins, b, c int) {
		bc := rc.barcode
		if rc.suffix > 0 {
			bc = fmt.Sprintf("%s-%d", rc.barcode, rc.suffix)
		}
		w.WriteString(rc.chromosome)
		w.WriteByte('\t')
		w.WriteString(bc)
		w.WriteByte('\t')
		w.WriteString(statsLine(s, e, ins, b, c))
	}

	for idx := 1; idx < len(rc.positions); idx++ {
		prevEnd := rc.positions[idx-1].end
		currStart := rc.positions[idx].start
		currEnd := rc.positions[idx].end

		if cutoff == 0 || (currStart-prevEnd) <= cutoff {
			if currEnd > end {
				end = currEnd
			}
			insert += rc.inserts[idx]
			bp += rc.bp[idx]
			count += boolToInt(rc.count[idx])
		} else {
			writeRow(start, end, insert, bp, count)
			rc.suffix++
			start = currStart
			end = currEnd
			insert = rc.inserts[idx]
			bp = rc.bp[idx]
			count = boolToInt(rc.count[idx])
		}
	}
	writeRow(start, end, insert, bp, count)
	rc.reset()
}

// ---------------------------------------------------------------------------
// Stats line builder  (avoids fmt.Sprintf per call)
// ---------------------------------------------------------------------------

var statsBuf [256]byte // goroutine-local scratch; safe because we have one writer goroutine

func statsLine(start, end, insert, bp, count int) string {
	inferred := end - start
	var covBP, covIns float64
	if inferred > 0 {
		covBP = math.Min(1.0, float64(bp)/float64(inferred))
		if covBP < 0 {
			covBP = 0
		}
		covIns = math.Min(1.0, float64(insert)/float64(inferred))
		if covIns < 0 {
			covIns = 0
		}
	}
	if count < 1 {
		count = 1
	}
	return fmt.Sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n",
		count, start, end, inferred, bp, insert, covBP, covIns)
}

// ---------------------------------------------------------------------------
// Insert size (mirrors Python logic)
// ---------------------------------------------------------------------------

func insertSize(rec *sam.Record) int {
	isPaired := rec.Flags&sam.Paired != 0
	isSupp := rec.Flags&sam.Supplementary != 0

	if isPaired {
		if isSupp {
			return queryAlignmentLength(rec)
		}
		tlen := rec.TempLen
		if tlen < 0 {
			tlen = -tlen
		}
		if tlen > 0 {
			return tlen
		}
		return 0 // max(0, tlen) where tlen came from rec.TempLen
	}
	// unpaired
	tlen := rec.TempLen
	if tlen < 0 {
		tlen = -tlen
	}
	ql := inferQueryLength(rec)
	if tlen > ql {
		return tlen
	}
	return ql
}

// queryAlignmentLength sums CIGAR consumed-on-reference ops (M/=/>X/D).
// biogo/hts/sam does not expose query_alignment_length directly.
func queryAlignmentLength(rec *sam.Record) int {
	n := 0
	for _, op := range rec.Cigar {
		switch op.Type() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch,
			sam.CigarInsertion, sam.CigarSoftClipped:
			n += op.Len()
		}
	}
	return n
}

// inferQueryLength counts bases that consume the query (M/I/S/=/X).
func inferQueryLength(rec *sam.Record) int {
	n := 0
	for _, op := range rec.Cigar {
		switch op.Type() {
		case sam.CigarMatch, sam.CigarInsertion, sam.CigarSoftClipped,
			sam.CigarEqual, sam.CigarMismatch:
			n += op.Len()
		}
	}
	return n
}

// ---------------------------------------------------------------------------
// Tag helpers
// ---------------------------------------------------------------------------

func bxTag(rec *sam.Record) (string, error) {
	aux, ok := rec.Tag([]byte("BX"))
	if !ok {
		return "", fmt.Errorf("no BX tag")
	}
	v := aux.Value()
	s, ok := v.(string)
	if !ok {
		return "", fmt.Errorf("BX tag not a string")
	}
	return s, nil
}

func vxTag(rec *sam.Record) int {
	aux, ok := rec.Tag([]byte("VX"))
	if !ok {
		return -1
	}
	switch v := aux.Value().(type) {
	case int8:
		return int(v)
	case uint8:
		return int(v)
	case int16:
		return int(v)
	case uint16:
		return int(v)
	case int32:
		return int(v)
	case uint32:
		return int(v)
	case int:
		return v
	}
	return -1
}

// ---------------------------------------------------------------------------
// Pool for ReadCloud objects
// ---------------------------------------------------------------------------

var cloudPool = sync.Pool{
	New: func() any { return newReadCloud(true) },
}

func acquireCloud(valid bool) *ReadCloud {
	rc := cloudPool.Get().(*ReadCloud)
	rc.valid = valid
	rc.suffix = 0
	rc.barcode = ""
	rc.reset()
	return rc
}

func releaseCloud(rc *ReadCloud) {
	cloudPool.Put(rc)
}

// ---------------------------------------------------------------------------
// writeStats flushes and frees all clouds in the map.
// ---------------------------------------------------------------------------

func writeStats(d map[string]*ReadCloud, cutoff int, w *bufio.Writer) {
	for key, cloud := range d {
		cloud.deconvolve(cutoff, w)
		releaseCloud(cloud)
		delete(d, key)
	}
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

func sumInts(s []int) int {
	n := 0
	for _, v := range s {
		n += v
	}
	return n
}

func countTrue(s []bool) int {
	n := 0
	for _, v := range s {
		if v {
			n++
		}
	}
	return n
}

func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

func main() {
	dist := flag.Int("d", 0, "distance threshold for linking alignments (0 = no deconvolution)")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: bxstats [-d <distance>] <input.bam>\n\n")
		fmt.Fprintf(os.Stderr, "Linked-read metrics from a coordinate-sorted BAM/SAM file.\n\n")
		flag.PrintDefaults()
	}
	flag.Parse()

	if flag.NArg() != 1 {
		flag.Usage()
		os.Exit(1)
	}

	inputPath := flag.Arg(0)
	f, err := os.Open(inputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "error opening file: %v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	br, err := bam.NewReader(f, 0) // 0 = single-threaded; increase for parallel decompression
	if err != nil {
		fmt.Fprintf(os.Stderr, "error creating BAM reader: %v\n", err)
		os.Exit(1)
	}
	defer br.Close()

	bw := bufio.NewWriterSize(os.Stdout, 1<<20) // 1 MiB output buffer
	defer bw.Flush()

	bw.WriteString("contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tinsert_len\tcoverage_bp\tcoverage_inserts\n")

	d := make(map[string]*ReadCloud, 1024)
	var lastContig string

	for {
		rec, err := br.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintf(os.Stderr, "error reading record: %v\n", err)
			os.Exit(1)
		}

		chrom := rec.Ref.Name()

		// When the contig changes, flush and free everything accumulated so far.
		if lastContig != "" && chrom != lastContig {
			writeStats(d, *dist, bw)
		}
		lastContig = chrom

		flags := rec.Flags
		// Skip duplicates, unmapped, secondary.
		if flags&sam.Duplicate != 0 || flags&sam.Unmapped != 0 || flags&sam.Secondary != 0 {
			continue
		}
		// Skip chimeric supplementary alignments (different contigs).
		if flags&sam.Supplementary != 0 && rec.Ref != rec.MateRef {
			continue
		}
		// Skip records with no CIGAR (unaligned).
		if len(rec.Cigar) == 0 {
			continue
		}

		bx, bxErr := bxTag(rec)
		vx := vxTag(rec)

		if bxErr != nil || vx == 0 {
			// invalid / missing BX or VX==0
			if _, ok := d["invalid"]; !ok {
				d["invalid"] = acquireCloud(false)
			}
			d["invalid"].add(rec)
			continue
		}

		cloud, ok := d[bx]
		if !ok {
			cloud = acquireCloud(true)
			d[bx] = cloud
		}
		cloud.add(rec)
	}

	// Flush the final contig.
	writeStats(d, *dist, bw)
}

package main

import (
	"bytes"
	"strings"
	"testing"
)

const cfgPath = "../../templates/pheniqs.config.json"

// read builds a stagger-01 (7 T pad) R1: padded[7:58] + ME + tail, with the given A/B/C segments.
func read(t *testing.T, a, b, c string, q byte) (readPair, string) {
	t.Helper()
	const me = "AGATGTGTATAAGAGACAG"
	padded := "TTTTTTT" + "GGG" + a + "ACA" + b + "ACA" + c + "ACAGTA" // 7+3+12+3+12+3+12+6 = 58
	if len(padded) != 58 {
		t.Fatalf("bad test read, padded is %d", len(padded))
	}
	r1 := padded[7:] + me + strings.Repeat("ACGT", 20)
	r2 := strings.Repeat("TGCA", 20)
	qual := func(n int) []byte { return bytes.Repeat([]byte{q}, n) }
	return readPair{
		fastqRecord{[]byte("r1"), []byte(r1), qual(len(r1))},
		fastqRecord{[]byte("r1"), []byte(r2), qual(len(r2))},
	}, padded
}

func run(t *testing.T, p readPair) (string, stats) {
	t.Helper()
	decs, r1Start, err := loadConfig(cfgPath)
	if err != nil {
		t.Fatal(err)
	}
	w := &worker{decs: decs, r1Start: r1Start, me: []byte("AGATGTGTATAAGAGACAG"), maxMismatch: 2, bxSet: map[string]struct{}{}}
	var o1, o2 bytes.Buffer
	w.process(&p, &o1, &o2)
	return o1.String(), w.st
}

func TestExactAndCorrected(t *testing.T) {
	d, _, _ := loadConfig(cfgPath)
	a, b, c := d[1].seqs[0], d[2].seqs[1], d[3].seqs[2] // segment @01, @02, @03

	p, padded := read(t, a, b, c, 'I')
	out, st := run(t, p)
	if !strings.Contains(out, "\tVX:i:1\tBX:Z:A01C02B03D01\n") || st.qcfail != 0 || st.corrected != 0 {
		t.Fatalf("exact barcode: %q %+v", out, st)
	}
	// R1 sequence is padded[58:], the 6 bases after the barcode region plus ME-trimmed tail
	if want := strings.Repeat("ACGT", 20) + "\n"; !strings.Contains(out, "\n"+want+"+\n") {
		t.Fatalf("R1 trim: %q (padded %s)", out, padded)
	}

	// one low-quality mismatch in A is corrected
	mut := []byte(a)
	mut[3] = map[byte]byte{'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[3]]
	p, _ = read(t, string(mut), b, c, 'I')
	p.fq1.qual[3+3] = '#' // Q2 at the mutated base (A starts at padded index 10 -> raw index 3)
	out, st = run(t, p)
	if !strings.Contains(out, "BX:Z:A01C02B03D01\n") || st.corrected != 2 || st.qcfail != 0 {
		t.Fatalf("corrected barcode: %q %+v", out, st)
	}
}

func TestFailKeepsObserved(t *testing.T) {
	d, _, _ := loadConfig(cfgPath)
	// six high-quality mismatches in A: no confident call, observed (invalid) barcode kept
	mut := []byte(d[1].seqs[0])
	for i := 0; i < 6; i++ {
		mut[i] = map[byte]byte{'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[i]]
	}
	p, _ := read(t, string(mut), d[2].seqs[1], d[3].seqs[2], 'I')
	out, st := run(t, p)
	if !strings.Contains(out, "\tVX:i:0\tBX:Z:A00C02B03D01\n") || st.qcfail != 2 || st.valid != 0 {
		t.Fatalf("failed barcode: %q %+v", out, st)
	}
}

func TestNoME(t *testing.T) {
	p, _ := read(t, strings.Repeat("A", 12), strings.Repeat("C", 12), strings.Repeat("G", 12), 'I')
	p.fq1.seq = bytes.Repeat([]byte("A"), 150)
	out, st := run(t, p)
	if out != "" || st.discarded != 1 {
		t.Fatalf("no ME: %q %+v", out, st)
	}
}

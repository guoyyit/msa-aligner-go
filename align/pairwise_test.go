package align

import (
	"strings"
	"testing"
)

func ungap(b []byte) string {
	return strings.ReplaceAll(string(b), "-", "")
}

func scoreOf(a, b []byte, p Params) int {
	s := 0
	for k := range a {
		if a[k] != '-' && b[k] != '-' {
			s += p.Subst.Score(a[k], b[k])
		}
	}
	// A gap run of length L costs open + L*extend (Gotoh convention).
	run := func(row []byte) {
		for k := 0; k < len(row); {
			if row[k] != '-' {
				k++
				continue
			}
			start := k
			for k < len(row) && row[k] == '-' {
				k++
			}
			s += p.GapOpen + p.GapExtend*(k-start)
		}
	}
	run(a)
	run(b)
	return s
}

func TestPairwiseIdentical(t *testing.T) {
	p := DefaultParams()
	al, bl, score := Pairwise([]byte("ACGT"), []byte("ACGT"), p)
	if string(al) != "ACGT" || string(bl) != "ACGT" {
		t.Fatalf("identical seqs must align without gaps: %q vs %q", al, bl)
	}
	if score != 4*2 {
		t.Fatalf("score = %d, want 8", score)
	}
}

func TestPairwisePreservesResidueOrder(t *testing.T) {
	p := DefaultParams()
	cases := [][2]string{
		{"ACGT", "AGT"},
		{"AAAA", "AATAA"},
		{"GATTACA", "GCATGCU"},
		{"", "ACGT"},
		{"ACGT", ""},
	}
	for _, c := range cases {
		al, bl, score := Pairwise([]byte(c[0]), []byte(c[1]), p)
		if len(al) != len(bl) {
			t.Fatalf("lengths differ for %v", c)
		}
		if ungap(al) != c[0] || ungap(bl) != c[1] {
			t.Fatalf("residues changed: %q->%q, %q->%q", c[0], al, c[1], bl)
		}
		if got := scoreOf(al, bl, p); got != score {
			t.Fatalf("reported score %d != recomputed %d for %v", score, got, c)
		}
	}
}

func TestPairwiseDeterminism(t *testing.T) {
	p := DefaultParams()
	var refA, refB string
	for run := 0; run < 5; run++ {
		al, bl, _ := Pairwise([]byte("GATTACA"), []byte("GCATGCA"), p)
		if run == 0 {
			refA, refB = string(al), string(bl)
			continue
		}
		if string(al) != refA || string(bl) != refB {
			t.Fatalf("non-deterministic alignment: %q/%q vs %q/%q", al, bl, refA, refB)
		}
	}
}

func TestPairwiseParameterChange(t *testing.T) {
	a := []byte("ACGTACGT")
	b := []byte("ACGTTCGT")
	strict := Params{Subst: NewSubstMatrix(2, -1), GapOpen: -20, GapExtend: -10}
	loose := Params{Subst: NewSubstMatrix(2, -3), GapOpen: -1, GapExtend: 0}
	as, bs, _ := Pairwise(a, b, strict)
	al, bl, _ := Pairwise(a, b, loose)
	// With cheap gaps, a gap should be preferred over the central mismatch.
	if strings.Count(string(al), "-") == 0 {
		t.Fatalf("cheap-gap params should prefer gaps, got %q vs %q", al, bl)
	}
	// With expensive gaps, the alignment must stay gap-free here.
	if strings.Count(string(as), "-") != 0 {
		t.Fatalf("expensive-gap params should avoid gaps, got %q vs %q", as, bs)
	}
}

package msa

import (
	"testing"
)

func testConfig() Config {
	return Config{
		Matrix:    BLOSUM62(),
		GapOpen:   -10,
		GapExtend: -1,
		KmerSize:  3,
	}
}

func TestAlignProfilesIdentical(t *testing.T) {
	cfg := testConfig()
	a := [][]byte{[]byte("ACGTACGT")}
	b := [][]byte{[]byte("ACGTACGT")}
	out, err := alignProfiles(a, b, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(out) != 2 {
		t.Fatalf("rows = %d, want 2", len(out))
	}
	if string(out[0]) != "ACGTACGT" || string(out[1]) != "ACGTACGT" {
		t.Errorf("got %q %q", out[0], out[1])
	}
}

func TestAlignProfilesWithGap(t *testing.T) {
	cfg := testConfig()
	a := [][]byte{[]byte("ACGT")}
	b := [][]byte{[]byte("CGT")}
	out, err := alignProfiles(a, b, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(out) != 2 {
		t.Fatalf("rows = %d", len(out))
	}
	if len(out[0]) != len(out[1]) {
		t.Errorf("ragged: %d vs %d", len(out[0]), len(out[1]))
	}
	// Should introduce a gap to align the missing 'A'.
	hasGap := false
	for _, r := range out {
		for _, c := range r {
			if c == '-' {
				hasGap = true
			}
		}
	}
	if !hasGap {
		t.Errorf("expected a gap, got %q %q", out[0], out[1])
	}
}

func TestAlignProfilesRagged(t *testing.T) {
	cfg := testConfig()
	a := [][]byte{[]byte("AC"), []byte("ACG")}
	b := [][]byte{[]byte("GT")}
	if _, err := alignProfiles(a, b, cfg); err == nil {
		t.Error("expected ragged-profile error")
	}
}

func TestAlignProgressiveDeterministic(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "ACGTACGTAC"},
		{ID: "c", Seq: "GGGTTTAAAC"},
		{ID: "d", Seq: "TTTTGGGGCC"},
	}
	cfg := testConfig()
	dm, _ := ComputeDistances(seqs, 3, 1, nil)
	order, err := BuildMergeOrder(dm)
	if err != nil {
		t.Fatal(err)
	}
	a1, err := AlignProgressive(seqs, order, cfg, nil)
	if err != nil {
		t.Fatal(err)
	}
	a2, err := AlignProgressive(seqs, order, cfg, nil)
	if err != nil {
		t.Fatal(err)
	}
	if a1.Len != a2.Len {
		t.Fatalf("lengths differ: %d vs %d", a1.Len, a2.Len)
	}
	for i := range a1.Seq {
		if a1.Seq[i] != a2.Seq[i] {
			t.Errorf("row %d differs:\n%q\n%q", i, a1.Seq[i], a2.Seq[i])
		}
	}
}

func TestAlignProgressiveSingle(t *testing.T) {
	seqs := []Sequence{{ID: "a", Seq: "ACGT"}}
	order := &MergeOrder{N: 1}
	a, err := AlignProgressive(seqs, order, testConfig(), nil)
	if err != nil {
		t.Fatal(err)
	}
	if a.Len != 4 || a.Seq[0] != "ACGT" {
		t.Errorf("got len=%d seq=%q", a.Len, a.Seq[0])
	}
}

func TestAlignProgressiveAllRowsEqualLength(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGT"},
		{ID: "b", Seq: "ACGT"},
		{ID: "c", Seq: "GGGTTT"},
	}
	cfg := testConfig()
	dm, _ := ComputeDistances(seqs, 3, 1, nil)
	order, _ := BuildMergeOrder(dm)
	a, err := AlignProgressive(seqs, order, cfg, nil)
	if err != nil {
		t.Fatal(err)
	}
	for i, row := range a.Seq {
		if len(row) != a.Len {
			t.Errorf("row %d len = %d, want %d", i, len(row), a.Len)
		}
	}
}

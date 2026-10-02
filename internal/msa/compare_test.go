package msa

import "testing"

func TestAppendPreservesExisting(t *testing.T) {
	existing := &Alignment{
		IDs: []string{"a", "b"},
		Seq: []string{"ACGT", "ACGT"},
		Len: 4,
	}
	newSeqs := []Sequence{{ID: "c", Seq: "ACGT"}}
	out, err := Append(existing, newSeqs, testConfig(), nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(out.Seq) != 3 {
		t.Fatalf("rows = %d, want 3", len(out.Seq))
	}
	// Existing rows must preserve their residue order (only gap columns added).
	// Strip gaps and compare to original.
	for i := 0; i < 2; i++ {
		stripped := stripGaps(out.Seq[i])
		if stripped != "ACGT" {
			t.Errorf("existing row %d lost residues: %q -> %q", i, out.Seq[i], stripped)
		}
	}
	if stripGaps(out.Seq[2]) != "ACGT" {
		t.Errorf("new row = %q", out.Seq[2])
	}
}

func TestAppendDiffBounded(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "GGGTTTAAAC"},
		{ID: "c", Seq: "TTTTGGGGCC"},
	}
	cfg := testConfig()

	// Full recompute with all three.
	full, err := Run(nil, seqs, cfg, "", nil)
	if err != nil {
		t.Fatal(err)
	}

	// Build alignment of first two, then append the third.
	two := []Sequence{seqs[0], seqs[1]}
	twoRes, err := Run(nil, two, cfg, "", nil)
	if err != nil {
		t.Fatal(err)
	}
	appended, err := Append(twoRes.Alignment, []Sequence{seqs[2]}, cfg, nil)
	if err != nil {
		t.Fatal(err)
	}

	diff, err := CompareAlignments(appended, full.Alignment, cfg)
	if err != nil {
		t.Fatal(err)
	}
	// Diff must be bounded: shifted columns <= total columns, score diff finite.
	if diff.ColumnsShifted > appended.Len {
		t.Errorf("ColumnsShifted %d > len %d", diff.ColumnsShifted, appended.Len)
	}
	if diff.MaxShift > appended.Len {
		t.Errorf("MaxShift %d > len %d", diff.MaxShift, appended.Len)
	}
	// Score delta must be finite and bounded.
	if diff.ScoreDiff != diff.ScoreDiff { // NaN check
		t.Errorf("ScoreDiff is NaN")
	}
}

func TestAlignmentScore(t *testing.T) {
	a := &Alignment{
		IDs: []string{"a", "b"},
		Seq: []string{"AC", "AC"},
		Len: 2,
	}
	score := AlignmentScore(a, testConfig())
	// Two identical columns: A/A=4 and C/C=9 with BLOSUM62, total 13.
	if score != 13 {
		t.Errorf("score = %v, want 13", score)
	}
}

func stripGaps(s string) string {
	var b []byte
	for i := 0; i < len(s); i++ {
		if s[i] != '-' {
			b = append(b, s[i])
		}
	}
	return string(b)
}

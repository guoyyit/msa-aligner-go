package msa

import (
	"math"
	"testing"
)

func TestPairDistanceIdentical(t *testing.T) {
	d, err := PairDistance("ACGTACGT", "ACGTACGT", 3)
	if err != nil {
		t.Fatal(err)
	}
	if d != 0 {
		t.Errorf("identical distance = %v, want 0", d)
	}
}

func TestPairDistanceDifferent(t *testing.T) {
	// Completely non-overlapping k-mers -> distance near 1.
	d, err := PairDistance("AAAAAA", "CCCCCC", 2)
	if err != nil {
		t.Fatal(err)
	}
	if d < 0.9 {
		t.Errorf("different distance = %v, want ~1", d)
	}
}

func TestPairDistanceSymmetric(t *testing.T) {
	a := "ACGTACGTAC"
	b := "GGGTTTAAAC"
	d1, _ := PairDistance(a, b, 3)
	d2, _ := PairDistance(b, a, 3)
	if math.Abs(d1-d2) > 1e-12 {
		t.Errorf("asymmetric: %v vs %v", d1, d2)
	}
}

func TestPairDistanceShortSequence(t *testing.T) {
	// Sequence shorter than k degrades to k=1, should not error.
	d, err := PairDistance("AC", "ACGTACGT", 5)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if d < 0 || d > 1 {
		t.Errorf("distance out of range: %v", d)
	}
}

func TestComputeDistancesMatrix(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGT"},
		{ID: "b", Seq: "ACGTACGT"},
		{ID: "c", Seq: "GGGTTTAA"},
	}
	dm, errs := ComputeDistances(seqs, 3, 2, nil)
	if len(errs) != 0 {
		t.Fatalf("unexpected errors: %v", errs)
	}
	if dm.N != 3 {
		t.Fatalf("N = %d, want 3", dm.N)
	}
	if dm.Dist[0][1] != 0 {
		t.Errorf("Dist[0][1] = %v, want 0", dm.Dist[0][1])
	}
	if dm.Dist[0][1] != dm.Dist[1][0] {
		t.Error("matrix not symmetric")
	}
	if dm.Dist[0][0] != 0 {
		t.Error("diagonal not zero")
	}
	if dm.Dist[0][2] < 0.5 {
		t.Errorf("Dist[0][2] = %v, want large", dm.Dist[0][2])
	}
}

func TestComputeDistancesDeterministic(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "GGGTTTAAAC"},
		{ID: "c", Seq: "TTTTGGGGCC"},
		{ID: "d", Seq: "ACGTACGTAC"},
	}
	dm1, _ := ComputeDistances(seqs, 3, 1, nil)
	dm2, _ := ComputeDistances(seqs, 3, 8, nil)
	for i := 0; i < dm1.N; i++ {
		for j := 0; j < dm1.N; j++ {
			if dm1.Dist[i][j] != dm2.Dist[i][j] {
				t.Errorf("Dist[%d][%d] differs: %v vs %v", i, j, dm1.Dist[i][j], dm2.Dist[i][j])
			}
		}
	}
}

func TestComputeDistancesProgress(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGT"},
		{ID: "b", Seq: "ACGT"},
		{ID: "c", Seq: "GGGT"},
	}
	var calls int
	dm, errs := ComputeDistances(seqs, 2, 2, func(p Progress) {
		calls++
		if p.Phase != "distance" {
			t.Errorf("phase = %q", p.Phase)
		}
	})
	if len(errs) != 0 {
		t.Fatalf("errors: %v", errs)
	}
	if calls != 3 {
		t.Errorf("progress calls = %d, want 3", calls)
	}
	_ = dm
}

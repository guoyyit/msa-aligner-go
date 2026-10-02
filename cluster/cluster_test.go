package cluster

import (
	"math"
	"testing"

	"github.com/guoyyit/msa-aligner-go/distance"
)

func matrixFrom(ids []string, cells []float64) *distance.Matrix {
	return &distance.Matrix{N: len(ids), IDs: ids, Cells: cells, Errors: map[[2]int]string{}}
}

func TestBuildTrivial(t *testing.T) {
	m := matrixFrom([]string{"a"}, []float64{0})
	tree, err := Build(m)
	if err != nil {
		t.Fatal(err)
	}
	if len(tree.Steps) != 0 || tree.N != 1 {
		t.Fatalf("single leaf must have no merges, got %+v", tree.Steps)
	}
}

func TestBuildDeterministicOrder(t *testing.T) {
	// Distances: a-b=0.1 (closest), a-c=0.9, b-c=0.8.
	m := matrixFrom([]string{"a", "b", "c"}, []float64{
		0, 0.1, 0.9,
		0.1, 0, 0.8,
		0.9, 0.8, 0,
	})
	tree, err := Build(m)
	if err != nil {
		t.Fatal(err)
	}
	if len(tree.Steps) != 2 {
		t.Fatalf("want 2 steps, got %d", len(tree.Steps))
	}
	if s := tree.Steps[0]; s.Left != 0 || s.Right != 1 || s.NewIndex != 3 {
		t.Fatalf("first merge wrong: %+v", s)
	}
	// Cluster {a,b} vs c = (0.9+0.8)/2 = 0.85.
	if d := tree.Steps[1].Distance; math.Abs(d-0.85) > 1e-12 {
		t.Fatalf("average linkage distance = %v", d)
	}
	if got := tree.Members[tree.Steps[1].NewIndex]; len(got) != 3 {
		t.Fatalf("root must contain all leaves: %v", got)
	}
}

func TestBuildTieBreaksByIndex(t *testing.T) {
	m := matrixFrom([]string{"a", "b", "c", "d"}, []float64{
		0, 1, 1, 1,
		1, 0, 1, 1,
		1, 1, 0, 1,
		1, 1, 1, 0,
	})
	tree, err := Build(m)
	if err != nil {
		t.Fatal(err)
	}
	s := tree.Steps[0]
	if s.Left != 0 || s.Right != 1 {
		t.Fatalf("tie must break toward smallest pair, got (%d,%d)", s.Left, s.Right)
	}

	// Re-running must produce byte-identical steps.
	tree2, _ := Build(m)
	for i := range tree.Steps {
		if tree.Steps[i] != tree2.Steps[i] {
			t.Fatalf("non-deterministic build at step %d", i)
		}
	}
}

func TestBuildNaNDeferred(t *testing.T) {
	nan := math.NaN()
	m := matrixFrom([]string{"a", "b", "c"}, []float64{
		0, nan, 0.5,
		nan, 0, nan,
		0.5, nan, 0,
	})
	tree, err := Build(m)
	if err != nil {
		t.Fatal(err)
	}
	s0 := tree.Steps[0]
	if s0.Left != 0 || s0.Right != 2 {
		t.Fatalf("first merge must be finite pair (0,2), got %+v", s0)
	}
	if math.IsNaN(s0.Distance) {
		t.Fatalf("first merge distance must be finite, got NaN")
	}
	if !math.IsNaN(tree.Steps[1].Distance) {
		t.Fatalf("merge with unavailable pairs must record NaN, got %v", tree.Steps[1].Distance)
	}
}

func TestBuildAllFailed(t *testing.T) {
	nan := math.NaN()
	m := matrixFrom([]string{"a", "b"}, []float64{0, nan, nan, 0})
	if _, err := Build(m); err == nil {
		t.Fatal("all-NaN matrix must be rejected")
	}
}

func TestBuildEmpty(t *testing.T) {
	if _, err := Build(nil); err == nil {
		t.Fatal("nil matrix must be rejected")
	}
}

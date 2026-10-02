package distance

import (
	"math"
	"testing"

	"github.com/guoyyit/msa-aligner-go/fasta"
)

func seqsFromStrings(idData ...string) []fasta.Sequence {
	var out []fasta.Sequence
	for i := 0; i+1 < len(idData); i += 2 {
		out = append(out, fasta.Sequence{ID: idData[i], Data: []byte(idData[i+1])})
	}
	return out
}

func TestComputeBasic(t *testing.T) {
	seqs := seqsFromStrings("a", "ACGTACGT", "b", "ACGTACGT", "c", "TTTTAAAA")
	m, err := Compute(seqs, Options{Concurrency: 2}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if d := m.At(0, 1); d != 0 {
		t.Fatalf("identical sequences distance = %v, want 0", d)
	}
	if d := m.At(1, 0); d != m.At(0, 1) {
		t.Fatal("matrix not symmetric")
	}
	if d := m.At(0, 2); d <= 0 || d > 1 {
		t.Fatalf("distance out of range: %v", d)
	}
	if len(m.Errors) != 0 {
		t.Fatalf("unexpected errors: %v", m.Errors)
	}
	if math.IsNaN(m.At(0, 1)) {
		t.Fatal("diagonal-adjacent NaN")
	}
}

func TestComputeDeterminismAcrossConcurrency(t *testing.T) {
	seqs := seqsFromStrings(
		"a", "ACGTACGTAC", "b", "ACGTTCGTAC",
		"c", "TTTTAAAAGG", "d", "ACGTACGTAA",
	)
	var refs []float64
	for _, workers := range []int{1, 2, 8} {
		m, err := Compute(seqs, Options{Concurrency: workers}, nil)
		if err != nil {
			t.Fatal(err)
		}
		if refs == nil {
			refs = append([]float64(nil), m.Cells...)
			continue
		}
		for i := range refs {
			if m.Cells[i] != refs[i] {
				t.Fatalf("cell %d differs with concurrency %d: %v vs %v", i, workers, m.Cells[i], refs[i])
			}
		}
	}
}

func TestComputeOneBadPairDoesNotBlock(t *testing.T) {
	// Sequence "AC" is shorter than k=3, so pairs involving it fail,
	// but all other pairs still complete.
	seqs := seqsFromStrings("good1", "ACGTACGT", "tiny", "AC", "good2", "TTTTAAAA")
	m, err := Compute(seqs, Options{Concurrency: 2}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := m.Failed(0, 1); !ok {
		t.Fatal("expected failure for pair (good1,tiny)")
	}
	if _, ok := m.Failed(1, 2); !ok {
		t.Fatal("expected failure for pair (tiny,good2)")
	}
	if _, ok := m.Failed(0, 2); ok {
		t.Fatal("pair of two good sequences must not fail")
	}
	if math.IsNaN(m.At(0, 2)) {
		t.Fatal("good pair distance must be computed")
	}
	if d01 := m.At(0, 1); !math.IsNaN(d01) {
		t.Fatalf("failed cell must be NaN, got %v", d01)
	}
	failed := m.OrderedErrors()
	if len(failed) != 2 || failed[0].I != 0 || failed[1].I != 1 {
		t.Fatalf("unexpected ordered failures: %+v", failed)
	}
}

func TestComputeProgress(t *testing.T) {
	seqs := seqsFromStrings("a", "ACGTACGT", "b", "ACGTACGT", "c", "ACGTACGT")
	var lastDone, sawTotal int
	m, err := Compute(seqs, Options{Concurrency: 1}, func(done, total int) {
		lastDone = done
		sawTotal = total
	})
	if err != nil {
		t.Fatal(err)
	}
	if lastDone != 3 || sawTotal != 3 {
		t.Fatalf("progress ended at %d/%d", lastDone, sawTotal)
	}
	_ = m
}

func TestComputeEmpty(t *testing.T) {
	if _, err := Compute(nil, Options{}, nil); err == nil {
		t.Fatal("expected error for empty input")
	}
}

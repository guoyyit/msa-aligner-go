package msa

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"testing"
)

func TestRunEndToEnd(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "ACGTACGTAC"},
		{ID: "c", Seq: "GGGTTTAAAC"},
		{ID: "d", Seq: "TTTTGGGGCC"},
	}
	cfg := testConfig()
	res, err := Run(context.Background(), seqs, cfg, "", nil)
	if err != nil {
		t.Fatal(err)
	}
	if res.Alignment == nil || res.Alignment.Len == 0 {
		t.Fatal("empty alignment")
	}
	if len(res.Alignment.Seq) != 4 {
		t.Errorf("rows = %d, want 4", len(res.Alignment.Seq))
	}
	for i, row := range res.Alignment.Seq {
		if len(row) != res.Alignment.Len {
			t.Errorf("row %d len = %d, want %d", i, len(row), res.Alignment.Len)
		}
	}
	if len(res.Consensus) != res.Alignment.Len {
		t.Errorf("consensus len = %d, want %d", len(res.Consensus), res.Alignment.Len)
	}
	if len(res.Columns) != res.Alignment.Len {
		t.Errorf("columns = %d, want %d", len(res.Columns), res.Alignment.Len)
	}
}

func TestRunDeterministic(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "GGGTTTAAAC"},
		{ID: "c", Seq: "TTTTGGGGCC"},
	}
	cfg := testConfig()
	r1, _ := Run(context.Background(), seqs, cfg, "", nil)
	r2, _ := Run(context.Background(), seqs, cfg, "", nil)
	for i := range r1.Alignment.Seq {
		if r1.Alignment.Seq[i] != r2.Alignment.Seq[i] {
			t.Errorf("row %d differs:\n%q\n%q", i, r1.Alignment.Seq[i], r2.Alignment.Seq[i])
		}
	}
}

func TestRunResourceLimits(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGT"},
		{ID: "b", Seq: "ACGT"},
	}
	cfg := testConfig()
	cfg.MaxSequences = 1
	_, err := Run(context.Background(), seqs, cfg, "", nil)
	if !errors.Is(err, ErrLimitExceeded) {
		t.Errorf("err = %v, want ErrLimitExceeded", err)
	}

	cfg = testConfig()
	cfg.MaxLength = 2
	_, err = Run(context.Background(), seqs, cfg, "", nil)
	if !errors.Is(err, ErrLimitExceeded) {
		t.Errorf("err = %v, want ErrLimitExceeded", err)
	}

	cfg = testConfig()
	cfg.MaxTotalRes = 3
	_, err = Run(context.Background(), seqs, cfg, "", nil)
	if !errors.Is(err, ErrLimitExceeded) {
		t.Errorf("err = %v, want ErrLimitExceeded", err)
	}

	cfg = testConfig()
	cfg.MaxMemory = 1
	_, err = Run(context.Background(), seqs, cfg, "", nil)
	if !errors.Is(err, ErrLimitExceeded) {
		t.Errorf("err = %v, want ErrLimitExceeded", err)
	}
}

func TestRunCheckpointResume(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGTACGTAC"},
		{ID: "b", Seq: "GGGTTTAAAC"},
		{ID: "c", Seq: "TTTTGGGGCC"},
	}
	cfg := testConfig()
	dir := t.TempDir()

	r1, err := Run(context.Background(), seqs, cfg, dir, nil)
	if err != nil {
		t.Fatal(err)
	}
	cpPath := filepath.Join(dir, "alignment.json")
	if _, err := os.Stat(cpPath); err != nil {
		t.Fatalf("checkpoint not written: %v", err)
	}

	// Second run resumes from checkpoint; result must be identical.
	r2, err := Run(context.Background(), seqs, cfg, dir, nil)
	if err != nil {
		t.Fatal(err)
	}
	for i := range r1.Alignment.Seq {
		if r1.Alignment.Seq[i] != r2.Alignment.Seq[i] {
			t.Errorf("resumed row %d differs:\n%q\n%q", i, r1.Alignment.Seq[i], r2.Alignment.Seq[i])
		}
	}
}

func TestRunCheckpointInputMismatch(t *testing.T) {
	seqs1 := []Sequence{
		{ID: "a", Seq: "ACGT"},
		{ID: "b", Seq: "ACGT"},
	}
	seqs2 := []Sequence{
		{ID: "a", Seq: "ACGT"},
		{ID: "b", Seq: "GGGT"},
	}
	cfg := testConfig()
	dir := t.TempDir()
	if _, err := Run(context.Background(), seqs1, cfg, dir, nil); err != nil {
		t.Fatal(err)
	}
	_, err := Run(context.Background(), seqs2, cfg, dir, nil)
	if !errors.Is(err, ErrCheckpointInput) {
		t.Errorf("err = %v, want ErrCheckpointInput", err)
	}
}

func TestInputHashDeterministic(t *testing.T) {
	seqs := []Sequence{
		{ID: "a", Seq: "ACGT"},
		{ID: "b", Seq: "GGGT"},
	}
	h1 := InputHash(seqs)
	h2 := InputHash(seqs)
	if h1 != h2 {
		t.Errorf("hash not deterministic: %s vs %s", h1, h2)
	}
	seqs[1].Seq = "GGGG"
	if InputHash(seqs) == h1 {
		t.Error("hash did not change with input")
	}
}

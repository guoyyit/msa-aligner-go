package pipeline

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/guoyyit/msa-aligner-go/checkpoint"
)

func cfgWithCheckpoint(t *testing.T) (Config, string) {
	cfg := DefaultConfig()
	cfg.Concurrency = 2
	cfg.CheckpointPath = filepath.Join(t.TempDir(), "state.ckpt")
	cfg.CheckpointInterval = 20 * time.Millisecond
	return cfg, cfg.CheckpointPath
}

func TestResumeSkipsCompletedDistancePairs(t *testing.T) {
	cfg, path := cfgWithCheckpoint(t)

	// First run: cancel during the distance stage.
	ctx, cancel := context.WithCancel(context.Background())
	svc := New(cfg)
	cancelled := false
	_, err := svc.Align(ctx, strings.NewReader(fastaInput), func(p Progress) {
		// Cancel once at least one pair is done but not all, so there is
		// both persisted and remaining work.
		if p.Stage == stageDistance && p.Done >= 1 && p.Done < p.Total && !cancelled {
			cancelled = true
			cancel()
		}
	})
	if err == nil {
		t.Fatal("expected cancellation error")
	}
	if _, statErr := os.Stat(path); statErr != nil {
		t.Fatalf("checkpoint not written: %v", statErr)
	}

	// Capture how many pairs were persisted.
	st, err := checkpoint.Load(path)
	if err != nil {
		t.Fatal(err)
	}
	var d jobData
	if err := st.DecodeData(&d); err != nil {
		t.Fatal(err)
	}
	persistedBefore := countPersistedPairs(d)
	if persistedBefore == 0 {
		t.Fatal("expected at least one persisted pair before crash")
	}
	total := d.MatrixN * (d.MatrixN - 1) / 2

	// Second run resumes: progress must begin above zero and finish.
	svc2 := New(cfg)
	minDistanceDone := total
	var events []Progress
	res, err := svc2.Align(context.Background(), strings.NewReader(fastaInput), func(p Progress) {
		events = append(events, p)
		if p.Stage == stageDistance && p.Total > 0 && p.Done < minDistanceDone {
			minDistanceDone = p.Done
		}
	})
	if err != nil {
		t.Fatalf("resume failed: %v", err)
	}
	if minDistanceDone < persistedBefore {
		t.Fatalf("resume recomputed completed work: first seen done=%d, persisted=%d",
			minDistanceDone, persistedBefore)
	}
	if res.Alignment.Width == 0 {
		t.Fatal("alignment missing after resume")
	}
}

func countPersistedPairs(d jobData) int {
	n := d.MatrixN
	c := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if d.Matrix[i][j] == d.Matrix[i][j] { // not NaN
				c++
			} else if i < len(d.MatrixErrs) && d.MatrixErrs[i][j] != "" {
				c++
			}
		}
	}
	return c
}

func TestResumeAfterFullRunIsComplete(t *testing.T) {
	cfg, path := cfgWithCheckpoint(t)
	svc := New(cfg)
	res1, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil)
	if err != nil {
		t.Fatal(err)
	}

	// Second run loads a "done" checkpoint and must return the same
	// alignment without recomputation.
	svc2 := New(cfg)
	res2, err := svc2.Align(context.Background(), strings.NewReader(fastaInput), nil)
	if err != nil {
		t.Fatal(err)
	}
	for i, id := range res1.Alignment.IDs {
		if string(res1.Alignment.Rows[i]) != string(res2.Alignment.Rows[i]) {
			t.Fatalf("resumed alignment differs for %s", id)
		}
	}
	_ = path
}

func TestCorruptCheckpointRejected(t *testing.T) {
	cfg, path := cfgWithCheckpoint(t)
	svc := New(cfg)
	if _, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil); err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	// Flip a payload byte to invalidate the checksum.
	raw[len("MSACHECKPOINT-V1\n")+5] ^= 0xFF
	if err := os.WriteFile(path, raw, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil); err == nil {
		t.Fatal("corrupt checkpoint must be rejected")
	} else if !strings.Contains(err.Error(), "corrupt") {
		t.Fatalf("want corrupt error, got: %v", err)
	}
}

func TestCheckpointInputMismatchRejected(t *testing.T) {
	cfg, _ := cfgWithCheckpoint(t)
	svc := New(cfg)
	if _, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil); err != nil {
		t.Fatal(err)
	}
	different := strings.ReplaceAll(fastaInput, "TTTTAAAAGGGG", "GGGGTTTTAAAA")
	if _, err := svc.Align(context.Background(), strings.NewReader(different), nil); err == nil {
		t.Fatal("checkpoint for a different input must be rejected")
	}
}

func TestCheckpointConfigMismatchRejected(t *testing.T) {
	cfg, _ := cfgWithCheckpoint(t)
	svc := New(cfg)
	if _, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil); err != nil {
		t.Fatal(err)
	}
	cfg2 := cfg
	cfg2.Scoring.GapOpen = -99
	svc2 := New(cfg2)
	if _, err := svc2.Align(context.Background(), strings.NewReader(fastaInput), nil); err == nil {
		t.Fatal("checkpoint for different scoring must be rejected")
	}
}

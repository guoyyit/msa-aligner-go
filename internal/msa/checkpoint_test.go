package msa

import (
	"errors"
	"os"
	"path/filepath"
	"testing"
)

func sampleCheckpoint() *Checkpoint {
	return &Checkpoint{
		InputHash:   "abc123",
		KmerSize:    3,
		Concurrency: 4,
		GapOpen:     -10,
		GapExtend:   -1,
		MatrixName:  "BLOSUM62",
		SeqIDs:      []string{"a", "b"},
		SeqData:     []string{"ACGT", "ACGT"},
		DistanceDone: true,
		DistMatrix:  [][]float64{{0, 0.5}, {0.5, 0}},
		MergeDone:   true,
		Order:       &MergeOrder{N: 2, Events: []MergeEvent{{0, 1}}},
		NextMerge:   1,
		ClusterIDs:  map[int][]string{2: {"a", "b"}},
		ClusterSeq:  map[int][][]byte{2: {[]byte("ACGT"), []byte("ACGT")}},
	}
}

func TestCheckpointRoundTrip(t *testing.T) {
	path := filepath.Join(t.TempDir(), "ckpt.json")
	cp := sampleCheckpoint()
	if err := SaveCheckpoint(path, cp); err != nil {
		t.Fatal(err)
	}
	loaded, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	if loaded.InputHash != cp.InputHash || loaded.KmerSize != cp.KmerSize {
		t.Errorf("mismatch: %+v vs %+v", loaded, cp)
	}
	if loaded.NextMerge != cp.NextMerge {
		t.Errorf("NextMerge = %d, want %d", loaded.NextMerge, cp.NextMerge)
	}
	if loaded.Order == nil || loaded.Order.Events[0] != (MergeEvent{0, 1}) {
		t.Errorf("Order = %+v", loaded.Order)
	}
	if loaded.ClusterSeq[2][0][0] != 'A' {
		t.Errorf("ClusterSeq mismatch: %q", loaded.ClusterSeq[2][0])
	}
}

func TestCheckpointCorruptedPayload(t *testing.T) {
	path := filepath.Join(t.TempDir(), "ckpt.json")
	if err := SaveCheckpoint(path, sampleCheckpoint()); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	// Flip a byte in the payload region (after the checksum field).
	for i := len(data) - 1; i >= 0; i-- {
		if data[i] == 'A' || data[i] == 'C' || data[i] == 'G' || data[i] == 'T' {
			data[i] = 'X'
			break
		}
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	_, err = LoadCheckpoint(path)
	if !errors.Is(err, ErrCheckpointCorrupt) {
		t.Errorf("err = %v, want ErrCheckpointCorrupt", err)
	}
}

func TestCheckpointCorruptedChecksum(t *testing.T) {
	path := filepath.Join(t.TempDir(), "ckpt.json")
	if err := SaveCheckpoint(path, sampleCheckpoint()); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	// Flip the first hex char of the checksum.
	if data[14] == '0' {
		data[14] = '1'
	} else {
		data[14] = '0'
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	_, err = LoadCheckpoint(path)
	if !errors.Is(err, ErrCheckpointCorrupt) {
		t.Errorf("err = %v, want ErrCheckpointCorrupt", err)
	}
}

func TestConfigFromCheckpoint(t *testing.T) {
	cp := sampleCheckpoint()
	cfg := ConfigFromCheckpoint(cp)
	if cfg.Matrix == nil || cfg.Matrix.Name != "BLOSUM62" {
		t.Errorf("matrix = %+v", cfg.Matrix)
	}
	if cfg.GapOpen != -10 || cfg.GapExtend != -1 {
		t.Errorf("gaps = %v %v", cfg.GapOpen, cfg.GapExtend)
	}
}

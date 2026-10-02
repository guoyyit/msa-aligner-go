package msa

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"time"
)

// Checkpoint is the serializable resume state. The substitution matrix is
// stored by name (reconstructed on load) to keep the checkpoint compact.
type Checkpoint struct {
	Version     int
	InputHash   string
	KmerSize    int
	Concurrency int
	GapOpen     float64
	GapExtend   float64
	MatrixName  string

	SeqIDs  []string
	SeqData []string

	DistanceDone bool
	DistMatrix   [][]float64
	MergeDone    bool
	Order        *MergeOrder

	NextMerge  int
	ClusterIDs map[int][]string
	ClusterSeq map[int][][]byte

	UpdatedAt time.Time
}

// checkpointVersion is the current checkpoint format version.
const checkpointVersion = 1

// checkpointFile wraps the payload with a checksum for corruption detection.
type checkpointFile struct {
	Checksum string          `json:"checksum"`
	Payload  json.RawMessage `json:"payload"`
}

// SaveCheckpoint writes the checkpoint to path atomically with a checksum.
func SaveCheckpoint(path string, cp *Checkpoint) error {
	cp.Version = checkpointVersion
	cp.UpdatedAt = time.Now().UTC()
	payload, err := json.Marshal(cp)
	if err != nil {
		return fmt.Errorf("msa: checkpoint marshal: %w", err)
	}
	sum := sha256.Sum256(payload)
	wrapper := checkpointFile{
		Checksum: hex.EncodeToString(sum[:]),
		Payload:  payload,
	}
	data, err := json.Marshal(wrapper)
	if err != nil {
		return fmt.Errorf("msa: checkpoint marshal wrapper: %w", err)
	}

	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return fmt.Errorf("msa: checkpoint mkdir: %w", err)
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0o644); err != nil {
		return fmt.Errorf("msa: checkpoint write: %w", err)
	}
	if err := os.Rename(tmp, path); err != nil {
		return fmt.Errorf("msa: checkpoint rename: %w", err)
	}
	return nil
}

// LoadCheckpoint reads and verifies a checkpoint. A checksum mismatch returns
// ErrCheckpointCorrupt and never silently resumes.
func LoadCheckpoint(path string) (*Checkpoint, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("msa: checkpoint read: %w", err)
	}
	var wrapper checkpointFile
	if err := json.Unmarshal(data, &wrapper); err != nil {
		return nil, fmt.Errorf("%w: malformed file: %v", ErrCheckpointCorrupt, err)
	}
	sum := sha256.Sum256(wrapper.Payload)
	if hex.EncodeToString(sum[:]) != wrapper.Checksum {
		return nil, fmt.Errorf("%w: checksum mismatch", ErrCheckpointCorrupt)
	}
	var cp Checkpoint
	if err := json.Unmarshal(wrapper.Payload, &cp); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrCheckpointCorrupt, err)
	}
	if cp.Version != checkpointVersion {
		return nil, fmt.Errorf("%w: version %d", ErrCheckpointCorrupt, cp.Version)
	}
	return &cp, nil
}

// ConfigFromCheckpoint reconstructs a Config from checkpoint fields.
func ConfigFromCheckpoint(cp *Checkpoint) Config {
	var mat *SubstitutionMatrix
	switch cp.MatrixName {
	case "DNA":
		mat = DNAMatrix(2, -1)
	case "identity":
		mat = IdentityMatrix(5, -4)
	default:
		mat = BLOSUM62()
	}
	return Config{
		Matrix:      mat,
		GapOpen:     cp.GapOpen,
		GapExtend:   cp.GapExtend,
		KmerSize:    cp.KmerSize,
		Concurrency: cp.Concurrency,
	}
}

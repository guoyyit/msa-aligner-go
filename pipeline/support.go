package pipeline

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"time"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

// ErrLimitExceeded reports an input rejected before any heavy work ran.
type ErrLimitExceeded struct {
	What  string
	Value int64
	Limit int64
}

func (e *ErrLimitExceeded) Error() string {
	return fmt.Sprintf("pipeline: %s %d exceeds limit %d", e.What, e.Value, e.Limit)
}

// DefaultLimits returns the built-in resource ceilings.
func DefaultLimits() Limits {
	return Limits{
		MaxSequences:      5000,
		MaxSequenceLength: 200_000,
		MaxTotalResidues:  50_000_000,
		MaxMemoryBytes:    2 * 1024 * 1024 * 1024,
	}
}

// DefaultConfig returns the fully reproducible default configuration.
func DefaultConfig() Config {
	return Config{
		Alphabet:           fasta.AlphabetDNA,
		K:                  3,
		Concurrency:        1,
		Scoring:            align.DefaultParams(),
		Limits:             DefaultLimits(),
		CheckpointInterval: 5 * time.Second,
	}
}

// withDefaults fills zero-valued configuration fields with defaults.
func (c Config) withDefaults() Config {
	d := DefaultConfig()
	if c.Alphabet == "" {
		c.Alphabet = d.Alphabet
	}
	if c.K == 0 {
		c.K = d.K
	}
	if c.Concurrency == 0 {
		c.Concurrency = d.Concurrency
	}
	if c.Scoring.Subst == nil {
		c.Scoring = d.Scoring
	}
	if c.Scoring.GapOpen == 0 && c.Scoring.GapExtend == 0 {
		c.Scoring.GapOpen, c.Scoring.GapExtend = d.Scoring.GapOpen, d.Scoring.GapExtend
	}
	if c.Limits.MaxSequences == 0 {
		c.Limits = d.Limits
	}
	if c.CheckpointInterval == 0 {
		c.CheckpointInterval = d.CheckpointInterval
	}
	return c
}

// hashInput hashes the parsed sequences in their deterministic input order.
func hashInput(seqs []fasta.Sequence) string {
	h := sha256.New()
	for _, s := range seqs {
		fmt.Fprintf(h, "%s\x00%s\x00", s.ID, s.Data)
	}
	return hex.EncodeToString(h.Sum(nil))
}

// hashConfig captures every parameter that can change the output.
func hashConfig(cfg Config) (string, error) {
	// Reconstruct substitution matrix through a fixed enumeration so the
	// hash never depends on internal map iteration order.
	type substEntry struct {
		A, B  byte
		Score int
	}
	// Snapshot via Score over the supported alphabet, which is exhaustive
	// for the documented input alphabets.
	var letters []byte
	switch cfg.Alphabet {
	case fasta.AlphabetProtein:
		letters = []byte("ACDEFGHIKLMNPQRSTVWYBXZ*")
	default:
		letters = []byte("ACGTN")
		if cfg.Alphabet == fasta.AlphabetRNA {
			letters = []byte("ACGUN")
		}
	}
	entries := make([]substEntry, 0, len(letters)*len(letters))
	for _, a := range letters {
		for _, b := range letters {
			entries = append(entries, substEntry{a, b, cfg.Scoring.Subst.Score(a, b)})
		}
	}
	doc := struct {
		Alphabet string
		K        int
		GapOpen  int
		GapExt   int
		Subst    []substEntry
	}{cfg.Alphabet, cfg.K, cfg.Scoring.GapOpen, cfg.Scoring.GapExtend, entries}
	raw, err := json.Marshal(doc)
	if err != nil {
		return "", err
	}
	sum := sha256.Sum256(raw)
	return hex.EncodeToString(sum[:]), nil
}

// checkLimits rejects oversized inputs before distance/alignment work runs.
func checkLimits(seqs []fasta.Sequence, lim Limits) error {
	if len(seqs) > lim.MaxSequences {
		return &ErrLimitExceeded{What: "sequence count", Value: int64(len(seqs)), Limit: int64(lim.MaxSequences)}
	}
	var total int64
	for _, s := range seqs {
		if len(s.Data) > lim.MaxSequenceLength {
			return &ErrLimitExceeded{What: "sequence length (" + s.ID + ")", Value: int64(len(s.Data)), Limit: int64(lim.MaxSequenceLength)}
		}
		total += int64(len(s.Data))
	}
	if total > lim.MaxTotalResidues {
		return &ErrLimitExceeded{What: "total residues", Value: total, Limit: lim.MaxTotalResidues}
	}
	// Matrix estimate: N*N float64 cells. Alignment DP peaks near
	// N*width, which is bounded by MaxTotalResidues; check the matrix.
	n := int64(len(seqs))
	matrixBytes := n * n * 8
	if matrixBytes > lim.MaxMemoryBytes {
		return &ErrLimitExceeded{What: "distance matrix bytes", Value: matrixBytes, Limit: lim.MaxMemoryBytes}
	}
	return nil
}

// jobData is the resumable payload embedded in a checkpoint.
type jobData struct {
	Stage      string        `json:"stage"` // parsing|distance|cluster|align|done
	IDs        []string      `json:"ids"`
	Data       [][]byte      `json:"data"`
	Descs      []string      `json:"descriptions,omitempty"`
	MatrixN    int           `json:"matrix_n"`
	Matrix     [][]jsonFloat `json:"matrix"`
	MatrixErrs [][]string    `json:"matrix_errors,omitempty"`
	Steps      []stepJSON    `json:"steps,omitempty"`
	Members    [][]int       `json:"members,omitempty"`
	AlignedIDs []string      `json:"aligned_ids,omitempty"`
	Rows       []string      `json:"rows,omitempty"`
	Appended   []int         `json:"appended_columns,omitempty"`
}

// jsonFloat serializes NaN/Inf as JSON null instead of failing.
type jsonFloat float64

func (f jsonFloat) MarshalJSON() ([]byte, error) {
	v := float64(f)
	if v != v || v > 1e308 || v < -1e308 {
		return []byte("null"), nil
	}
	return json.Marshal(v)
}

func (f *jsonFloat) UnmarshalJSON(b []byte) error {
	if string(b) == "null" {
		*f = jsonFloat(math.NaN())
		return nil
	}
	var v float64
	if err := json.Unmarshal(b, &v); err != nil {
		return err
	}
	*f = jsonFloat(v)
	return nil
}

type stepJSON struct {
	Step     int     `json:"step"`
	Left     int     `json:"left"`
	Right    int     `json:"right"`
	Distance float64 `json:"distance"`
	NewIndex int     `json:"new_index"`
}

// stage names reported in progress events.
const (
	stageParsing  = "parsing"
	stageDistance = "distance"
	stageCluster  = "cluster"
	stageAlign    = "align"
	stageDone     = "done"
)

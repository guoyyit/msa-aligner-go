// Package checkpoint persists resumable MSA job state with an explicit
// SHA-256 integrity trailer. Corrupted, truncated, or tampered files are
// rejected with ErrCorrupt rather than silently recovered.
package checkpoint

import (
	"bytes"
	"crypto/sha256"
	"encoding/base64"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"time"
)

// ErrCorrupt reports a checkpoint that failed integrity verification.
var ErrCorrupt = errors.New("checkpoint corrupt")

// CurrentVersion is the only envelope version Load accepts.
const CurrentVersion = 1

const (
	magicLine = "MSACHECKPOINT-V1"
	sepLine   = "---"
)

// State is the versioned, integrity-protected checkpoint envelope.
type State struct {
	Version    int             `json:"version"`
	JobID      string          `json:"job_id"`
	InputHash  string          `json:"input_hash"`
	ConfigHash string          `json:"config_hash"`
	UpdatedAt  time.Time       `json:"updated_at"`
	Data       json.RawMessage `json:"data"`
}

// Checksum returns the hex SHA-256 over the exact bytes that must appear
// between the magic line and the checksum trailer.
func checksum(body []byte) string {
	sum := sha256.Sum256(body)
	return hex.EncodeToString(sum[:])
}

// Save atomically writes the state to path: payload is written to a
// sibling temp file, fsynced, and renamed over the destination.
func Save(path string, st *State) error {
	if st == nil {
		return errors.New("checkpoint: nil state")
	}
	st.Version = CurrentVersion
	payload, err := json.MarshalIndent(st, "", "  ")
	if err != nil {
		return fmt.Errorf("checkpoint: marshal: %w", err)
	}
	sum := checksum(payload)
	var buf bytes.Buffer
	buf.WriteString(magicLine)
	buf.WriteByte('\n')
	buf.WriteString(base64.StdEncoding.EncodeToString(payload))
	buf.WriteByte('\n')
	buf.WriteString(sepLine)
	buf.WriteByte('\n')
	buf.WriteString(sum)
	buf.WriteByte('\n')

	dir := filepath.Dir(path)
	if dir == "" {
		dir = "."
	}
	tmp, err := os.CreateTemp(dir, ".msa-checkpoint-*")
	if err != nil {
		return fmt.Errorf("checkpoint: create temp: %w", err)
	}
	tmpName := tmp.Name()
	cleanup := true
	defer func() {
		if cleanup {
			_ = os.Remove(tmpName)
		}
	}()
	if _, err := tmp.Write(buf.Bytes()); err != nil {
		_ = tmp.Close()
		return fmt.Errorf("checkpoint: write: %w", err)
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		return fmt.Errorf("checkpoint: sync: %w", err)
	}
	if err := tmp.Close(); err != nil {
		return fmt.Errorf("checkpoint: close: %w", err)
	}
	if err := os.Rename(tmpName, path); err != nil {
		return fmt.Errorf("checkpoint: rename: %w", err)
	}
	cleanup = false
	return nil
}

// Load reads, structurally verifies, and checksum-verifies a checkpoint.
func Load(path string) (*State, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("checkpoint: read: %w", err)
	}
	return Parse(raw)
}

// Parse verifies the serialized checkpoint form without touching the disk.
func Parse(raw []byte) (*State, error) {
	lines := bytes.Split(bytes.TrimRight(raw, "\n"), []byte("\n"))
	if len(lines) != 4 {
		return nil, fmt.Errorf("%w: expected 4 lines, got %d", ErrCorrupt, len(lines))
	}
	if string(lines[0]) != magicLine {
		return nil, fmt.Errorf("%w: bad magic", ErrCorrupt)
	}
	if string(lines[2]) != sepLine {
		return nil, fmt.Errorf("%w: missing separator", ErrCorrupt)
	}
	wantSum := string(bytes.TrimSpace(lines[3]))
	payloadEnc := bytes.TrimSpace(lines[1])
	payload, err := base64.StdEncoding.DecodeString(string(payloadEnc))
	if err != nil {
		return nil, fmt.Errorf("%w: payload not base64: %v", ErrCorrupt, err)
	}
	if got := checksum(payload); got != wantSum {
		return nil, fmt.Errorf("%w: checksum mismatch (want %s, got %s)", ErrCorrupt, wantSum, got)
	}
	var st State
	dec := json.NewDecoder(bytes.NewReader(payload))
	if err := dec.Decode(&st); err != nil {
		return nil, fmt.Errorf("%w: invalid json: %v", ErrCorrupt, err)
	}
	if st.Version != CurrentVersion {
		return nil, fmt.Errorf("%w: unsupported version %d", ErrCorrupt, st.Version)
	}
	return &st, nil
}

// DecodeData unmarshals the embedded job data into out.
func (s *State) DecodeData(out any) error {
	if len(s.Data) == 0 {
		return errors.New("checkpoint: empty data payload")
	}
	if err := json.Unmarshal(s.Data, out); err != nil {
		return fmt.Errorf("%w: data: %v", ErrCorrupt, err)
	}
	return nil
}

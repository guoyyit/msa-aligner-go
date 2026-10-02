package checkpoint

import (
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func sampleState() *State {
	return &State{
		JobID:      "job-1",
		InputHash:  "in",
		ConfigHash: "cfg",
		UpdatedAt:  time.Unix(1700000000, 0).UTC(),
		Data:       json.RawMessage(`{"stage":"distance","done":3}`),
	}
}

func TestSaveLoadRoundTrip(t *testing.T) {
	path := filepath.Join(t.TempDir(), "cp")
	if err := Save(path, sampleState()); err != nil {
		t.Fatal(err)
	}
	got, err := Load(path)
	if err != nil {
		t.Fatal(err)
	}
	if got.JobID != "job-1" || got.Version != CurrentVersion {
		t.Fatalf("bad state: %+v", got)
	}
	var data map[string]any
	if err := got.DecodeData(&data); err != nil {
		t.Fatal(err)
	}
	if data["stage"] != "distance" {
		t.Fatalf("data mismatch: %v", data)
	}
}

func TestLoadMissingFile(t *testing.T) {
	_, err := Load(filepath.Join(t.TempDir(), "nope"))
	if err == nil || errors.Is(err, ErrCorrupt) {
		t.Fatalf("missing file must be a plain read error, got %v", err)
	}
}

func TestCorruptionVariants(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "cp")
	if err := Save(path, sampleState()); err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}

	// Flip one payload byte: checksum must fail.
	tampered := append([]byte(nil), raw...)
	tampered[len(magicLine)+3] ^= 0xFF
	if _, err := Parse(tampered); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("payload tamper not detected: %v", err)
	}

	// Truncate the file.
	if _, err := Parse(raw[:len(raw)-5]); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("truncation not detected: %v", err)
	}

	// Corrupt the trailing checksum line.
	sumLine := raw[:len(raw)-1]
	sumLine[len(sumLine)-1] = '0'
	if _, err := Parse(sumLine); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("checksum tamper not detected: %v", err)
	}

	// Replace magic.
	badMagic := append([]byte(nil), raw...)
	copy(badMagic[:8], "XXXXXXXX")
	if _, err := Parse(badMagic); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("bad magic not detected: %v", err)
	}

	// Garbage content.
	if _, err := Parse([]byte("garbage")); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("garbage not detected: %v", err)
	}
}

func TestSaveRejectsNil(t *testing.T) {
	if err := Save(filepath.Join(t.TempDir(), "cp"), nil); err == nil {
		t.Fatal("nil state must be rejected")
	}
}

func TestSaveOverwriteAtomic(t *testing.T) {
	path := filepath.Join(t.TempDir(), "cp")
	st1 := sampleState()
	if err := Save(path, st1); err != nil {
		t.Fatal(err)
	}
	st2 := sampleState()
	st2.JobID = "job-2"
	if err := Save(path, st2); err != nil {
		t.Fatal(err)
	}
	got, err := Load(path)
	if err != nil {
		t.Fatal(err)
	}
	if got.JobID != "job-2" {
		t.Fatalf("overwrite failed: %+v", got)
	}
	if len(strings.TrimSpace(string(mustRead(t, path)))) == 0 {
		t.Fatal("empty checkpoint file")
	}
}

func mustRead(t *testing.T, p string) []byte {
	t.Helper()
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	return b
}

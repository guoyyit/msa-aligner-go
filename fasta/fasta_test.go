package fasta

import (
	"strings"
	"testing"
)

func parseMust(t *testing.T, input string) ([]Sequence, error) {
	t.Helper()
	return Parse(strings.NewReader(input))
}

func kindSet(err error) map[ErrorKind]bool {
	out := map[ErrorKind]bool{}
	for _, e := range err.(Errors) {
		out[e.Kind] = true
	}
	return out
}

func TestParseOK(t *testing.T) {
	seqs, err := parseMust(t, ">a some desc\nACGT\nacgt u\n>b\nTTTT\n")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) != 2 {
		t.Fatalf("got %d sequences", len(seqs))
	}
	if seqs[0].ID != "a" || seqs[0].Description != "some desc" {
		t.Fatalf("bad header: %+v", seqs[0])
	}
	if string(seqs[0].Data) != "ACGTACGTT" {
		t.Fatalf("normalization failed: %q", seqs[0].Data)
	}
	if string(seqs[1].Data) != "TTTT" {
		t.Fatalf("bad second seq: %q", seqs[1].Data)
	}
}

func TestParseEmptyInput(t *testing.T) {
	_, err := parseMust(t, "  \n\n")
	if err == nil {
		t.Fatal("expected error for empty input")
	}
	if !kindSet(err)[ErrKindMalformed] {
		t.Fatalf("expected malformed, got %v", err)
	}
}

func TestParseEmptySequence(t *testing.T) {
	_, err := parseMust(t, ">a\n>b\nACGT\n")
	if err == nil || !kindSet(err)[ErrKindEmptySequence] {
		t.Fatalf("expected empty_sequence, got %v", err)
	}
}

func TestParseIllegalCharacter(t *testing.T) {
	_, err := parseMust(t, ">a\nACGZ\n")
	if err == nil || !kindSet(err)[ErrKindIllegalCharacter] {
		t.Fatalf("expected illegal_character, got %v", err)
	}
}

func TestParseDuplicateID(t *testing.T) {
	_, err := parseMust(t, ">a\nACGT\n>a\nTTTT\n")
	if err == nil || !kindSet(err)[ErrKindDuplicateID] {
		t.Fatalf("expected duplicate_id, got %v", err)
	}
}

func TestParseMissingID(t *testing.T) {
	_, err := parseMust(t, ">\nACGT\n")
	if err == nil || !kindSet(err)[ErrKindMissingID] {
		t.Fatalf("expected missing_id, got %v", err)
	}
}

func TestParseDataBeforeHeader(t *testing.T) {
	_, err := parseMust(t, "ACGT\n>a\nACGT\n")
	if err == nil || !kindSet(err)[ErrKindMalformed] {
		t.Fatalf("expected malformed, got %v", err)
	}
}

func TestParseProteinAlphabet(t *testing.T) {
	p := Parser{Alphabet: AlphabetProtein}
	seqs, err := p.Parse(strings.NewReader(">x\nMAK*\n"))
	if err != nil {
		t.Fatalf("protein parse failed: %v", err)
	}
	if string(seqs[0].Data) != "MAK*" {
		t.Fatalf("bad data: %q", seqs[0].Data)
	}
	if _, err := p.Parse(strings.NewReader(">x\nMAK1\n")); err == nil {
		t.Fatal("expected illegal character in protein alphabet")
	}
}

func TestParsePartialReturnOnErrors(t *testing.T) {
	seqs, err := parseMust(t, ">good\nACGT\n>bad\nZZZ\n")
	if err == nil {
		t.Fatal("expected error")
	}
	if len(seqs) != 1 || seqs[0].ID != "good" {
		t.Fatalf("expected valid records returned alongside error, got %+v", seqs)
	}
}

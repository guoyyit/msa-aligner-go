package fasta

import (
	"errors"
	"testing"
)

func TestParseValid(t *testing.T) {
	data := []byte(">seq1 description here\nACGT acgt\n>seq2\nGGGTTT\n")
	seqs, err := Parse(data)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(seqs) != 2 {
		t.Fatalf("got %d sequences, want 2", len(seqs))
	}
	if seqs[0].ID != "seq1" || seqs[0].Desc != "description here" {
		t.Errorf("seq1 header = %q %q", seqs[0].ID, seqs[0].Desc)
	}
	if seqs[0].Seq != "ACGTACGT" {
		t.Errorf("seq1 = %q, want ACGTACGT", seqs[0].Seq)
	}
	if seqs[1].ID != "seq2" || seqs[1].Seq != "GGGTTT" {
		t.Errorf("seq2 = %+v", seqs[1])
	}
}

func TestParseEmptyInput(t *testing.T) {
	if _, err := Parse(nil); !errors.Is(err, ErrEmptyInput) {
		t.Errorf("err = %v, want ErrEmptyInput", err)
	}
	if _, err := Parse([]byte("\n\n")); !errors.Is(err, ErrEmptyInput) {
		t.Errorf("err = %v, want ErrEmptyInput", err)
	}
}

func TestParseEmptySequence(t *testing.T) {
	data := []byte(">a\nACGT\n>b\n>c\nGGG\n")
	_, err := Parse(data)
	if !errors.Is(err, ErrEmptySequence) {
		t.Errorf("err = %v, want ErrEmptySequence", err)
	}
}

func TestParseInvalidCharacter(t *testing.T) {
	data := []byte(">a\nACGTZ\n>b\nACGX1\n")
	_, err := Parse(data)
	if !errors.Is(err, ErrInvalidCharacter) {
		t.Errorf("err = %v, want ErrInvalidCharacter", err)
	}
}

func TestParseDuplicateID(t *testing.T) {
	data := []byte(">a\nACGT\n>a\nGGGG\n")
	_, err := Parse(data)
	if !errors.Is(err, ErrDuplicateID) {
		t.Errorf("err = %v, want ErrDuplicateID", err)
	}
}

func TestParseEmptyHeader(t *testing.T) {
	data := []byte(">\nACGT\n")
	_, err := Parse(data)
	if !errors.Is(err, ErrEmptyHeader) {
		t.Errorf("err = %v, want ErrEmptyHeader", err)
	}
}

func TestParseSeqWithoutHeader(t *testing.T) {
	data := []byte("ACGT\n>a\nGGGG\n")
	_, err := Parse(data)
	if !errors.Is(err, ErrSeqWithoutHeader) {
		t.Errorf("err = %v, want ErrSeqWithoutHeader", err)
	}
}

func TestParseLowercaseAndWhitespace(t *testing.T) {
	data := []byte(">a\n  acgt \n>b\na c g t\n")
	seqs, err := Parse(data)
	if err != nil {
		t.Fatalf("unexpected: %v", err)
	}
	if seqs[0].Seq != "ACGT" || seqs[1].Seq != "ACGT" {
		t.Errorf("got %q %q", seqs[0].Seq, seqs[1].Seq)
	}
}

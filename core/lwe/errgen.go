package lwe

import (
	"math/rand"
)

type ErrGen struct {
	stdev float64
	// bound int64
}

func NewErrorGenerator(stdev float64) *ErrGen {
	return &ErrGen{stdev: stdev}
}

func (erg *ErrGen) GenErr() int64 {
	e := int64(rand.NormFloat64() * erg.stdev)
	return e
}

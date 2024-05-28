package ngsw

import "github.com/tuneinsight/lattigo/v5/core/rlwe"

type Ciphertext struct {
	Value rlwe.GadgetCiphertext
} // Only use one part of

// LevelQ returns the level of the modulus Q of the target.
func (ct Ciphertext) LevelQ() int {
	return ct.Value.LevelQ()
}

// LevelP returns the level of the modulus P of the target.
func (ct Ciphertext) LevelP() int {
	return ct.Value.LevelP()
}

func NewCiphertext(params rlwe.ParameterProvider, levelQ, levelP, BaseTwoDecomposition int) (ct *Ciphertext) {
	return &Ciphertext{
		Value: *rlwe.NewGadgetCiphertext(params, 1, levelQ, levelP, BaseTwoDecomposition),
	}
}

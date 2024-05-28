package ngsw

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

// Automorphism computes phi(ct), where phi is the map X -> X^galEl. The method requires
// that the corresponding RotationKey has been added to the Evaluator. The method will
// return an error if either ctIn or opOut degree is not equal to 1.
func (eval Evaluator) Automorphism(ctIn *rlwe.Ciphertext, galEl uint64, opOut *rlwe.Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("cannot apply Automorphism: input and output Ciphertext must be of degree 1")
	}

	if galEl == 1 {
		if opOut != ctIn {
			opOut.Copy(ctIn)
		}
		return
	}

	level := utils.Min(ctIn.Level(), opOut.Level())
	ringQ := eval.GetRLWEParameters().RingQ().AtLevel(level)

	ctTmp := &rlwe.Ciphertext{Element: rlwe.Element[ring.Poly]{Value: []ring.Poly{eval.BuffQP[0].Q, eval.BuffQP[1].Q}}}
	ctTmp.MetaData = ctIn.MetaData

	// first do the automorphism on the encrypted polynomial, should be given in NTT form
	if ctIn.IsNTT {
		opOut.Value[0].Zero() // ringQ.AutomorphismNTTWithIndex(ctIn.Value[0], eval.automorphismIndex[galEl], opOut.Value[0])
		ringQ.AutomorphismNTTWithIndex(ctIn.Value[1], eval.AutomorphismIndex(galEl), ctTmp.Value[1])
	} else {
		return fmt.Errorf("Automorphism Cipher should be in NTT form")
	}

	// get the key-switching key
	var evk *rlwe.GaloisKey
	if evk, err = eval.CheckAndGetGaloisKey(galEl); err != nil {
		return fmt.Errorf("cannot apply Automorphism: %w", err)
	}

	// multiply (key switch to original key)
	eval.ExternalProduct(ctTmp, &Ciphertext{Value: evk.GadgetCiphertext}, opOut)

	return
}

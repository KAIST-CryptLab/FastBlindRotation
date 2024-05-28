package ngsw

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
)

type KeyGenerator struct {
	*Encryptor
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as EvaluationKeys.
func NewKeyGenerator(params rlwe.Parameters) *KeyGenerator {
	return &KeyGenerator{
		Encryptor: NewEncryptor(params, nil),
	}
}

// generating automorphism keys

/*

 */

/*
func (kgen KeyGenerator) GenAutoKeysNew(galEls []uint64, sk *SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (aks []*AutoKey) {
	levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(kgen.params, evkParams)
	aks = make([]*AutoKey, len(galEls))
	for i, galEl := range galEls {
		aks[i] = newAutoKey(kgen.params, levelQ, levelP, BaseTwoDecomposition)
		kgen.GenAutoKey(galEl, sk, aks[i])
	}
	return
}
*/

func newAutoKey(params rlwe.Parameters, levelQ, levelP, BaseTwoDecomposition int) *rlwe.GaloisKey {
	return &rlwe.GaloisKey{
		EvaluationKey: rlwe.EvaluationKey{GadgetCiphertext: *rlwe.NewGadgetCiphertext(params, 1, levelQ, levelP, BaseTwoDecomposition)},
		NthRoot:       params.GetRLWEParameters().RingQ().NthRoot(),
	}
}

func (kgen KeyGenerator) GenAutoKey(k uint64, sk *rlwe.SecretKey, sk_inv *rlwe.SecretKey, ak *rlwe.GaloisKey) {

	skIn := sk.Value
	skOut := kgen.buffQP

	ringQP := kgen.Encryptor.GetRLWEParameters().RingQP().AtLevel(ak.LevelQ(), ak.LevelP())
	ringQ := ringQP.RingQ

	// get f(X^k)
	index, err := ring.AutomorphismNTTIndex(ringQ.N(), ringQ.NthRoot(), k)
	if err != nil {
		panic(err)
	}

	ringQ.AutomorphismNTTWithIndex(skIn.Q, index, skOut.Q)
	ringQ.MulCoeffsMontgomery(skOut.Q, sk_inv.Value.Q, skOut.Q) // f(X^k)/f(X)

	kgen.genEvaluationKey(skOut.Q, sk_inv.Value, &ak.EvaluationKey)
	ak.GaloisElement = k
	ak.NthRoot = ringQ.NthRoot()
}

func (kgen KeyGenerator) GenAutoKeysNew(ks []uint64, sk *rlwe.SecretKey, sk_inv *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (aks []*rlwe.GaloisKey) {
	param := *kgen.Encryptor.GetRLWEParameters()
	levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(param, evkParams)
	aks = make([]*rlwe.GaloisKey, len(ks))

	for i, k := range ks {
		aks[i] = newAutoKey(param, levelQ, levelP, BaseTwoDecomposition)
		kgen.GenAutoKey(k, sk, sk_inv, aks[i])
	}

	return
}

func (kgen KeyGenerator) genEvaluationKey(skIn ring.Poly, skOut ringqp.Poly, evk *rlwe.EvaluationKey) {

	enc := kgen.WithKey(&rlwe.SecretKey{Value: skOut})

	// Samples an encryption of zero for each element of the EvaluationKey.
	for i := 0; i < len(evk.Value); i++ {
		for j := 0; j < len(evk.Value[i]); j++ {
			if err := enc.EncryptZero(rlwe.Element[ringqp.Poly]{MetaData: &rlwe.MetaData{CiphertextMetaData: rlwe.CiphertextMetaData{IsNTT: true, IsMontgomery: true}}, Value: []ringqp.Poly(evk.Value[i][j])}); err != nil {
				// Sanity check, this error should not happen.
				panic(err)
			}
		}
	}

	// Adds the plaintext (input-key) to the EvaluationKey.
	if err := AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, &Ciphertext{Value: evk.GadgetCiphertext}, *kgen.Encryptor.GetRLWEParameters().RingQP(), kgen.buffQP.Q); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}
}

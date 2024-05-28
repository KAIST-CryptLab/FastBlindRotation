package lwe

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v5/ring"
)

type Decryptor struct {
	lweParam Parameters
	skCoeffs []*big.Int
}

func NewDecryptor(lweParam Parameters, ringQ *ring.Ring, sk ring.Poly) *Decryptor {
	skCoeffs := make([]*big.Int, ringQ.N())
	for i := range skCoeffs {
		skCoeffs[i] = new(big.Int)
	}
	ringQ.AtLevel(0).PolyToBigintCentered(sk, 1, skCoeffs)
	return &Decryptor{
		lweParam: lweParam,
		skCoeffs: skCoeffs,
	}
}

func (dec Decryptor) DecryptNew(ct *Ciphertext) uint64 {

	lweParam := dec.lweParam
	n := lweParam.N()
	q := lweParam.Q()

	var dotProd int64
	dotProd = 0
	for i := 0; i < int(n); i++ {
		dotProd += int64(ct.A[i]) * dec.skCoeffs[i].Int64()
	}

	// calculate m = (b - a * s)
	m := int64(ct.B) - dotProd
	m = ((m % int64(q)) + int64(q)) % int64(q)

	return uint64(m)
}

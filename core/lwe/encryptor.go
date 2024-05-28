package lwe

import (
	"math/big"
	"math/rand"

	"github.com/tuneinsight/lattigo/v5/ring"
)

type Encryptor struct {
	lweParam Parameters
	skCoeffs []*big.Int
}

func NewEncryptor(lweParam Parameters, ringQ *ring.Ring, sk ring.Poly) *Encryptor {
	skCoeffs := make([]*big.Int, ringQ.N())
	for i := range skCoeffs {
		skCoeffs[i] = new(big.Int)
	}
	ringQ.AtLevel(0).PolyToBigintCentered(sk, 1, skCoeffs)

	return &Encryptor{
		lweParam: lweParam,
		skCoeffs: skCoeffs,
	}
}

// m includes error term
func (enc Encryptor) EncryptNew(m uint64) (ct *Ciphertext) {

	lweParam := enc.lweParam
	q := lweParam.Q()
	n := lweParam.N()

	// generate new coefficient for encrypting
	a := make([]uint64, n)
	for i := range a {
		a[i] = uint64(rand.Intn(int(q))) // (uint64(rand.Intn(int(q/2)))*2 + 1) % q
	}

	// calculate b = a * s + m + e
	var dotProd int64
	dotProd = 0
	for i := 0; i < int(n); i++ {
		si := enc.skCoeffs[i].Int64()
		if si != 0 {
			dotProd += int64(a[i]) * si
		}
	}

	dotProd = ((dotProd % int64(q)) + int64(q)) % int64(q)
	b := uint64(dotProd) + m
	b %= q

	return &Ciphertext{A: a, B: b}
}

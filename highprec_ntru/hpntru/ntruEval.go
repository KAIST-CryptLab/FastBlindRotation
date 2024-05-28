package hpntru

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/ngsw"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
)

type Evaluator struct {
	*ngsw.Evaluator
	paramsBR  rlwe.Parameters
	paramsLWE lwe.Parameters

	accumulator *rlwe.Ciphertext
	auto_buffer *rlwe.Ciphertext
}

func InitTestPolynomial(p, q uint64, ringQ *ring.Ring) (F ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.ModuliChain()[:ringQ.Level()+1]

	N := ringQ.N()

	for j, qi := range Q {
		for i := 0; i < (N >> 1); i++ {
			F.Coeffs[j][2*i] = uint64(math.Round(float64(qi)/float64(2*p))) * (uint64(i) / (q / (2 * p)))
		}
	}

	ringQ.NTT(F, F)

	return
}

func InitTestPolynomialDebug(p, q uint64, ringQ *ring.Ring) (F ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.ModuliChain()[:ringQ.Level()+1]

	for j := range Q {
		F.Coeffs[j][2*96] = 8888 // uint64(math.Round(float64(qi)/float64(2*p))) * (uint64(i) / (q / (2 * p)))
	}

	ringQ.NTT(F, F)

	return
}

// NewEvaluator instantiates a new Evaluator.
func NewEvaluator(paramsBR rlwe.Parameters, paramsLWE lwe.Parameters, evk *rlwe.MemEvaluationKeySet) (eval *Evaluator) {
	eval = new(Evaluator)
	eval.Evaluator = ngsw.NewEvaluator(paramsBR, evk)
	eval.paramsBR = paramsBR
	eval.paramsLWE = paramsLWE

	// assure that q = N
	if paramsLWE.Q() != uint64(paramsBR.N()) {
		panic("q != N")
	}

	eval.accumulator = rlwe.NewCiphertext(paramsBR, 1, eval.paramsBR.MaxLevel())
	eval.auto_buffer = rlwe.NewCiphertext(paramsBR, 1, eval.paramsBR.MaxLevel())

	return
}

// returns NTRU cipher
// test poly should come as r(X^2) format (q = N)
func (eval *Evaluator) BlindRotate(ct *lwe.Ciphertext, testPoly *ring.Poly, BRK BlindRotationEvaluationKeySet) (res *rlwe.Ciphertext, err error) {
	var a = ct.GetA()
	var b = ct.GetB()

	var n = len(a)
	var N = eval.paramsBR.N()
	var twoN = N << 1
	twoN_big := new(big.Int).SetUint64(uint64(twoN))

	w := make([]uint64, n)
	w_inv := make([]uint64, n+1)

	for i := 0; i < n; i++ {
		w[i] = a[i]*2 + 1
		wInv_big := new(big.Int).SetUint64(w[i])
		wInv_big.ModInverse(wInv_big, twoN_big)
		w_inv[i] = uint64(wInv_big.Int64())
	}
	w_inv[n] = 1

	acc := eval.accumulator
	// buff := eval.auto_buffer

	ringQ := eval.paramsBR.RingQ()
	ringQ.AutomorphismNTT(*testPoly, w_inv[0], acc.Value[1])
	ringQ.INTT(acc.Value[1], acc.Value[1])

	k := -2 * int(b) * int(w_inv[0])
	k = ((k % twoN) + twoN) % twoN
	acc.Value[1] = MulByMonomial(ringQ, acc.Value[1], k)

	ringQ.NTT(acc.Value[1], acc.Value[1])

	/*
		// automorphism
		var brk_i *ngsw.Ciphertext
		if brk_i, err = BRK.GetBlindRotationKey(0); err != nil {
			panic(err)
		}
		eval.ExternalProduct(acc, brk_i, acc)

		exp := (w[0] * w_inv[1]) % uint64(2*N)
		if exp != 1 {
			eval.Automorphism(acc, exp, acc)
		}
	*/
	var brk_i *ngsw.Ciphertext
	for i := 0; i < n; i++ {
		if brk_i, err = BRK.GetBlindRotationKey(i); err != nil {
			panic(err)
		}
		eval.ExternalProduct(acc, brk_i, acc) // evk_i
		exp := (w[i] * w_inv[i+1]) % uint64(2*N)
		if k != 1 {
			eval.Automorphism(acc, exp, acc)
		}
	}
	if brk_i, err = BRK.GetBlindRotationKey(n); err != nil {
		panic(err)
	}
	eval.ExternalProduct(acc, brk_i, acc) // evk_n

	res = acc.CopyNew()
	return
}

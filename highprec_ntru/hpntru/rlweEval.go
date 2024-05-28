package hpntru

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/rgsw"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
)

func UNUSED_T(x ...interface{}) {}

type RLWEEvaluator struct {
	*rgsw.Evaluator
	paramsBR  rlwe.Parameters
	paramsLWE lwe.Parameters

	accumulator *rlwe.Ciphertext
	auto_buffer *rlwe.Ciphertext

	galoisGenDiscreteLog map[uint64]int
}

// NewEvaluator instantiates a new Evaluator.
func NewRLWEEvaluator(paramsBR rlwe.Parameters, paramsLWE lwe.Parameters) (eval *RLWEEvaluator) {
	eval = new(RLWEEvaluator)
	eval.Evaluator = rgsw.NewEvaluator(paramsBR, nil)
	eval.paramsBR = paramsBR
	eval.paramsLWE = paramsLWE

	// assure that q = N
	if paramsLWE.Q() != uint64(paramsBR.N()) {
		panic("q != N")
	}

	eval.accumulator = rlwe.NewCiphertext(paramsBR, 1, eval.paramsBR.MaxLevel())
	eval.auto_buffer = rlwe.NewCiphertext(paramsBR, 1, eval.paramsBR.MaxLevel())

	eval.galoisGenDiscreteLog = getGaloisElementInverseMap(ring.GaloisGen, eval.paramsBR.N())

	return
}

// returns NTRU cipher
// test poly should come as r(X^2) format (q = N)
func (eval *RLWEEvaluator) BlindRotate(ct *lwe.Ciphertext, testPoly *ring.Poly, BRK RLWEBlindRotationEvaluationKeySet, extkey *rgsw.Ciphertext) (res *rlwe.Ciphertext, err error) {
	evk, err := BRK.GetEvaluationKeySet()
	if err != nil {
		panic(err)
	}
	eval.Evaluator = eval.Evaluator.WithKey(evk)

	var a = ct.GetA()
	var b = ct.GetB()

	var n = len(a)
	var N = eval.paramsBR.N()
	var twoN = uint64(N << 1)

	w := make([]uint64, n)
	for i := 0; i < n; i++ {
		w[i] = (a[i]*2 + 1) % twoN
	}

	acc := eval.accumulator
	buff := eval.auto_buffer

	ringQ := eval.paramsBR.RingQ()

	ringQ.INTT(*testPoly, acc.Value[1])
	acc.Value[1] = MulByMonomial(ringQ, acc.Value[1], int(twoN)-2*int(b))
	ringQ.NTT(acc.Value[1], acc.Value[1])
	ringQ.AutomorphismNTT(acc.Value[1], ringQ.NthRoot()-ring.GaloisGen, acc.Value[0])
	acc.Value[1].Zero()

	if err = eval.BlindRotateCore(w, acc, BRK); err != nil {
		return nil, fmt.Errorf("BlindRotateCore: %s", err)
	}

	eval.ExternalProduct(acc, extkey, acc)
	res = acc.CopyNew()

	UNUSED_T(b, buff)
	return
}

// BlindRotateCore implements Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *RLWEEvaluator) BlindRotateCore(a []uint64, acc *rlwe.Ciphertext, BRK RLWEBlindRotationEvaluationKeySet) (err error) {

	evk, err := BRK.GetEvaluationKeySet()

	if err != nil {
		return err
	}

	eval.Evaluator = eval.Evaluator.WithKey(evk)

	// GaloisElement(k) = GaloisGen^{k} mod 2N
	GaloisElement := eval.paramsBR.GaloisElement

	// Maps a[i] to (+/-) g^{k} mod 2N
	discreteLogSets := eval.getDiscreteLogSets(a, eval.paramsBR.N())

	Nhalf := eval.paramsBR.N() >> 1
	twoN := eval.paramsBR.N() << 1

	// Algorithm 3 of https://eprint.iacr.org/2022/198
	var v int = 0
	// Lines 3 to 9 (negative set of a[i] = -g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, -i, v, acc, BRK); err != nil {
			return
		}
	}

	// Line 10 (0 in the negative set is 2N)
	if _, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, twoN, 0, acc, BRK); err != nil {
		return
	}

	// Line 12
	// acc = acc(X^{-g})
	if err = eval.Automorphism(acc, eval.paramsBR.RingQ().NthRoot()-ring.GaloisGen, acc); err != nil {
		return
	}

	// Lines 13 - 19 (positive set of a[i] = g^{k} mod 2N)
	for i := Nhalf - 1; i > 0; i-- {
		if v, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, i, v, acc, BRK); err != nil {
			return
		}
	}

	// Lines 20 - 21 (0 in the positive set is 0)
	if _, err = eval.evaluateFromDiscreteLogSets(GaloisElement, discreteLogSets, 0, 0, acc, BRK); err != nil {
		return
	}

	return
}

// evaluateFromDiscreteLogSets loops of Algorithm 3 of https://eprint.iacr.org/2022/198
func (eval *RLWEEvaluator) evaluateFromDiscreteLogSets(GaloisElement func(k int) (galEl uint64), sets map[int][]int, k, v int, acc *rlwe.Ciphertext, BRK RLWEBlindRotationEvaluationKeySet) (int, error) {

	// Checks if k is in the discrete log sets
	if set, ok := sets[k]; ok {
		// First condition of line 7 or 17
		if v != 0 {

			if err := eval.Automorphism(acc, GaloisElement(v), acc); err != nil {
				return v, err
			}

			v = 0
		}

		for _, j := range set {

			brk, err := BRK.GetBlindRotationKey(j)
			if err != nil {
				return v, err
			}

			// acc = acc * RGSW(X^{s[j]})
			eval.ExternalProduct(acc, brk, acc)
		}
	}

	v++

	// Second and third conditions of line 7 or 17
	if v == windowSize || k == 1 || k == -1 {

		if err := eval.Automorphism(acc, GaloisElement(v), acc); err != nil {
			return v, err
		}

		v = 0
	}

	return v, nil
}

// getGaloisElementInverseMap generates a map [(+/-) g^{k} mod 2N] = +/- k
func getGaloisElementInverseMap(GaloisGen uint64, N int) (GaloisGenDiscreteLog map[uint64]int) {

	twoN := N << 1
	NHalf := N >> 1
	mask := uint64(twoN - 1)

	GaloisGenDiscreteLog = map[uint64]int{}

	var pow uint64 = 1
	for i := 0; i < NHalf; i++ {
		GaloisGenDiscreteLog[pow] = i
		GaloisGenDiscreteLog[uint64(twoN)-pow] = -i
		pow *= GaloisGen
		pow &= mask
	}

	return
}

// getDiscreteLogSets returns map[+/-k] = [i...] for a[0 <= i < N] = {(+/-) g^{k} mod 2N for +/- k}
func (eval *RLWEEvaluator) getDiscreteLogSets(a []uint64, N int) (discreteLogSets map[int][]int) {

	GaloisGenDiscreteLog := eval.galoisGenDiscreteLog
	twoN := N << 1

	// Maps (2*N*a[i]/QLWE) to -N/2 < k <= N/2 for a[i] = (+/- 1) * g^{k}
	discreteLogSets = map[int][]int{}
	for i, ai := range a {

		if ai&1 != 1 && ai != 0 {
			panic("getDiscreteLogSets: a[i] is not odd and thus not an element of Z_{2N}^{*} -> a[i] = (+/- 1) * g^{k} does not exist.")
		}

		dlog := GaloisGenDiscreteLog[ai]
		if ai == uint64(twoN-1) {
			dlog = twoN
		}

		if _, ok := discreteLogSets[dlog]; !ok {
			discreteLogSets[dlog] = []int{i}
		} else {
			discreteLogSets[dlog] = append(discreteLogSets[dlog], i)
		}
	}

	return
}

/*

func (eval *RLWEEvaluator) BlindRotate(ct *lwe.Ciphertext, testPoly *ring.Poly, BRK RLWEBlindRotationEvaluationKeySet, ) (res *rlwe.Ciphertext, err error) {
	var a = ct.GetA()
	var b = ct.GetB()

	var n = len(a)
	var N = eval.paramsBR.N()
	var twoN = uint64(N << 1)
	twoN_big := new(big.Int).SetUint64(twoN)

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
	buff := eval.auto_buffer

	ringQ := eval.paramsBR.RingQ()
	ringQ.AutomorphismNTT(*testPoly, w_inv[0], acc.Value[1])
	ringQ.INTT(acc.Value[1], acc.Value[1])

	k := (4*N*N - 2*int(b)*int(w_inv[0])) % (2 * N)
	acc.Value[1] = MulByMonomial(ringQ, acc.Value[1], k)

	ringQ.NTT(acc.Value[1], acc.Value[1])
	ringQ.MForm(acc.Value[1], acc.Value[1])

	// automorphism
	var brk_i *rgsw.Ciphertext
	for i := 0; i < n; i++ {
		if brk_i, err = BRK.GetBlindRotationKey(i); err != nil {
			panic(err)
		}
		eval.ExternalProduct(acc, brk_i, acc) // evk_i
		k := (w[i] * w_inv[i+1]) % uint64(2*N)
		if k != 1 {
			eval.Automorphism(acc, k, buff)
			acc.Copy(buff)
		}
	}
	if brk_i, err = BRK.GetBlindRotationKey(n); err != nil {
		panic(err)
	}
	eval.ExternalProduct(acc, brk_i, acc) // evk_n

	res = acc.CopyNew()
	return
}

func (eval *Evaluator) BlindRotate(ct *lwe.Ciphertext, testPoly *ring.Poly, BRK BlindRotationEvaluationKeySet, extkey *rgsw.Ciphertext) (res *rlwe.Ciphertext, err error) {
	acc := eval.accumulator

	a := ct.A
	b := ct.B

	var n = len(a)
	var N = eval.paramsBR.N()
	var twoN = uint64(N << 1)
	w := make([]uint64, n)

	for i := 0; i < n; i++ {
		w[i] = (a[i]*2 + 1) % twoN
	}

	// Line 2 of Algorithm 7 of https://eprint.iacr.org/2022/198
	// Acc = (f(X^{-g}) * X^{-g * b}, 0)

	ringQBR := eval.paramsBR.RingQ().AtLevel(0)

	ringQBR.INTT(*testPoly, acc.Value[1])
	acc.Value[1] = MulByMonomial(ringQBR, acc.Value[1], -2*int(b))
	ringQBR.NTT(acc.Value[1], acc.Value[1])
	ringQBR.MForm(acc.Value[1], acc.Value[1])
	ringQBR.AutomorphismNTT(acc.Value[1], ringQBR.NthRoot()-ring.GaloisGen, acc.Value[0])
	acc.Value[1].Zero()

	if err = eval.BlindRotateCore(w, acc, BRK); err != nil {
		return nil, fmt.Errorf("BlindRotateCore: %s", err)
	}

	eval.ExternalProduct(acc, extkey, acc)


	res = acc.CopyNew()
	return
}
*/

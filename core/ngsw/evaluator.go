package ngsw

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils"
)

type Evaluator struct {
	rlwe.Evaluator
}

func NewEvaluator(params rlwe.ParameterProvider, evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*rlwe.NewEvaluator(params, evk)}
}

// ShallowCopy creates a shallow copy of this Evaluator in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Evaluators can be used concurrently.
func (eval Evaluator) ShallowCopy() *Evaluator {
	return &Evaluator{*eval.Evaluator.ShallowCopy()}
}

// WithKey creates a shallow copy of the receiver Evaluator for which the new EvaluationKey is evaluationKey
// and where the temporary buffers are shared. The receiver and the returned Evaluators cannot be used concurrently.
func (eval Evaluator) WithKey(evk rlwe.EvaluationKeySet) *Evaluator {
	return &Evaluator{*eval.Evaluator.WithKey(evk)}
}

func (eval Evaluator) ExternalProduct(op0 *rlwe.Ciphertext, op1 *Ciphertext, opOut *rlwe.Ciphertext) {

	levelQ, levelP := op1.LevelQ(), op1.LevelP()

	var c0QP, c1QP ringqp.Poly
	if op0 == opOut {
		c0QP, c1QP = eval.BuffQP[1], eval.BuffQP[2]
	} else {
		c0QP, c1QP = ringqp.Poly{Q: opOut.Value[0], P: eval.BuffQP[1].P}, ringqp.Poly{Q: opOut.Value[1], P: eval.BuffQP[2].P}
	}

	if levelP < 1 {
		params := eval.GetRLWEParameters()
		if ringQ := params.RingQ(); levelQ == 0 && levelP == -1 {
			eval.externalProductInPlaceNoPSingleQ(op0, op1, c1QP.Q)
			ringQ.AtLevel(0).IMForm(c1QP.Q, opOut.Value[1])
		} else {
			eval.externalProductInPlaceSinglePAndBitDecomp(op0, op1, c1QP)
			if levelP == 0 {
				eval.BasisExtender.ModDownQPtoQNTT(levelQ, levelP, c1QP.Q, c1QP.P, opOut.Value[1])
			} else {
				opOut.Value[0].CopyLvl(levelQ, c0QP.Q)
				opOut.Value[1].CopyLvl(levelQ, c1QP.Q)
			}
		}
	} else {
		fmt.Printf("Implement Later\n")
	}
}

func (eval Evaluator) externalProductInPlaceNoPSingleQ_temp(ct0 *rlwe.Ciphertext, ngsw *Ciphertext, c1 ring.Poly) {
	params := eval.GetRLWEParameters()
	ringQ := params.RingQ().AtLevel(0)
	subRing := ringQ.SubRings[0]
	pw2 := ngsw.Value.BaseTwoDecomposition
	mask := uint64(((1 << pw2) - 1))

	cw := eval.BuffQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	acc1 := c1.Coeffs[0]

	// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
	// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])

	ringQ.INTTLazy(ct0.Value[1], eval.BuffInvNTT)

	// d := len(ngsw.Value.Value[0])
	// decomp_res := decompose_centered(eval.BuffInvNTT.Coeffs[0], uint64(pw2), d, ringQ.Modulus().Int64())

	for j := range ngsw.Value.Value[0] {
		ring.MaskVec(eval.BuffInvNTT.Coeffs[0], j*pw2, mask, cw)
		if j == 0 {
			subRing.NTTLazy(cw, cwNTT)
			subRing.MulCoeffsLazy(ngsw.Value.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
		} else {
			subRing.NTTLazy(cw, cwNTT)
			subRing.MulCoeffsLazyThenAddLazy(ngsw.Value.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
		}
	}
}

func (eval Evaluator) externalProductInPlaceNoPSingleQ(ct0 *rlwe.Ciphertext, ngsw *Ciphertext, c1 ring.Poly) {
	params := eval.GetRLWEParameters()
	ringQ := params.RingQ().AtLevel(0)
	subRing := ringQ.SubRings[0]
	pw2 := ngsw.Value.BaseTwoDecomposition
	// mask := uint64(((1 << pw2) - 1))

	// cw := eval.BuffQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	acc1 := c1.Coeffs[0]

	// (a, b) + (c0 * rgsw[0][0], c0 * rgsw[0][1])
	// (a, b) + (c1 * rgsw[1][0], c1 * rgsw[1][1])

	ringQ.INTTLazy(ct0.Value[1], eval.BuffInvNTT)

	d := len(ngsw.Value.Value[0])
	decomp_res := decompose_centered(eval.BuffInvNTT.Coeffs[0], uint64(pw2), d, ringQ.Modulus().Int64())

	for j := range ngsw.Value.Value[0] {
		if j == 0 {
			subRing.NTTLazy(decomp_res[j], cwNTT)
			subRing.MulCoeffsLazy(ngsw.Value.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
		} else {
			subRing.NTTLazy(decomp_res[j], cwNTT)
			subRing.MulCoeffsLazyThenAddLazy(ngsw.Value.Value[0][j][1].Q.Coeffs[0], cwNTT, acc1)
		}
	}
}

func decompose_centered(coeffs []uint64, base_log uint64, d int, Q int64) [][]uint64 {
	base := 1 << base_log
	digit_mask := uint64(base - 1)
	base_over_2_threshold := int64(1 << (base_log - 1))
	carry := 0

	var unsigned_digit uint64
	var signed_digit int64

	N := len(coeffs)
	decomp_res := make([][]uint64, d)
	for i := range decomp_res {
		decomp_res[i] = make([]uint64, N)
	}

	for i, v := range coeffs {
		for j := 0; j < d; j++ {
			unsigned_digit = (v >> (uint64(j) * base_log)) & digit_mask
			if carry == 1 {
				unsigned_digit += uint64(carry)
				carry = 0
			}

			signed_digit = int64(unsigned_digit)
			if signed_digit >= base_over_2_threshold {
				signed_digit -= int64(base)
				carry = 1
			}

			decomp_res[j][i] = uint64((signed_digit + Q) % Q)
		}
	}

	return decomp_res
}

func (eval Evaluator) externalProductInPlaceSinglePAndBitDecomp(ct0 *rlwe.Ciphertext, ntruv *Ciphertext, c1QP ringqp.Poly) {

	// rgsw = [(-as + P*w*m1 + e, a), (-bs + e, b + P*w*m1)]
	// ct = [0, g/f + u/f]
	// opOut = [<ct, rgsw[0]>, <ct, rgsw[1]>] = [ct[0] * rgsw[0][0] + ct[1] * rgsw[0][1], ct[0] * rgsw[1][0] + ct[1] * rgsw[1][1]]
	levelQ := ntruv.LevelQ()
	levelP := ntruv.LevelP()

	params := eval.GetRLWEParameters()
	ringQP := params.RingQP().AtLevel(levelQ, levelP)

	ringQ := ringQP.RingQ
	ringP := ringQP.RingP

	pw2 := ntruv.Value.BaseTwoDecomposition
	mask := uint64(((1 << pw2) - 1))
	if mask == 0 {
		mask = 0xFFFFFFFFFFFFFFFF
	}

	BaseRNSDecompositionVectorSize := ntruv.Value.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := ntruv.Value.BaseTwoDecompositionVectorSize()

	ringQ.INTT(ct0.Value[1], eval.BuffInvNTT) // BuffInvNTT = Mont(g/f + u/f)
	cw := eval.BuffQP[0].Q.Coeffs[0]
	cwNTT := eval.BuffBitDecomp

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {
			// TODO: center values if mask == 0
			ring.MaskVec(eval.BuffInvNTT.Coeffs[i], j*pw2, mask, cw)

			if i == 0 && j == 0 {
				for u, s := range ringQ.SubRings[:levelQ+1] {
					s.NTTLazy(cw, cwNTT)
					s.MulCoeffsMontgomery(ntruv.Value.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u])
				}

				if ringP != nil {
					for u, s := range ringP.SubRings[:levelP+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomery(ntruv.Value.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u])
					}
				}

			} else {
				for u, s := range ringQ.SubRings[:levelQ+1] {
					s.NTTLazy(cw, cwNTT)
					s.MulCoeffsMontgomeryThenAdd(ntruv.Value.Value[i][j][1].Q.Coeffs[u], cwNTT, c1QP.Q.Coeffs[u])
				}

				if ringP != nil {
					for u, s := range ringP.SubRings[:levelP+1] {
						s.NTTLazy(cw, cwNTT)
						s.MulCoeffsMontgomeryThenAdd(ntruv.Value.Value[i][j][1].P.Coeffs[u], cwNTT, c1QP.P.Coeffs[u])
					}
				}
			}
		}
	}
}

// use the algorithm in 2020-015
func (eval Evaluator) CombineTestPoly(cts map[int]*rlwe.Ciphertext, n, N int) *rlwe.Ciphertext {
	BR_cts := NewQueue() // make it in queue structure
	ord := BitRevGen(uint(n))
	ringQ := eval.GetRLWEParameters().RingQ()
	invN := new(big.Int).ModInverse(new(big.Int).SetUint64(uint64(N)), ringQ.ModulusAtLevel[0])

	for _, v := range ord {
		ct_temp := cts[int(v)].CopyNew()
		ringQ.IMForm(ct_temp.Value[1], ct_temp.Value[1])
		ringQ.INTT(ct_temp.Value[1], ct_temp.Value[1])
		ringQ.MulScalarBigint(ct_temp.Value[1], invN, ct_temp.Value[1])
		ringQ.NTT(ct_temp.Value[1], ct_temp.Value[1])
		ringQ.MForm(ct_temp.Value[1], ct_temp.Value[1])
		BR_cts.Push(ct_temp)
	}

	var cur_n int
	var idx int
	var cur_len int

	idx = 1
	xPow2 := GenXPow2(ringQ, eval.GetRLWEParameters().LogN(), false)

	for BR_cts.Len() > 1 {
		cur_len = BR_cts.Len()
		cur_n = 1 << idx

		for i := 0; i < cur_len; i += 2 {
			// pop 2 elements
			ct_even := BR_cts.Pop().(*rlwe.Ciphertext)
			ct_odd := BR_cts.Pop().(*rlwe.Ciphertext)

			// mult ct_odd by X^(N/cur_n)
			ringQ.MulCoeffsMontgomery(ct_odd.Value[1], xPow2[len(xPow2)-idx], ct_odd.Value[1])

			ct1 := eval.BuffCt
			ringQ.Add(ct_even.Value[1], ct_odd.Value[1], ct1.Value[1])     // ct_even + X^N/n_cur * ct_odd
			ringQ.Sub(ct_even.Value[1], ct_odd.Value[1], ct_even.Value[1]) // ct_even - X^N/n_cur * ct_odd

			// ringQ.IMForm(ct_even.Value[1], ct_even.Value[1])
			eval.Automorphism(ct_even, uint64(cur_n+1), ct_odd)
			// ringQ.MForm(ct_odd.Value[1], ct_odd.Value[1])
			ringQ.Add(ct_odd.Value[1], ct1.Value[1], ct_odd.Value[1])

			// push the element to BR_cts
			BR_cts.Push(ct_odd)
		}

		idx += 1
	}

	fin := BR_cts.Pop().(*rlwe.Ciphertext)
	eval.Trace(fin, int(math.Logb(float64(n)))-1, fin)

	// After trace, mult by X^(2^k), add, repeat
	for i := 1; i < int(math.Log2(float64(N/n))); i++ {
		ringQ.MulCoeffsMontgomeryThenAdd(fin.Value[1], xPow2[i], fin.Value[1])
	}

	// ringQ.IMForm(fin.Value[1], fin.Value[1])
	return fin
}

func GenXPow2(r *ring.Ring, logN int, div bool) (xPow []ring.Poly) {

	// Compute X^{-n} from 0 to LogN
	xPow = make([]ring.Poly, logN)

	moduli := r.ModuliChain()[:r.Level()+1]
	BRC := r.BRedConstants()

	var idx int
	for i := 0; i < logN; i++ {

		idx = 1 << i

		if div {
			idx = r.N() - idx
		}

		xPow[i] = r.NewPoly()

		if i == 0 {

			for j := range moduli {
				xPow[i].Coeffs[j][idx] = ring.MForm(1, moduli[j], BRC[j])
			}

			r.NTT(xPow[i], xPow[i])

		} else {
			r.MulCoeffsMontgomery(xPow[i-1], xPow[i-1], xPow[i]) // X^{n} = X^{1} * X^{n-1}
		}
	}

	if div {
		r.Neg(xPow[0], xPow[0])
	}

	return
}

func (eval Evaluator) Trace(ctIn *rlwe.Ciphertext, logN int, opOut *rlwe.Ciphertext) (err error) {

	if ctIn.Degree() != 1 || opOut.Degree() != 1 {
		return fmt.Errorf("ctIn.Degree() != 1 or opOut.Degree() != 1")
	}

	params := eval.GetRLWEParameters()

	level := utils.Min(ctIn.Level(), opOut.Level())

	opOut.Resize(opOut.Degree(), level)

	*opOut.MetaData = *ctIn.MetaData

	gap := 1 << (params.LogN() - logN - 1)

	if logN == 0 {
		gap <<= 1
	}

	if gap > 1 {

		ringQ := params.RingQ().AtLevel(level)

		if ringQ.Type() == ring.ConjugateInvariant {
			gap >>= 1 // We skip the last step that applies phi(5^{-1})
		}

		// NInv := new(big.Int).SetUint64(uint64(gap))
		// NInv.ModInverse(NInv, ringQ.ModulusAtLevel[level])

		// pre-multiplication by (N/n)^-1
		// ringQ.MulScalarBigint(ctIn.Value[0], NInv, opOut.Value[0])
		// ringQ.MulScalarBigint(ctIn.Value[1], NInv, opOut.Value[1])

		if !ctIn.IsNTT {
			ringQ.NTT(opOut.Value[0], opOut.Value[0])
			ringQ.NTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = true
		}

		buff, err := rlwe.NewCiphertextAtLevelFromPoly(level, []ring.Poly{eval.BuffQP[3].Q, eval.BuffQP[4].Q})

		// Sanity check, this error should not happen unless the
		// evaluator's buffer thave been improperly tempered with.
		if err != nil {
			panic(err)
		}

		buff.IsNTT = true

		for i := logN; i < params.LogN()-1; i++ {

			if err = eval.Automorphism(opOut, params.GaloisElement(1<<i), buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if logN == 0 && ringQ.Type() == ring.Standard {

			if err = eval.Automorphism(opOut, ringQ.NthRoot()-1, buff); err != nil {
				return err
			}

			ringQ.Add(opOut.Value[0], buff.Value[0], opOut.Value[0])
			ringQ.Add(opOut.Value[1], buff.Value[1], opOut.Value[1])
		}

		if !ctIn.IsNTT {
			ringQ.INTT(opOut.Value[0], opOut.Value[0])
			ringQ.INTT(opOut.Value[1], opOut.Value[1])
			opOut.IsNTT = false
		}

	} else {
		if ctIn != opOut {
			opOut.Copy(ctIn)
		}
	}

	return
}

// Package bootstrapping implement the bootstrapping for the CKKS scheme.
package bootstrapping

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// Bootstrap re-encrypts a ciphertext to a ciphertext at MaxLevel - k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller than Q[0]/MessageRatio
// (it can't be equal since Q[0] is not a power of two).
// The message ratio is an optional field in the bootstrapping parameters, by default it set to 2^{LogMessageRatio = 8}.
// See the bootstrapping parameters for more information about the message ratio or other parameters related to the bootstrapping.
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) Bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	// Pre-processing
	ctDiff := ctIn.CopyNew()

	// Drops the level to 1
	for ctDiff.Level() > 1 {
		btp.DropLevel(ctDiff, 1)
	}

	// Brings the ciphertext scale to Q0/MessageRatio
	if ctDiff.Level() == 1 {

		// If one level is available, then uses it to match the scale
		btp.SetScale(ctDiff, rlwe.NewScale(btp.q0OverMessageRatio))

		// Then drops to level 0
		for ctDiff.Level() != 0 {
			btp.DropLevel(ctDiff, 1)
		}

	} else {

		// Does an integer constant mult by round((Q0/Delta_m)/ctscale)
		if scale := ctDiff.PlaintextScale.Float64(); scale != math.Exp2(math.Round(math.Log2(scale))) || btp.q0OverMessageRatio < scale {
			msgRatio := btp.EvalModParameters.LogMessageRatio
			panic(fmt.Sprintf("ciphertext scale must be a power of two smaller than Q[0]/2^{LogMessageRatio=%d} = %f but is %f", msgRatio, float64(btp.params.Q()[0])/math.Exp2(float64(msgRatio)), scale))
		}

		btp.ScaleUp(ctDiff, rlwe.NewScale(math.Round(btp.q0OverMessageRatio/ctDiff.PlaintextScale.Float64())), ctDiff)
	}

	// Scales the message to Q0/|m|, which is the maximum possible before ModRaise to avoid plaintext overflow.
	if scale := math.Round((float64(btp.params.Q()[0]) / btp.evalModPoly.MessageRatio()) / ctDiff.PlaintextScale.Float64()); scale > 1 {
		btp.ScaleUp(ctDiff, rlwe.NewScale(scale), ctDiff)
	}

	// 2^d * M + 2^(d-n) * e
	ctOut = btp.bootstrap(ctDiff.CopyNew())

	for i := 1; i < btp.Iterations; i++ {
		// 2^(d-n)*e <- [2^d * M + 2^(d-n) * e] - [2^d * M]
		tmp := btp.SubNew(ctDiff, ctOut)

		// 2^d * e
		btp.Mul(tmp, 1<<16, tmp)

		// 2^d * e + 2^(d-n) * e'
		tmp = btp.bootstrap(tmp)

		// 2^(d-n) * e + 2^(d-2n) * e'
		btp.Mul(tmp, 1/float64(uint64(1<<16)), tmp)

		if err := btp.Rescale(tmp, btp.params.PlaintextScale(), tmp); err != nil {
			panic(err)
		}

		// [2^d * M + 2^(d-2n) * e'] <- [2^d * M + 2^(d-n) * e] - [2^(d-n) * e + 2^(d-2n) * e']
		btp.Add(ctOut, tmp, ctOut)
	}

	return
}

func (btp *Bootstrapper) bootstrap(ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {

	// Step 1 : Extend the basis from q to Q
	ctOut = btp.modUpFromQ0(ctIn)

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (btp.evalModPoly.ScalingFactor().Float64() / btp.evalModPoly.MessageRatio()) / ctOut.PlaintextScale.Float64(); scale > 1 {
		btp.ScaleUp(ctOut, rlwe.NewScale(scale), ctOut)
	}

	//SubSum X -> (N/dslots) * Y^dslots
	btp.Trace(ctOut, ctOut.PlaintextLogDimensions[1], ctOut)

	// Step 2 : CoeffsToSlots (Homomorphic encoding)
	ctReal, ctImag := btp.CoeffsToSlotsNew(ctOut, btp.ctsMatrices)

	// Step 3 : EvalMod (Homomorphic modular reduction)
	// ctReal = Ecd(real)
	// ctImag = Ecd(imag)
	// If n < N/2 then ctReal = Ecd(real|imag)
	ctReal = btp.EvalModNew(ctReal, btp.evalModPoly)
	ctReal.PlaintextScale = btp.params.PlaintextScale()

	if ctImag != nil {
		ctImag = btp.EvalModNew(ctImag, btp.evalModPoly)
		ctImag.PlaintextScale = btp.params.PlaintextScale()
	}

	// Step 4 : SlotsToCoeffs (Homomorphic decoding)
	ctOut = btp.SlotsToCoeffsNew(ctReal, ctImag, btp.stcMatrices)

	return
}

func (btp *Bootstrapper) modUpFromQ0(ct *rlwe.Ciphertext) *rlwe.Ciphertext {

	if btp.EvkDtS != nil {
		btp.ApplyEvaluationKey(ct, btp.EvkDtS, ct)
	}

	ringQ := btp.params.RingQ().AtLevel(ct.Level())
	ringP := btp.params.RingP()

	for i := range ct.Value {
		ringQ.INTT(ct.Value[i], ct.Value[i])
	}

	// Extend the ciphertext with zero polynomials.
	ct.Resize(ct.Degree(), btp.params.MaxLevel())

	levelQ := btp.params.QCount() - 1
	levelP := btp.params.PCount() - 1

	ringQ = ringQ.AtLevel(levelQ)

	Q := ringQ.ModuliChain()
	P := ringP.ModuliChain()
	q := Q[0]
	BRCQ := ringQ.BRedConstants()
	BRCP := ringP.BRedConstants()

	var coeff, tmp, pos, neg uint64

	N := ringQ.N()

	// ModUp q->Q for ct[0] centered around q
	for j := 0; j < N; j++ {

		coeff = ct.Value[0].Coeffs[0][j]
		pos, neg = 1, 0
		if coeff >= (q >> 1) {
			coeff = q - coeff
			pos, neg = 0, 1
		}

		for i := 1; i < levelQ+1; i++ {
			tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
			ct.Value[0].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
		}
	}

	if btp.EvkStD != nil {

		ks := btp.Evaluator.Evaluator

		// ModUp q->QP for ct[1] centered around q
		for j := 0; j < N; j++ {

			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff > (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 0; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
				ks.BuffDecompQP[0].Q.Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg

			}

			for i := 0; i < levelP+1; i++ {
				tmp = ring.BRedAdd(coeff, P[i], BRCP[i])
				ks.BuffDecompQP[0].P.Coeffs[i][j] = tmp*pos + (P[i]-tmp)*neg
			}
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringQ.NTT(ks.BuffDecompQP[0].Q, ks.BuffDecompQP[i].Q)
		}

		for i := len(ks.BuffDecompQP) - 1; i >= 0; i-- {
			ringP.NTT(ks.BuffDecompQP[0].P, ks.BuffDecompQP[i].P)
		}

		ringQ.NTT(ct.Value[0], ct.Value[0])

		ctTmp := &rlwe.Ciphertext{}
		ctTmp.Value = []ring.Poly{ks.BuffQP[1].Q, ct.Value[1]}
		ctTmp.MetaData = ct.MetaData

		ks.GadgetProductHoisted(levelQ, ks.BuffDecompQP, &btp.EvkStD.GadgetCiphertext, ctTmp)
		ringQ.Add(ct.Value[0], ctTmp.Value[0], ct.Value[0])

	} else {

		for j := 0; j < N; j++ {

			coeff = ct.Value[1].Coeffs[0][j]
			pos, neg = 1, 0
			if coeff >= (q >> 1) {
				coeff = q - coeff
				pos, neg = 0, 1
			}

			for i := 1; i < levelQ+1; i++ {
				tmp = ring.BRedAdd(coeff, Q[i], BRCQ[i])
				ct.Value[1].Coeffs[i][j] = tmp*pos + (Q[i]-tmp)*neg
			}
		}

		ringQ.NTT(ct.Value[0], ct.Value[0])
		ringQ.NTT(ct.Value[1], ct.Value[1])
	}

	return ct
}

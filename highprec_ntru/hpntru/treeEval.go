package hpntru

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/ngsw"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type TreeEvaluator struct {
	*ngsw.Evaluator
	paramsBR  rlwe.Parameters
	paramsLWE lwe.Parameters

	accumulator *rlwe.Ciphertext
	auto_buffer *rlwe.Ciphertext
}

func InitReLUPolynomial(p, q uint64, ringQ *ring.Ring) (F_high ring.Poly, F_low map[int]ring.Poly) {
	F_high = ringQ.NewPoly()
	F_low = make(map[int]ring.Poly)

	Q := ringQ.ModuliChain()[:ringQ.Level()+1]

	N := ringQ.N()

	// F_high : 0 <= x < p/2 -> x, p/2 <= x -> 0
	// F_low[k] : 0 <= x < p/2 -> k, p/2 <= x -> 0 (0 <= k < p)
	for k := 0; k < int(p); k++ {
		F_low[k] = ringQ.NewPoly()
		for j := range Q {
			for i := 0; i < (N >> 1); i++ {
				m := (uint64(i) / (q / (2 * p)))
				if m < p/2 {
					F_low[k].Coeffs[j][2*i] = uint64(k) // uint64(math.Round(float64(qi)/float64(2*p))) * uint64(k)
				} else {
					F_low[k].Coeffs[j][2*i] = 0
				}
			}
		}
		ringQ.NTT(F_low[k], F_low[k])
	}

	for j := range Q {
		for i := 0; i < (N >> 1); i++ {
			m := (uint64(i) / (q / (2 * p)))
			if m < p/2 {
				F_high.Coeffs[j][2*i] = uint64(m) // uint64(math.Round(float64(qi)/float64(2*p))) * uint64(m)
				// F_high.Coeffs[j][2*i] = uint64(math.Round(float64(qi)/float64(2*p))) * uint64(m)
			} else {
				F_high.Coeffs[j][2*i] = 0
			}
		}
	}
	ringQ.NTT(F_high, F_high)

	return
}

// (Q+1)/2 * (1 + X^2 + ... X^(N-2)) * Q/2p
func InitMVPolynomial(p uint64, ringQ *ring.Ring) (F ring.Poly) {
	F = ringQ.NewPoly()
	Qbig := ringQ.Modulus()
	Q := ringQ.ModuliChain()[:ringQ.Level()+1]

	N := ringQ.N()
	twoInvbig := bignum.NewInt(2)
	twoInvbig.ModInverse(twoInvbig, Qbig)
	twoInv := twoInvbig.Uint64()

	for j, qi := range Q {
		for i := 0; i < (N >> 1); i++ {
			F.Coeffs[j][2*i] = (twoInv * uint64(math.Round(float64(qi)/float64(2*p)))) % qi
		}
	}

	ringQ.NTT(F, F)

	return
}

func InitMVDebug(p uint64, ringQ *ring.Ring) (F ring.Poly) {
	F = ringQ.NewPoly()
	Q := ringQ.ModuliChain()[:ringQ.Level()+1]

	for j := range Q {
		F.Coeffs[j][2*96] = 87777 // uint64(math.Round(float64(qi)/float64(2*p))) * (uint64(i) / (q / (2 * p)))
	}

	ringQ.NTT(F, F)

	return
}

var norm_squares []float64

// assume that f is given in NTT format, will be used in MVBR
// f only has even degree coefficients
func ConvertPoly(f ring.Poly, ringQ *ring.Ring) (f_conv ring.Poly) {
	f_copy := ringQ.NewPoly()
	f_conv = ringQ.NewPoly()
	f_copy.Copy(f)
	ringQ.INTT(f_copy, f_copy)

	Q := int64(ringQ.ModuliChain()[0])
	N := ringQ.N()
	for i := 0; i < (N >> 1); i++ {
		if i == 0 {
			f_conv.Coeffs[0][0] = f_copy.Coeffs[0][0] + f_copy.Coeffs[0][N-2]
		} else {
			sub := int64(f_copy.Coeffs[0][2*i]) - int64(f_copy.Coeffs[0][2*i-2])
			f_conv.Coeffs[0][2*i] = uint64(((sub % Q) + Q) % Q)
		}
	}

	var sum float64 = 0
	TVCoeffs := make([]*big.Int, N)
	for i := range TVCoeffs {
		TVCoeffs[i] = new(big.Int)
	}
	ringQ.PolyToBigintCentered(f_conv, 1, TVCoeffs)
	for i := range TVCoeffs {
		sum += float64(TVCoeffs[i].Int64()) * float64(TVCoeffs[i].Int64())
	}
	if true {
		// fmt.Printf("TV norm square : %f\n", sum)
		norm_squares = append(norm_squares, sum)
	}

	ringQ.NTT(f_conv, f_conv)
	return
}

func ConvertPoly2(f ring.Poly, ringQ *ring.Ring) (f_conv ring.Poly) {
	f_copy := ringQ.NewPoly()
	f_conv = ringQ.NewPoly()
	f_copy.Copy(f)
	ringQ.INTT(f_copy, f_copy)

	Q := int64(ringQ.ModuliChain()[0])
	N := ringQ.N()
	for i := 0; i < (N >> 1); i++ {
		if i == 0 {
			f_conv.Coeffs[0][0] = f_copy.Coeffs[0][0] + f_copy.Coeffs[0][N-2]
		} else {
			sub := int64(f_copy.Coeffs[0][2*i]) - int64(f_copy.Coeffs[0][2*i-2])
			f_conv.Coeffs[0][2*i] = uint64(((sub % Q) + Q) % Q)
		}
	}

	ringQ.NTT(f_conv, f_conv)
	return
}

// NewEvaluator instantiates a new Evaluator.
func NewTreeEvaluator(paramsBR rlwe.Parameters, paramsLWE lwe.Parameters, evk *rlwe.MemEvaluationKeySet) (eval *TreeEvaluator) {
	eval = new(TreeEvaluator)
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

func (eval *TreeEvaluator) MultiValueBR(ct *lwe.Ciphertext, testPoly *ring.Poly, testMultiPoly map[int]ring.Poly, BRK BlindRotationEvaluationKeySet) (res *rlwe.Ciphertext, multires map[int]*rlwe.Ciphertext, err error) {
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

	/*
		wmap := map[int][]int{}
		for i := 0; i < n; i++ {

			w_val := int(w[i])
			if _, ok := wmap[w_val]; !ok {
				wmap[w_val] = []int{i}
			} else {
				wmap[w_val] = append(wmap[w_val], i)
			}
		}

		for i, v := range wmap {
			fmt.Printf("%d, %v\n", i, v)
		}
	*/

	// Get the index mapping

	acc := eval.accumulator
	// buff := eval.auto_buffer

	ringQ := eval.paramsBR.RingQ()
	basePoly := InitMVPolynomial(eval.paramsLWE.P(), ringQ)
	convPoly := ConvertPoly2(*testPoly, ringQ)

	ringQ.AutomorphismNTT(basePoly, w_inv[0], acc.Value[1])
	ringQ.INTT(acc.Value[1], acc.Value[1])

	k := -2 * int(b) * int(w_inv[0])
	k = ((k % twoN) + twoN) % twoN
	acc.Value[1] = MulByMonomial(ringQ, acc.Value[1], k)

	ringQ.NTT(acc.Value[1], acc.Value[1])

	var brk_i *ngsw.Ciphertext
	for i := 0; i < n; i++ {
		if brk_i, err = BRK.GetBlindRotationKey(i); err != nil {
			panic(err)
		}
		eval.ExternalProduct(acc, brk_i, acc) // evk_i
		exp := (w[i] * w_inv[i+1]) % uint64(2*N)
		if exp != 1 {
			eval.Automorphism(acc, exp, acc)
		}
	}
	if brk_i, err = BRK.GetBlindRotationKey(n); err != nil {
		panic(err)
	}
	eval.ExternalProduct(acc, brk_i, acc) // evk_n

	res = acc.CopyNew()
	ringQ.MForm(convPoly, convPoly)
	ringQ.MulCoeffsMontgomery(convPoly, res.Value[1], res.Value[1])
	res.Value[0].Zero()

	if testMultiPoly != nil {
		multires = make(map[int]*rlwe.Ciphertext)
		for i, v := range testMultiPoly {
			multires[i] = acc.CopyNew()
			convPoly := ConvertPoly(v, ringQ)
			ringQ.MForm(convPoly, convPoly)
			ringQ.MulCoeffsMontgomery(convPoly, multires[i].Value[1], multires[i].Value[1])
		}
	}

	return
}

func (eval *TreeEvaluator) EncBlindRotate(ct *lwe.Ciphertext, testPoly *rlwe.Ciphertext, BRK BlindRotationEvaluationKeySet) (res *rlwe.Ciphertext, err error) {
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

	eval.Automorphism(testPoly, w_inv[0], acc)
	ringQ.INTT(acc.Value[1], acc.Value[1])

	k := (4*N*N - 2*int(b)*int(w_inv[0])) % (2 * N)
	acc.Value[1] = MulByMonomial(ringQ, acc.Value[1], k)

	ringQ.NTT(acc.Value[1], acc.Value[1])

	// automorphism
	var brk_i *ngsw.Ciphertext
	for i := 0; i < n; i++ {
		if i == 0 {
			if brk_i, err = BRK.GetBlindRotationKey(n + 1); err != nil {
				panic(err)
			}
		} else {
			if brk_i, err = BRK.GetBlindRotationKey(i); err != nil {
				panic(err)
			}
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

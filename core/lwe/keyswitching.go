package lwe

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/structs"
)

// define structure of LKSK (union of LWE ciphers)
type LKSK struct {
	KeySet map[int]structs.Matrix[*Ciphertext]
}

// f -> skLWE keyswitching key generation
// encryptor is based on skLWE, QKS

// need fixing
func LWE_KSKG(BKS int, f *rlwe.SecretKey, enc *Encryptor, ringQBR *ring.Ring, errgen *ErrGen) *LKSK {
	var KeySet map[int]structs.Matrix[*Ciphertext] = make(map[int]structs.Matrix[*Ciphertext])

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)
	fCoeffs := make([]*big.Int, ringQBR.N())
	for i := range fCoeffs {
		fCoeffs[i] = new(big.Int)
	}
	ringQBR.PolyToBigintCentered(fCopy, 1, fCoeffs)

	QKS := int(enc.lweParam.Q())
	dKS := int(math.Ceil(math.Log2(float64(QKS))/math.Log2(float64(BKS)))) + 1 // len(Decomp(uint64(QKS)-1, uint64(BKS), uint64(QKS)))

	var cur_val int
	for i := 0; i < ringQBR.N(); i++ {
		KeySet[i] = make(structs.Matrix[*Ciphertext], dKS)

		cur_val = int(fCoeffs[i].Int64())
		for j := 0; j < dKS; j++ { // dks
			KeySet[i][j] = make([]*Ciphertext, BKS-1)

			for v := 1; v < BKS/2; v++ {
				e := int(errgen.GenErr())
				m := ((cur_val*v+e)%QKS + QKS) % QKS // positive value modulo QKS (vB^jf_i + e)
				KeySet[i][j][v-1] = enc.EncryptNew(uint64(m))
			}
			for v := BKS / 2; v < BKS; v++ {
				e := int(errgen.GenErr())
				m := ((cur_val*(QKS-(v-BKS/2+1))+e)%QKS + QKS) % QKS // positive value modulo QKS (vB^jf_i + e)
				KeySet[i][j][v-1] = enc.EncryptNew(uint64(m))
			}

			cur_val *= BKS
		}
	}

	return &LKSK{KeySet: KeySet}
}

func LWE_KS(ct *rlwe.Ciphertext, BKS int, enc *Encryptor, ringQBR *ring.Ring, lksk *LKSK) (ct_ks *Ciphertext) {
	// Extract, ModSwitch, Decompose on the fly

	n := int(enc.lweParam.N())
	ct_ks = NewCiphertext(n)

	c_Poly := ringQBR.NewPoly()
	c_Poly.Copy(ct.Value[1])
	ringQBR.INTT(c_Poly, c_Poly)
	// ringQBR.IMForm(c_Poly, c_Poly)

	QKS := enc.lweParam.Q()
	QKSBig := bignum.NewInt(QKS)
	QBig := ringQBR.ModulusAtLevel[0]

	Coeffs := make([]*big.Int, ringQBR.N())
	for i := range Coeffs {
		Coeffs[i] = new(big.Int)
	}
	ringQBR.PolyToBigint(c_Poly, 1, Coeffs)

	N := ringQBR.N()
	A := make([]uint64, N)
	dKS := int(math.Ceil(math.Log2(float64(QKS))/math.Log2(float64(BKS)))) + 1

	for i := 0; i < N; i++ {
		Coeffs[i].Mul(Coeffs[i], QKSBig)
		bignum.DivRound(Coeffs[i], QBig, Coeffs[i])
	}

	for i := 0; i < N; i++ {
		if i == 0 {
			A[i] = Coeffs[i].Uint64()
		} else {
			A[i] = QKS - Coeffs[N-i].Uint64()
		}

		decomp := Decomp(A[i], uint64(BKS), QKS, dKS)

		for j, v := range decomp {
			if v != 0 {
				if v < uint64(BKS)/2 {
					tmp := lksk.KeySet[i][j][v-1].A
					for w := 0; w < n; w++ {
						ct_ks.A[w] = ring.CRed(ct_ks.A[w]+tmp[w], QKS)
					}
					ct_ks.B = ring.CRed(ct_ks.B+lksk.KeySet[i][j][v-1].B, QKS)
				} else {
					tmp := lksk.KeySet[i][j][int64(QKS)-int64(v)+int64(BKS)/2-2].A
					for w := 0; w < n; w++ {
						ct_ks.A[w] = ring.CRed(ct_ks.A[w]+tmp[w], QKS)
					}
					ct_ks.B = ring.CRed(ct_ks.B+lksk.KeySet[i][j][int64(QKS)-int64(v)+int64(BKS)/2-2].B, QKS)
				}
			}
		}
	}

	// ct_ks.B = ring.CRed(ct_ks.B+uint64(math.Round(float64(QKS)/float64((enc.lweParam.P()*4)))), QKS)
	return
}

func Decomp(n, B, QKS uint64, dKS int) []uint64 {
	base_log := uint64(math.Log2(float64(B)))
	digit_mask := uint64(B - 1)
	base_over_2_threshold := int64(B >> 1)
	carry := 0

	var unsigned_digit uint64
	var signed_digit int64

	var digits []uint64 = []uint64{}

	for j := 0; j < dKS; j++ {
		unsigned_digit = (n >> (uint64(j) * base_log)) & digit_mask
		if carry == 1 {
			unsigned_digit += uint64(carry)
			carry = 0
		}

		signed_digit = int64(unsigned_digit)
		if signed_digit >= base_over_2_threshold {
			signed_digit -= int64(B)
			carry = 1
		}

		digits = append(digits, uint64(signed_digit+int64(QKS))%QKS)
	}

	return digits
}

func LWE_MS(ct *Ciphertext, paramsLWE Parameters, Q uint64) (ct_switch *Ciphertext) {
	n := int(paramsLWE.N())
	q := paramsLWE.Q()
	ct_switch = NewCiphertext(n)

	for i := 0; i < n; i++ {
		ct_switch.A[i] = uint64(math.Round(float64(ct.A[i]*q) / float64(Q)))
	}
	ct_switch.B = uint64(math.Round(float64(ct.B*q) / float64(Q)))
	return
}

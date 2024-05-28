package ngsw

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/ntru"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func UNUSED(x ...interface{}) {}

func TestEncryption(t *testing.T) {

	var err error
	var params hefloat.Parameters
	if params, err = hefloat.NewParametersFromLiteral(testParamsLiteral[0]); err != nil {
		panic(err)
	}

	ringQP := params.RingQP()
	ringQ := ringQP.RingQ
	levelQ := ringQ.Level()
	Q := ringQ.ModuliChain()[:levelQ+1]
	N := ringQ.N()

	// plaintext to encrypt
	p_poly := ringQ.NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	UniSampler := ringqp.NewUniformSampler(prng, *ringQP)
	UniSampler.AtLevel(p_poly.Level(), -1).Read(ringqp.Poly{Q: p_poly})

	var pt *rlwe.Plaintext
	if pt, err = rlwe.NewPlaintextAtLevelFromPoly(p_poly.Level(), p_poly); err != nil {
		panic(err)
	}

	for j := range Q {
		for i := 0; i < N; i++ {
			fmt.Printf("%d ", pt.Value.Coeffs[j][i])
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n\n")

	ringQ.NTT(pt.Value, pt.Value)
	ringQ.MForm(pt.Value, pt.Value)

	pt.MetaData.IsNTT = true
	pt.MetaData.IsMontgomery = true

	kgen := ntru.NewKeyGenerator(params)

	// key pair generation using invert
	sk, sk_inv := kgen.GenSecretKeyPairNew()
	encryptor := ntru.NewEncryptor(params, sk_inv)
	decryptor := ntru.NewDecryptor(params, sk)

	var ct *rlwe.Ciphertext
	if ct, err = encryptor.EncryptNew(pt); err != nil {
		panic(err)
	}

	pt_dec := decryptor.DecryptNew(ct)
	ringQ.INTT(pt_dec.Value, pt_dec.Value)
	ringQ.IMForm(pt_dec.Value, pt_dec.Value)

	for j := range Q {
		for i := 0; i < N; i++ {
			fmt.Printf("%d ", pt_dec.Value.Coeffs[j][i])
		}
		fmt.Printf("\n")
	}
}

func TestExternalProduct(t *testing.T) {
	var err error
	var params hefloat.Parameters
	if params, err = hefloat.NewParametersFromLiteral(testParamsLiteral[0]); err != nil {
		panic(err)
	}

	ringQP := params.RingQP()
	ringQ := ringQP.RingQ
	N := ringQ.N()

	kgen := ntru.NewKeyGenerator(params)

	// key pair generation using invert
	sk, sk_inv := kgen.GenSecretKeyPairNew()

	encryptor := ntru.NewEncryptor(params, sk_inv)
	decryptor := ntru.NewDecryptor(params, sk)

	p1_poly := ringQ.NewPoly()
	p2_poly := ringQ.NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	UniSampler := ringqp.NewUniformSampler(prng, *ringQP)
	smallSamplerMontgomeryQ, err := ring.NewSampler(prng, params.RingQ(), params.Xe(), false)
	if err != nil {
		panic(err)
	}

	opt := true
	// case 1 : non-small NTRU w/ small NTRUvec -> fine
	// case 2 : small NTRU w/ non-small NTRUvec -> not fine

	if opt {
		UniSampler.AtLevel(p1_poly.Level(), -1).Read(ringqp.Poly{Q: p1_poly})
		smallSamplerMontgomeryQ.AtLevel(p2_poly.Level()).Read(p2_poly)
	} else {
		smallSamplerMontgomeryQ.AtLevel(p1_poly.Level()).Read(p1_poly)
		UniSampler.AtLevel(p2_poly.Level(), -1).Read(ringqp.Poly{Q: p2_poly})
	}

	var pt1, pt2 *rlwe.Plaintext
	if pt1, err = rlwe.NewPlaintextAtLevelFromPoly(p1_poly.Level(), p1_poly); err != nil {
		panic(err)
	}
	if pt2, err = rlwe.NewPlaintextAtLevelFromPoly(p2_poly.Level(), p2_poly); err != nil {
		panic(err)
	}

	ringQ.NTT(pt1.Value, pt1.Value)
	ringQ.NTT(pt2.Value, pt2.Value)
	ringQ.MForm(pt1.Value, pt1.Value)
	ringQ.MForm(pt2.Value, pt2.Value)

	pt1.MetaData.IsNTT = true
	pt1.MetaData.IsMontgomery = true
	pt2.MetaData.IsNTT = true
	pt2.MetaData.IsMontgomery = true

	T := ringQ.NewPoly() // save the value of pt1 * pt2
	ringQ.MulCoeffsMontgomery(pt1.Value, pt2.Value, T)
	ringQ.INTT(T, T)
	ringQ.IMForm(T, T)

	pRLWE := *params.GetRLWEParameters()
	NTRUVencryptor := NewEncryptor(pRLWE, sk_inv)
	ev := NewEvaluator(params, nil)
	var ct1 *rlwe.Ciphertext
	ct2 := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 4)

	if ct1, err = encryptor.EncryptNew(pt1); err != nil {
		panic(err)
	}
	if err = NTRUVencryptor.Encrypt(pt2, ct2); err != nil {
		panic(err)
	}

	ev.ExternalProduct(ct1, ct2, ct1)

	pt_dec := decryptor.DecryptNew(ct1)
	ringQ.INTT(pt_dec.Value, pt_dec.Value)
	ringQ.IMForm(pt_dec.Value, pt_dec.Value)

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", pt_dec.Value.Coeffs[0][i])
	}
	fmt.Printf("\n")
	fmt.Printf("\n")

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", T.Coeffs[0][i])
	}
}

// testing polynomial (not encrypted) * non-small NTRU-vec (includes f_inv term)
// should be shifted left by 1 in the test case below
func TestPolyMult(t *testing.T) {
	var err error
	var params hefloat.Parameters
	if params, err = hefloat.NewParametersFromLiteral(testParamsLiteral[0]); err != nil {
		panic(err)
	}

	ringQP := params.RingQP()
	ringQ := ringQP.RingQ
	N := ringQ.N()

	kgen := ntru.NewKeyGenerator(params)

	// key pair generation using invert
	sk, sk_inv := kgen.GenSecretKeyPairNew()

	encryptor := ntru.NewEncryptor(params, sk_inv)
	decryptor := ntru.NewDecryptor(params, sk)

	p1_poly := ringQ.NewPoly()
	p2_poly := ringQ.NewPoly()
	T := ringQ.NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	UniSampler := ringqp.NewUniformSampler(prng, *ringQP)
	UniSampler.AtLevel(p1_poly.Level(), -1).Read(ringqp.Poly{Q: p1_poly})
	T.Copy(p1_poly)

	// set the p2 value as X/f
	tmp := ringQ.NewMonomialXi(1)
	ringQ.MForm(tmp, tmp)
	ringQ.NTT(tmp, tmp)

	p2_poly.Copy(sk_inv.Value.Q)
	ringQ.MulCoeffsMontgomery(p2_poly, tmp, p2_poly)

	var pt1, pt2 *rlwe.Plaintext
	if pt1, err = rlwe.NewPlaintextAtLevelFromPoly(p1_poly.Level(), p1_poly); err != nil {
		panic(err)
	}
	if pt2, err = rlwe.NewPlaintextAtLevelFromPoly(p2_poly.Level(), p2_poly); err != nil {
		panic(err)
	}

	ringQ.NTT(pt1.Value, pt1.Value)
	ringQ.MForm(pt1.Value, pt1.Value)

	pt1.MetaData.IsNTT = true
	pt1.MetaData.IsMontgomery = true
	pt2.MetaData.IsNTT = true
	pt2.MetaData.IsMontgomery = true

	pRLWE := *params.GetRLWEParameters()
	NTRUVencryptor := NewEncryptor(pRLWE, sk_inv)
	ev := NewEvaluator(params, nil)
	ct1 := rlwe.NewCiphertext(params, 1, params.MaxLevel())
	ct2 := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), 16)

	if err = encryptor.EncryptPoly(pt1, ct1); err != nil {
		panic(err)
	}
	if err = NTRUVencryptor.Encrypt(pt2, ct2); err != nil {
		panic(err)
	}

	ev.ExternalProduct(ct1, ct2, ct1)

	pt_dec := decryptor.DecryptNew(ct1)
	ringQ.INTT(pt_dec.Value, pt_dec.Value)
	ringQ.IMForm(pt_dec.Value, pt_dec.Value)

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", pt_dec.Value.Coeffs[0][i])
	}
	fmt.Printf("\n")
	fmt.Printf("\n")

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", T.Coeffs[0][i])
	}
}

func TestKeySwitching(t *testing.T) {
	var err error
	var params hefloat.Parameters
	if params, err = hefloat.NewParametersFromLiteral(testParamsLiteral[0]); err != nil {
		panic(err)
	}

	ringQP := params.RingQP()
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP
	N := ringQ.N()

	// key pair generation using invert
	kgen := ntru.NewKeyGenerator(params)

	sk, sk_inv := kgen.GenSecretKeyPairNew()
	sk2, sk_inv2 := kgen.GenSecretKeyPairNew()

	encryptor := ntru.NewEncryptor(params, sk_inv)
	decryptor := ntru.NewDecryptor(params, sk2)

	p1_poly := ringQ.NewPoly()
	p2_poly := ringQ.NewPoly()
	T := ringQ.NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	UniSampler := ringqp.NewUniformSampler(prng, *ringQP)
	// ternarySamplerMontgomeryQ, err := ring.NewSampler(prng, ringQ, ring.Ternary{P: 1.0 / 3.0}, false)

	UniSampler.AtLevel(p1_poly.Level(), -1).Read(ringqp.Poly{Q: p1_poly})
	T.Copy(p1_poly)

	// set the p2 value as sk1 / sk2
	tmp := ringP.NewPoly()
	ringQP.MulCoeffsMontgomery(sk.Value, sk_inv2.Value, ringqp.Poly{Q: p2_poly, P: tmp})

	var pt1, pt2 *rlwe.Plaintext
	if pt1, err = rlwe.NewPlaintextAtLevelFromPoly(p1_poly.Level(), p1_poly); err != nil {
		panic(err)
	}
	if pt2, err = rlwe.NewPlaintextAtLevelFromPoly(p2_poly.Level(), p2_poly); err != nil {
		panic(err)
	}

	ringQ.NTT(pt1.Value, pt1.Value)
	ringQ.MForm(pt1.Value, pt1.Value)

	pt1.MetaData.IsNTT = true
	pt1.MetaData.IsMontgomery = true
	pt2.MetaData.IsNTT = true
	pt2.MetaData.IsMontgomery = true

	pRLWE := *params.GetRLWEParameters()
	NTRUVencryptor := NewEncryptor(pRLWE, sk_inv2)
	ev := NewEvaluator(params, nil)
	basetwo := 16

	var ct1 *rlwe.Ciphertext
	ct2 := NewCiphertext(params, params.MaxLevelQ(), params.MaxLevelP(), basetwo)

	if ct1, err = encryptor.EncryptNew(pt1); err != nil {
		panic(err)
	}
	if err = NTRUVencryptor.Encrypt(pt2, ct2); err != nil {
		panic(err)
	}

	ev.ExternalProduct(ct1, ct2, ct1)

	pt_dec := decryptor.DecryptNew(ct1)
	ringQ.INTT(pt_dec.Value, pt_dec.Value)
	ringQ.IMForm(pt_dec.Value, pt_dec.Value)

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", pt_dec.Value.Coeffs[0][i])
	}
	fmt.Printf("\n")
	fmt.Printf("\n")

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", T.Coeffs[0][i])
	}
}

// 2N과 서로소인 수 ex) 3
// X^i -> X^3i를 기존 key로 decrypt하도록 만들어야 한다.
func TestAutomorphism(t *testing.T) {
	var err error
	var params hefloat.Parameters
	if params, err = hefloat.NewParametersFromLiteral(testParamsLiteral[0]); err != nil {
		panic(err)
	}

	ringQP := params.RingQP()
	ringQ := ringQP.RingQ
	// ringP := ringQP.RingP

	kgen := ntru.NewKeyGenerator(params)

	pRLWE := *params.GetRLWEParameters()
	ntruv_kgen := NewKeyGenerator(pRLWE)

	N := params.N()
	galEl := make([]uint64, 2*N-1)

	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	sk, sk_inv := kgen.GenSecretKeyPairNew()
	encryptor := ntru.NewEncryptor(params, sk_inv)
	decryptor := ntru.NewDecryptor(params, sk)

	basetwo := 16
	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(params.MaxLevelQ()), LevelP: utils.Pointy(params.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(basetwo)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ntruv_kgen.GenAutoKeysNew(galEl, sk, sk_inv, evkParams)...)

	ev := NewEvaluator(params, evk)

	/// polynomial setup

	p1_poly := ringQ.NewPoly()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	UniSampler := ringqp.NewUniformSampler(prng, *ringQP)
	UniSampler.AtLevel(p1_poly.Level(), -1).Read(ringqp.Poly{Q: p1_poly})

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", p1_poly.Coeffs[0][i])
	}
	fmt.Printf("\n\n")

	var pt *rlwe.Plaintext
	if pt, err = rlwe.NewPlaintextAtLevelFromPoly(p1_poly.Level(), p1_poly); err != nil {
		panic(err)
	}

	ringQ.NTT(pt.Value, pt.Value)
	ringQ.MForm(pt.Value, pt.Value)

	pt.MetaData.IsNTT = true
	pt.MetaData.IsMontgomery = true

	var ct *rlwe.Ciphertext
	if ct, err = encryptor.EncryptNew(pt); err != nil {
		panic(err)
	}

	ct2 := ct.CopyNew()
	ev.Automorphism(ct, 5, ct2)

	pt_dec := decryptor.DecryptNew(ct2)
	ringQ.INTT(pt_dec.Value, pt_dec.Value)
	ringQ.IMForm(pt_dec.Value, pt_dec.Value)

	for i := 0; i < N; i++ {
		fmt.Printf("%d ", pt_dec.Value.Coeffs[0][i])
	}
}

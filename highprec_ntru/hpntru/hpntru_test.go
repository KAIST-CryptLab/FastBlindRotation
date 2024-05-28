package hpntru

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"

	"github.com/montanaflynn/stats"
	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/ngsw"
	"github.com/tuneinsight/lattigo/v5/core/ntru"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func UNUSED(x ...interface{}) {}

// log p, log q, n, QKS, logB, logBKS, stdlwe_lksk, sk_dist, tn, BR param idx
func TestNTRUP3(t *testing.T) {
	testHPNTRU(t, 3, 10, 512, 65536, 3, 7, 2., 1./2., 10000, 0)
}

func TestNTRUErrP3(t *testing.T) {
	testNTRUBRErrorAnalysis(t, 3, 10, 570, 65536, 3, 7, 2., 0)
}

// ////////////////////////////////////////////////////////////////////
func TestNTRUP4_1(t *testing.T) {
	testHPNTRU(t, 4, 10, 512, 65536, 3, 7, 2., 1./2., 10000, 1)
}

func TestNTRUP4_2(t *testing.T) {
	testHPNTRU(t, 5, 11, 1024, 524288, 4, 7, 3., 1./4., 5000, 2)
}

func TestNTRUErrP4_1(t *testing.T) {
	testNTRUBRErrorAnalysis(t, 4, 10, 570, 65536, 3, 7, 2., 1)
}

func TestNTRUErrP4_2(t *testing.T) {
	testNTRUBRErrorAnalysis(t, 4, 11, 570, 65536, 4, 7, 1.5, 2)
}

// ////////////////////////////////////////////////////////////////////
// log p, log q, d, n, QKS, logB, logBKS, stdlwe_lksk, bound_lksk, test case, BR param idx, verbose
/*
func TestNTRUTreeP3(t *testing.T) {
	testTreeBSNTRU(t, 3, 10, 2, 570, 65536, 2, 5, 1.5, 1000, 0, false)
}

// Rather Unstable with q = 1024, not useful (QKS could be 557057)
func TestNTRUTreeP4_1(t *testing.T) {
	testTreeBSNTRU(t, 4, 10, 3, 570, 65536, 2, 5, 1.5, 100, 1, false)
}
*/

// Stable, but slow
func TestNTRUTreeP4_2(t *testing.T) {
	testTreeBSNTRU(t, 4, 11, 4, 1024, 524288, 2, 7, 3., 100, 2, false)
}

func BenchmarkNTRUTreeP4(b *testing.B) {
	benchmarkTreeBSNTRU(b, 4, 11, 5, 1024, 524288, 2, 7, 3., 100, 2, false)
}

// ////////////////////////////////////////////////////////////////////
func BenchmarkCompRLWENTRUP3(b *testing.B) {
	benchmarkNTRURLWE(b, 3, 10, 3, 512, 1./2., 1000, 0)
}

func BenchmarkCompRLWENTRUP4(b *testing.B) {
	benchmarkNTRURLWE(b, 4, 11, 4, 1024, 1./4., 100, 2)
}

func testHPNTRU(t *testing.T, logp, logq int, n, QKS uint64, logB, logBKS int, stdev_lksk float64, sk_dist float64, tn, idx int) {
	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << logp // plaintext Modulus
	var q uint64 = 1 << logq // ciphertext Modulus

	var BKS = 1 << logBKS // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(t, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: sk_dist}, false)
	require.NoError(t, err)
	skLWE := ternarySamplerLWE.ReadNew()

	skLWECoeffs := make([]*big.Int, ringQLWE.N())
	skLWEfloat := make([]float64, ringQLWE.N())
	for i := range skLWECoeffs {
		skLWECoeffs[i] = new(big.Int)
	}
	ringQLWE.PolyToBigintCentered(skLWE, 1, skLWECoeffs)
	var sk_ss float64 = 0
	for i := range skLWECoeffs {
		skLWEfloat[i] = float64(skLWECoeffs[i].Int64())
		if i < int(n) {
			sk_ss += skLWEfloat[i] * skLWEfloat[i]
		}
	}

	skvar, err := stats.StandardDeviation(skLWEfloat)
	if err != nil {
		panic(err)
	}
	fmt.Printf("lwe sk variance : %f\n", skvar)
	fmt.Printf("lwe sk ss : %f\n", sk_ss)

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)
	decLWE := lwe.NewDecoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testNTRUParamsLiteral[idx]) // N, Q
	require.NoError(t, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsBR).GenSecretKeyPairNew()

	N := paramsBR.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsBR).GenAutoKeysNew(galEl, f, f_inv, evkParams)...)
	fmt.Printf("Key Finished\n")
	BRK := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, f_inv, evkParams)

	// Encryptor, Decryptor for NTRU BR
	// encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	// decryptorBR := ntru.NewDecryptor(paramsBR, f)

	// Evaluator for BR
	evalBR := NewEvaluator(paramsBR, paramsLWE, evk)

	// TestPolynomial Used for rotation
	F := InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ())

	//////////////////////////////////////////
	// KeySwitching Context, happens at QKS //
	//////////////////////////////////////////
	fmt.Printf("Key Switching Context Generation\n")
	paramsKSLWE := lwe.NewParameters(n, QKS, p)

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)

	encryptorKSLWE := lwe.NewEncryptor(paramsKSLWE, ringQLWE, skLWE)
	// decryptorKSLWE := lwe.NewDecryptor(paramsKSLWE, ringQLWE, skLWE)

	errGenLKSK := lwe.NewErrorGenerator(stdev_lksk)

	lksk := lwe.LWE_KSKG(BKS, f, encryptorKSLWE, paramsBR.RingQ(), errGenLKSK) // Generate KeySwitching Key

	///////////////////////////////////
	// Setup finished, testing start //
	///////////////////////////////////
	fmt.Printf("Setup finished, start bootstrapping\n")

	var cnt = 0
	var errs = make([]float64, tn)
	var Success = map[bool]string{true: "Success", false: "Failed"}

	fmt.Printf("Current Precision : %d\n", logp)
	for ti := 0; ti < tn; ti++ {

		var m uint64 = uint64(rand.Intn(int(p))) // message to encrypt/decrypt
		ct := encryptorLWE.EncryptNew(encLWE.Encode(m, errGenLWE.GenErr()))

		// Do BlindRotation
		ct_br, err := evalBR.BlindRotate(ct, &F, BRK) // ct_br is NTRU cipher
		require.NoError(t, err)

		// ModSwitch, KeySwitch
		ct_ks := lwe.LWE_KS(ct_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher

		// Modswitch again to original q
		ct_boot := lwe.LWE_MS(ct_ks, paramsLWE, QKS)

		m_enc_boot := decryptorLWE.DecryptNew(ct_boot)
		m_dec := decLWE.Decode(m_enc_boot)

		// fmt.Printf("Bootstrap : %2d -> %2d (%s)\n", m, m_dec, Success[m_dec == m])
		if m_dec == m {
			cnt++
		}

		// Error Calculation
		if m_dec == m {
			err_calc := int64(m_enc_boot - (q / (2 * p) * m))
			if err_calc > int64(q/2) {
				err_calc = err_calc - int64(q)
			}
			errs[ti] = float64(err_calc)
		}
		/*
			else {
				fmt.Printf("Error : \n%b\n%b\n%f\n", m_enc_boot, (q / (2 * p) * m), errs[ti])
			}
		*/
	}

	// fmt.Printf("%v\n", errs)

	em, es := GetMeanStdev(errs)
	fmt.Printf("\nError Stats : %f, %f\n", em, es)
	fmt.Printf("Dec Success : %d/%d (%f%%)\n", cnt, tn, float64(cnt)/float64(tn)*100)
	UNUSED(Success)
}

func testNTRUBRErrorAnalysis(t *testing.T, logp, logq int, n, QKS uint64, logB, logBKS int, stdev_lksk float64, idx int) {
	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << logp // plaintext Modulus
	var q uint64 = 1 << logq // ciphertext Modulus

	var BKS = 1 << logBKS // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2
	var tn = 1

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(t, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: 1.0 / 2.0}, false)
	require.NoError(t, err)
	skLWE := ternarySamplerLWE.ReadNew()

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)
	decLWE := lwe.NewDecoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testNTRUParamsLiteral[idx]) // N, Q
	require.NoError(t, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsBR).GenSecretKeyPairNew()

	N := paramsBR.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsBR).GenAutoKeysNew(galEl, f, f_inv, evkParams)...)
	fmt.Printf("Key Finished\n")
	BRK := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, f_inv, evkParams)

	// Encryptor, Decryptor for NTRU BR
	encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	decryptorBR := ntru.NewDecryptor(paramsBR, f)

	// Evaluator for BR
	evalBR := NewEvaluator(paramsBR, paramsLWE, evk)

	// TestPolynomial Used for rotation
	F := InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ())

	//////////////////////////////////////////
	// KeySwitching Context, happens at QKS //
	//////////////////////////////////////////
	fmt.Printf("Key Switching Context Generation\n")
	paramsKSLWE := lwe.NewParameters(n, QKS, p)

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)

	encryptorKSLWE := lwe.NewEncryptor(paramsKSLWE, ringQLWE, skLWE)
	decryptorKSLWE := lwe.NewDecryptor(paramsKSLWE, ringQLWE, skLWE)

	errGenLKSK := lwe.NewErrorGenerator(stdev_lksk)

	lksk := lwe.LWE_KSKG(BKS, f, encryptorKSLWE, paramsBR.RingQ(), errGenLKSK) // Generate KeySwitching Key

	///////////////////////////////////
	// Setup finished, testing start //
	///////////////////////////////////
	fmt.Printf("Setup finished, start bootstrapping\n")

	var errs = make([]float64, tn)
	var Success = map[bool]string{true: "Success", false: "Failed"}

	fmt.Printf("Current Precision : %d\n\n", logp)
	for ti := 0; ti < tn; ti++ {

		var m uint64 = uint64(rand.Intn(int(p))) // message to encrypt/decrypt
		ct := encryptorLWE.EncryptNew(encLWE.Encode(m, errGenLWE.GenErr()))

		// Do BlindRotation
		ct_br, err := evalBR.BlindRotate(ct, &F, BRK) // ct_br is NTRU cipher
		require.NoError(t, err)

		// Error after BR
		pt_br := decryptorBR.DecryptNew(ct_br)
		ringQBR.INTT(pt_br.Value, pt_br.Value)

		exp_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m)
		err_br := int64(pt_br.Value.Coeffs[0][0]) - exp_br
		if err_br > int64(paramsBR.Q()[0]/2) {
			err_br = -int64(paramsBR.Q()[0]) + err_br
		}
		limit_br := float64(paramsBR.Q()[0]) / float64(4*p)
		fmt.Printf("Expected    : %020b\nActual      : %020b\n", exp_br, int64(pt_br.Value.Coeffs[0][0]))
		fmt.Printf("BR Error    : %d/%d (%f%%)\n\n", err_br, int64(limit_br), math.Abs(float64(err_br))/limit_br*100.0)

		ct_ks := lwe.LWE_KS(ct_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk)
		pt_ks := decryptorKSLWE.DecryptNew(ct_ks)

		exp_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m)
		err_ks := int64(pt_ks) - exp_ks
		if err_ks > int64(QKS/2) {
			err_ks = -int64(QKS) + err_ks
		}
		limit_ks := float64(QKS) / float64(4*p)
		fmt.Printf("Expected    : %020b\nActual      : %020b\n", exp_ks, int64(pt_ks))
		fmt.Printf("MS-KS Error : %d/%d (%f%%)\n\n", err_ks, int64(limit_ks), math.Abs(float64(err_ks))/limit_ks*100.0)

		// Modswitch again to original q
		ct_boot := lwe.LWE_MS(ct_ks, paramsLWE, QKS)
		pt_boot := decryptorLWE.DecryptNew(ct_boot)

		exp_boot := int64(uint64(math.Round(float64(q)/float64(2*p))) * m)
		err_boot := int64(pt_boot) - exp_boot
		if err_boot > int64(q/2) {
			err_boot = -int64(q) + err_boot
		}
		limit_boot := float64(q) / float64(4*p)
		fmt.Printf("Expected    : %020b\nActual      : %020b\n", exp_boot, int64(pt_boot))
		fmt.Printf("MS Error    : %d/%d (%f%%)\n\n", err_boot, int64(limit_boot), math.Abs(float64(err_boot))/limit_boot*100.0)

		m_boot := decLWE.Decode(pt_boot)
		fmt.Printf("Bootstrap   : %d -> %d (%s)\n", m, m_boot, Success[m_boot == m])
	}

	UNUSED(decLWE, encryptorBR, decryptorLWE, decryptorKSLWE, lksk, errs, Success)
}

/*
func TestTreeBSNTRU(t *testing.T) {
	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << 3  // plaintext Modulus
	var q uint64 = 1 << 10 // ciphertext Modulus
	var n uint64 = 570     // TFHE dimension

	var QKS uint64 = 65536 //  // Key Switching Modulus

	var logB = 2     // Decomposition Base used for Blind Rotation B = 2^logB
	var BKS = 1 << 7 // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2
	var bound_lwe = int64(10)

	var stdev_lksk = 3.2
	var bound_lksk = int64(5)

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{0x3001})
	require.NoError(t, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: 1.0 / 2.0}, false)
	require.NoError(t, err)
	skLWE := ternarySamplerLWE.ReadNew()

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe, bound_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)
	decLWE := lwe.NewDecoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testNTRUParamsLiteral[0]) // N, Q
	require.NoError(t, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsBR).GenSecretKeyPairNew()

	N := paramsBR.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsBR).GenAutoKeysNew(galEl, f, f_inv, evkParams)...)
	BRK := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, f_inv, evkParams)

	// Encryptor, Decryptor for NTRU BR
	// encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	decryptorBR := ntru.NewDecryptor(paramsBR, f)

	//////////////////////////////////////////
	// KeySwitching Context, happens at QKS //
	//////////////////////////////////////////
	fmt.Printf("Key Switching Context Generation\n")
	paramsKSLWE := lwe.NewParameters(n, QKS, p)

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)

	encryptorKSLWE := lwe.NewEncryptor(paramsKSLWE, ringQLWE, skLWE)
	decryptorKSLWE := lwe.NewDecryptor(paramsKSLWE, ringQLWE, skLWE)

	errGenLKSK := lwe.NewErrorGenerator(stdev_lksk, bound_lksk)

	lksk := lwe.LWE_KSKG(BKS, f, encryptorKSLWE, paramsBR.RingQ(), errGenLKSK) // Generate KeySwitching Key

	// Evaluator for BR
	evalBR_Tree := NewTreeEvaluator(paramsBR, paramsLWE, evk)

	F_high, F_low := InitReLUPolynomial(p, q, ringQBR)

	var m_high = uint64(1) // uint64(rand.Intn(int(p)))
	var m_high_relu = m_high
	if m_high >= p/2 {
		m_high_relu = 0
	}
	ct_high := encryptorLWE.EncryptNew(encLWE.Encode(m_high, 0))

	// The high bits are processed using the same functions as single message
	// Along with it, the multi-value BR is calculated to get the low bit test polynomial
	ct_high_br, F_lows, err := evalBR_Tree.MultiValueBR(ct_high, &F_high, F_low, BRK) // ct_br is NTRU cipher
	require.NoError(t, err)

	// Error after BR
	pt_br := decryptorBR.DecryptNew(ct_high_br)
	ringQBR.INTT(pt_br.Value, pt_br.Value)

	exp_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_high_relu)
	fmt.Printf("Expected : %d\n", exp_br)
	// fmt.Printf("%v\n", pt_br.Value.Coeffs[0])


	F_low_enc := evalBR_Tree.CombineTestPoly(F_lows, int(p), N)
	pt_F := decryptorBR.DecryptNew(F_low_enc)
	ringQBR.INTT(pt_F.Value, pt_F.Value)

	// ringQBR.IMForm(pt_F.Value, pt_F.Value)

	fmt.Printf("%v\n", pt_F.Value.Coeffs[0])

	// Using the key for encrypted TP, we perform the blind rotation for low bits

	// ringQBR.MForm(F_low_enc.Value[1], F_low_enc.Value[1])

	var m_low = uint64(rand.Intn(int(p)))
	var m_low_relu = m_low
	if m_high >= p/2 {
		m_low_relu = 0
	}
	ct_low := encryptorLWE.EncryptNew(encLWE.Encode(m_low, errGenLWE.GenErr()))
	ct_low_br, err := evalBR_Tree.EncBlindRotate(ct_low, F_low_enc, BRK)
	require.NoError(t, err)

	pt_low_br := decryptorBR.DecryptNew(ct_low_br)
	ringQBR.INTT(pt_low_br.Value, pt_low_br.Value)

	exp_low_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_low_relu)
	fmt.Printf("Expected : %d\n", exp_low_br)
	fmt.Printf("%v\n", pt_low_br.Value.Coeffs[0])

	ct_low_ks := lwe.LWE_KS(ct_low_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher
	// ct_lows_boot[l] = lwe.LWE_MS(ct_low_ks, paramsLWE, QKS)

	pt_low_ks := decryptorKSLWE.DecryptNew(ct_low_ks)

	exp_low_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m_low_relu)
	err_low_ks := int64(pt_low_ks) - exp_low_ks
	if err_low_ks > int64(QKS/2) {
		err_low_ks = int64(QKS) - err_low_ks
	}
	limit_low_ks := float64(QKS) / float64(4*p)
	fmt.Printf("MS-KS Error : %8d among %8d (%f%%)\n", err_low_ks, int64(limit_low_ks), math.Abs(float64(err_low_ks))/limit_low_ks*100.0)

	// ct_lows_boot[l] = lwe.LWE_MS(ct_low_ks, paramsLWE, QKS)
	UNUSED(F_low, F_lows, decLWE, decryptorLWE, decryptorKSLWE, lksk, errGenLWE)
}
*/

// Tree Method, tested in ReLU
func testTreeBSNTRU(t *testing.T, logp, logq, d int, n, QKS uint64, logB, logBKS int, stdev_lksk float64, tn, idx int, debug bool) {
	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << logp // plaintext Modulus
	var q uint64 = 1 << logq // ciphertext Modulus

	var BKS = 1 << logBKS // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(t, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: 1. / 4.}, false)
	require.NoError(t, err)
	skLWE := ternarySamplerLWE.ReadNew()

	var sk_ss float64 = 0
	skLWECoeffs := make([]*big.Int, ringQLWE.N())
	for i := range skLWECoeffs {
		skLWECoeffs[i] = new(big.Int)
	}
	ringQLWE.PolyToBigintCentered(skLWE, 1, skLWECoeffs)
	for i := 0; i < int(n); i++ {
		sk_ss += float64(skLWECoeffs[i].Int64()) * float64(skLWECoeffs[i].Int64())
	}

	fmt.Printf("lwe sk ss : %f\n", sk_ss)

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)
	decLWE := lwe.NewDecoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testMVBRParamsLiteral[idx]) // N, Q
	require.NoError(t, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsBR).GenSecretKeyPairNew()

	N := paramsBR.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsBR).GenAutoKeysNew(galEl, f, f_inv, evkParams)...)
	BRK := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, f_inv, evkParams)

	// Encryptor, Decryptor for NTRU BR
	// encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	decryptorBR := ntru.NewDecryptor(paramsBR, f)

	// Evaluator for BR
	evalBR_Tree := NewTreeEvaluator(paramsBR, paramsLWE, evk)

	// TestPolynomial Used for rotation
	F_high, F_low := InitReLUPolynomial(p, q, ringQBR)

	//////////////////////////////////////////
	// KeySwitching Context, happens at QKS //
	//////////////////////////////////////////
	fmt.Printf("Key Switching Context Generation\n")
	paramsKSLWE := lwe.NewParameters(n, QKS, p)

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)

	encryptorKSLWE := lwe.NewEncryptor(paramsKSLWE, ringQLWE, skLWE)
	decryptorKSLWE := lwe.NewDecryptor(paramsKSLWE, ringQLWE, skLWE)

	errGenLKSK := lwe.NewErrorGenerator(stdev_lksk)

	lksk := lwe.LWE_KSKG(BKS, f, encryptorKSLWE, paramsBR.RingQ(), errGenLKSK) // Generate KeySwitching Key

	///////////////////////////////////
	// Setup finished, testing start //
	///////////////////////////////////
	fmt.Printf("Setup finished, start bootstrapping\n")

	var P = int(math.Pow(float64(p), float64(d)))
	fmt.Printf("Current Precision : %d\n", int(math.Log2(float64(p)))*d)

	var errs1 []float64 // make([]float64, tn)
	var errs2 []float64
	var Success = map[bool]string{true: "Success", false: "Failed"}

	var cnt = 0
	for ti := 0; ti < tn; ti++ {
		var m = uint64(rand.Intn(P)) // message to encrypt/decrypt
		var m_decomp = Decomp(m, p, d)

		m_high := m_decomp[0]
		m_lows := m_decomp[1:]

		var m_high_relu = m_high
		if m_high >= uint64(p)/2 {
			m_high_relu = 0
		}

		// test using ReLU (0 <= x < P/2 -> x, P/2 <= x -> 0)
		ct_high := encryptorLWE.EncryptNew(encLWE.Encode(m_high, errGenLWE.GenErr()))

		// The high bits are processed using the same functions as single message
		// Along with it, the multi-value BR is calculated to get the low bit test polynomial
		ct_high_br, F_lows, err := evalBR_Tree.MultiValueBR(ct_high, &F_high, F_low, BRK) // ct_br is NTRU cipher
		require.NoError(t, err)

		// Error after BR
		pt_high_br := decryptorBR.DecryptNew(ct_high_br)
		ringQBR.INTT(pt_high_br.Value, pt_high_br.Value)

		if debug {
			exp_high_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_high_relu)
			err_high_br := int64(pt_high_br.Value.Coeffs[0][0]) - exp_high_br
			if err_high_br > int64(paramsBR.Q()[0]/2) {
				err_high_br = -int64(paramsBR.Q()[0]) + err_high_br
			}
			limit_high_br := float64(paramsBR.Q()[0]) / float64(4*p)
			fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_br, int64(pt_high_br.Value.Coeffs[0][0]))
			fmt.Printf("(H) BR Error    : %d/%d (%f%%)\n\n", err_high_br, int64(limit_high_br), math.Abs(float64(err_high_br))/limit_high_br*100.0)
		}

		ct_high_ks := lwe.LWE_KS(ct_high_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher

		if debug {
			pt_high_ks := decryptorKSLWE.DecryptNew(ct_high_ks)
			exp_high_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m_high_relu)
			err_high_ks := int64(pt_high_ks) - exp_high_ks
			if err_high_ks > int64(QKS/2) {
				err_high_ks = -int64(QKS) + err_high_ks
			}
			limit_high_ks := float64(QKS) / float64(4*p)
			fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_ks, int64(pt_high_ks))
			fmt.Printf("(H) MS-KS Error : %d/%d (%f%%)\n\n", err_high_ks, int64(limit_high_ks), math.Abs(float64(err_high_ks))/limit_high_ks*100.0)
		}

		ct_high_boot := lwe.LWE_MS(ct_high_ks, paramsLWE, QKS)
		pt_high_boot := decryptorLWE.DecryptNew(ct_high_boot)
		m_high_dec := decLWE.Decode(pt_high_boot)

		if debug {
			exp_high_boot := int64(uint64(math.Round(float64(q)/float64(2*p))) * m_high_relu)
			err_high_boot := int64(pt_high_boot) - exp_high_boot
			if err_high_boot > int64(q/2) {
				err_high_boot = -int64(q) + err_high_boot
			}
			limit_high_boot := float64(q) / float64(4*p)

			fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_boot, int64(pt_high_boot))
			fmt.Printf("(H) MS Error    : %d/%d (%f%%)\n\n", err_high_boot, int64(limit_high_boot), math.Abs(float64(err_high_boot))/limit_high_boot*100.0)

			fmt.Printf("(H) %d -> %d\n\n", m_high, m_high_dec)
			fmt.Printf("=================================\n\n")
		}

		if m_high_dec == m_high_relu {
			err_calc := int64(pt_high_boot - (q / (2 * p) * m_high_relu))
			if err_calc > int64(q/2) {
				err_calc = err_calc - int64(q)
			}
			errs1 = append(errs1, float64(err_calc))
		}

		F_low_enc := evalBR_Tree.CombineTestPoly(F_lows, int(p), N)
		var m_dec = m_high_dec

		for _, m_low := range m_lows {
			m_low_relu := m_low
			if m_high >= uint64(p)/2 {
				m_low_relu = 0
			}
			ct_low := encryptorLWE.EncryptNew(encLWE.Encode(m_low, errGenLWE.GenErr()))

			// Using the key for encrypted TP, we perform the blind rotation for low bits
			ct_low_br, err := evalBR_Tree.EncBlindRotate(ct_low, F_low_enc, BRK)
			require.NoError(t, err)

			if debug {
				pt_low_br := decryptorBR.DecryptNew(ct_low_br)
				ringQBR.INTT(pt_low_br.Value, pt_low_br.Value)

				exp_low_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_low_relu)
				err_low_br := int64(pt_low_br.Value.Coeffs[0][0]) - exp_low_br
				if err_low_br > int64(paramsBR.Q()[0]/2) {
					err_low_br = -int64(paramsBR.Q()[0]) + err_low_br
				}
				limit_low_br := float64(paramsBR.Q()[0]) / float64(4*p)

				fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_br, int64(pt_low_br.Value.Coeffs[0][0]))
				fmt.Printf("(L) BR Error    : %d/%d (%f%%)\n\n", err_low_br, int64(limit_low_br), math.Abs(float64(err_low_br))/limit_low_br*100.0)
			}

			ct_low_ks := lwe.LWE_KS(ct_low_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher

			if debug {
				pt_low_ks := decryptorKSLWE.DecryptNew(ct_low_ks)
				exp_low_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m_low_relu)
				err_low_ks := int64(pt_low_ks) - exp_low_ks
				if err_low_ks > int64(QKS/2) {
					err_low_ks = -int64(QKS) + err_low_ks
				}
				limit_low_ks := float64(QKS) / float64(4*p)

				fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_ks, int64(pt_low_ks))
				fmt.Printf("(L) MS-KS Error : %d/%d (%f%%)\n\n", err_low_ks, int64(limit_low_ks), math.Abs(float64(err_low_ks))/limit_low_ks*100.0)
			}

			ct_low_boot := lwe.LWE_MS(ct_low_ks, paramsLWE, QKS)
			pt_low_boot := decryptorLWE.DecryptNew(ct_low_boot)
			m_low_dec := decLWE.Decode(pt_low_boot)

			if debug {
				exp_low_boot := int64(uint64(math.Round(float64(q)/float64(2*p))) * m_low_relu)
				err_low_boot := int64(pt_low_boot) - exp_low_boot
				if err_low_boot > int64(q/2) {
					err_low_boot = -int64(q) + err_low_boot
				}
				limit_low_boot := float64(q) / float64(4*p)
				fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_boot, int64(pt_low_boot))
				fmt.Printf("(L) MS Error    : %d/%d (%f%%)\n\n", err_low_boot, int64(limit_low_boot), math.Abs(float64(err_low_boot))/limit_low_boot*100.0)

				fmt.Printf("(L) %d -> %d\n\n", m_low, m_low_dec)
				fmt.Printf("=================================\n\n")
			}

			if m_low_dec == m_low_relu {
				err_calc := int64(pt_low_boot - (q / (2 * p) * m_low_relu))
				if err_calc > int64(q/2) {
					err_calc = err_calc - int64(q)
				}
				errs2 = append(errs2, float64(err_calc))
			}

			m_dec = m_dec*p + m_low_dec
		}

		if m >= uint64(P)/2 {
			// fmt.Printf("Bootstrap : %5d (>= %5d) -> %5d (%s)\n", m, uint64(P)/2, m_dec, Success[m_dec == 0])
			if m_dec == 0 {
				cnt++
			}
		} else {
			// fmt.Printf("Bootstrap : %5d (<  %5d) -> %5d (%s)\n", m, uint64(P)/2, m_dec, Success[m_dec == m])
			if m_dec == m {
				cnt++
			}
		}
	}

	em1, es1 := GetMeanStdev(errs1)
	fmt.Printf("\nError Stats : %f, %f\n", em1, es1)

	em2, es2 := GetMeanStdev(errs2)
	fmt.Printf("\nError Stats : %f, %f\n", em2, es2)

	// fmt.Printf("%v\n", norm_squares)
	nsm, _ := GetMeanStdev(norm_squares)
	fmt.Printf("\nNorm Square Mean : %f\n", nsm)

	fmt.Printf("Dec Success : %d/%d (%f%%)", cnt, tn, float64(cnt)/float64(tn)*100)
	UNUSED(errs1, errs2, Success)
}

// Tree Method, tested in ReLU
func benchmarkTreeBSNTRU(b *testing.B, logp, logq, d int, n, QKS uint64, logB, logBKS int, stdev_lksk float64, tn, idx int, debug bool) {
	prng, err := sampling.NewPRNG()
	require.NoError(b, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << logp // plaintext Modulus
	var q uint64 = 1 << logq // ciphertext Modulus

	var BKS = 1 << logBKS // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(b, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: 1. / 4.}, false)
	require.NoError(b, err)
	skLWE := ternarySamplerLWE.ReadNew()

	var sk_ss float64 = 0
	skLWECoeffs := make([]*big.Int, ringQLWE.N())
	for i := range skLWECoeffs {
		skLWECoeffs[i] = new(big.Int)
	}
	ringQLWE.PolyToBigintCentered(skLWE, 1, skLWECoeffs)
	for i := 0; i < int(n); i++ {
		sk_ss += float64(skLWECoeffs[i].Int64()) * float64(skLWECoeffs[i].Int64())
	}

	fmt.Printf("lwe sk ss : %f\n", sk_ss)

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)
	decLWE := lwe.NewDecoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testMVBRParamsLiteral[idx]) // N, Q
	require.NoError(b, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsBR).GenSecretKeyPairNew()

	N := paramsBR.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evk := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsBR).GenAutoKeysNew(galEl, f, f_inv, evkParams)...)
	BRK := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, f_inv, evkParams)

	// Encryptor, Decryptor for NTRU BR
	// encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	decryptorBR := ntru.NewDecryptor(paramsBR, f)

	// Evaluator for BR
	evalBR_Tree := NewTreeEvaluator(paramsBR, paramsLWE, evk)

	// TestPolynomial Used for rotation
	F_high, F_low := InitReLUPolynomial(p, q, ringQBR)

	//////////////////////////////////////////
	// KeySwitching Context, happens at QKS //
	//////////////////////////////////////////
	fmt.Printf("Key Switching Context Generation\n")
	paramsKSLWE := lwe.NewParameters(n, QKS, p)

	fCopy := ringQBR.NewPoly()
	ringQBR.INTT(f.Value.Q, fCopy)
	ringQBR.IMForm(fCopy, fCopy)

	encryptorKSLWE := lwe.NewEncryptor(paramsKSLWE, ringQLWE, skLWE)
	decryptorKSLWE := lwe.NewDecryptor(paramsKSLWE, ringQLWE, skLWE)

	errGenLKSK := lwe.NewErrorGenerator(stdev_lksk)

	lksk := lwe.LWE_KSKG(BKS, f, encryptorKSLWE, paramsBR.RingQ(), errGenLKSK) // Generate KeySwitching Key

	///////////////////////////////////
	// Setup finished, testing start //
	///////////////////////////////////
	fmt.Printf("Setup finished, start bootstrapping\n")

	var P = int(math.Pow(float64(p), float64(d)))
	fmt.Printf("Current Precision : %d\n", int(math.Log2(float64(p)))*d)

	var errs1 []float64 // make([]float64, tn)
	var errs2 []float64
	var Success = map[bool]string{true: "Success", false: "Failed"}

	var cnt = 0

	b.Run("NTRU Blind Rotation", func(b *testing.B) {

		for ti := 0; ti < tn; ti++ {
			var m = uint64(rand.Intn(P)) // message to encrypt/decrypt
			var m_decomp = Decomp(m, p, d)

			m_high := m_decomp[0]
			m_lows := m_decomp[1:]

			var m_high_relu = m_high
			if m_high >= uint64(p)/2 {
				m_high_relu = 0
			}

			// test using ReLU (0 <= x < P/2 -> x, P/2 <= x -> 0)
			ct_high := encryptorLWE.EncryptNew(encLWE.Encode(m_high, errGenLWE.GenErr()))

			// The high bits are processed using the same functions as single message
			// Along with it, the multi-value BR is calculated to get the low bit test polynomial
			ct_high_br, F_lows, err := evalBR_Tree.MultiValueBR(ct_high, &F_high, F_low, BRK) // ct_br is NTRU cipher
			require.NoError(b, err)

			// Error after BR
			pt_high_br := decryptorBR.DecryptNew(ct_high_br)
			ringQBR.INTT(pt_high_br.Value, pt_high_br.Value)

			if debug {
				exp_high_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_high_relu)
				err_high_br := int64(pt_high_br.Value.Coeffs[0][0]) - exp_high_br
				if err_high_br > int64(paramsBR.Q()[0]/2) {
					err_high_br = -int64(paramsBR.Q()[0]) + err_high_br
				}
				limit_high_br := float64(paramsBR.Q()[0]) / float64(4*p)
				fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_br, int64(pt_high_br.Value.Coeffs[0][0]))
				fmt.Printf("(H) BR Error    : %d/%d (%f%%)\n\n", err_high_br, int64(limit_high_br), math.Abs(float64(err_high_br))/limit_high_br*100.0)
			}

			ct_high_ks := lwe.LWE_KS(ct_high_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher

			if debug {
				pt_high_ks := decryptorKSLWE.DecryptNew(ct_high_ks)
				exp_high_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m_high_relu)
				err_high_ks := int64(pt_high_ks) - exp_high_ks
				if err_high_ks > int64(QKS/2) {
					err_high_ks = -int64(QKS) + err_high_ks
				}
				limit_high_ks := float64(QKS) / float64(4*p)
				fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_ks, int64(pt_high_ks))
				fmt.Printf("(H) MS-KS Error : %d/%d (%f%%)\n\n", err_high_ks, int64(limit_high_ks), math.Abs(float64(err_high_ks))/limit_high_ks*100.0)
			}

			ct_high_boot := lwe.LWE_MS(ct_high_ks, paramsLWE, QKS)
			pt_high_boot := decryptorLWE.DecryptNew(ct_high_boot)
			m_high_dec := decLWE.Decode(pt_high_boot)

			if debug {
				exp_high_boot := int64(uint64(math.Round(float64(q)/float64(2*p))) * m_high_relu)
				err_high_boot := int64(pt_high_boot) - exp_high_boot
				if err_high_boot > int64(q/2) {
					err_high_boot = -int64(q) + err_high_boot
				}
				limit_high_boot := float64(q) / float64(4*p)

				fmt.Printf("(H) Expected    : %020b\n(H) Actual      : %020b\n", exp_high_boot, int64(pt_high_boot))
				fmt.Printf("(H) MS Error    : %d/%d (%f%%)\n\n", err_high_boot, int64(limit_high_boot), math.Abs(float64(err_high_boot))/limit_high_boot*100.0)

				fmt.Printf("(H) %d -> %d\n\n", m_high, m_high_dec)
				fmt.Printf("=================================\n\n")
			}

			if m_high_dec == m_high_relu {
				err_calc := int64(pt_high_boot - (q / (2 * p) * m_high_relu))
				if err_calc > int64(q/2) {
					err_calc = err_calc - int64(q)
				}
				errs1 = append(errs1, float64(err_calc))
			}

			F_low_enc := evalBR_Tree.CombineTestPoly(F_lows, int(p), N)
			var m_dec = m_high_dec

			for _, m_low := range m_lows {
				m_low_relu := m_low
				if m_high >= uint64(p)/2 {
					m_low_relu = 0
				}
				ct_low := encryptorLWE.EncryptNew(encLWE.Encode(m_low, errGenLWE.GenErr()))

				// Using the key for encrypted TP, we perform the blind rotation for low bits
				ct_low_br, err := evalBR_Tree.EncBlindRotate(ct_low, F_low_enc, BRK)
				require.NoError(b, err)

				if debug {
					pt_low_br := decryptorBR.DecryptNew(ct_low_br)
					ringQBR.INTT(pt_low_br.Value, pt_low_br.Value)

					exp_low_br := int64(uint64(math.Round(float64(paramsBR.Q()[0])/float64(2*p))) * m_low_relu)
					err_low_br := int64(pt_low_br.Value.Coeffs[0][0]) - exp_low_br
					if err_low_br > int64(paramsBR.Q()[0]/2) {
						err_low_br = -int64(paramsBR.Q()[0]) + err_low_br
					}
					limit_low_br := float64(paramsBR.Q()[0]) / float64(4*p)

					fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_br, int64(pt_low_br.Value.Coeffs[0][0]))
					fmt.Printf("(L) BR Error    : %d/%d (%f%%)\n\n", err_low_br, int64(limit_low_br), math.Abs(float64(err_low_br))/limit_low_br*100.0)
				}

				ct_low_ks := lwe.LWE_KS(ct_low_br, BKS, encryptorKSLWE, paramsBR.RingQ(), lksk) // ct_ks is lwe cipher

				if debug {
					pt_low_ks := decryptorKSLWE.DecryptNew(ct_low_ks)
					exp_low_ks := int64(uint64(math.Round(float64(QKS)/float64(2*p))) * m_low_relu)
					err_low_ks := int64(pt_low_ks) - exp_low_ks
					if err_low_ks > int64(QKS/2) {
						err_low_ks = -int64(QKS) + err_low_ks
					}
					limit_low_ks := float64(QKS) / float64(4*p)

					fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_ks, int64(pt_low_ks))
					fmt.Printf("(L) MS-KS Error : %d/%d (%f%%)\n\n", err_low_ks, int64(limit_low_ks), math.Abs(float64(err_low_ks))/limit_low_ks*100.0)
				}

				ct_low_boot := lwe.LWE_MS(ct_low_ks, paramsLWE, QKS)
				pt_low_boot := decryptorLWE.DecryptNew(ct_low_boot)
				m_low_dec := decLWE.Decode(pt_low_boot)

				if debug {
					exp_low_boot := int64(uint64(math.Round(float64(q)/float64(2*p))) * m_low_relu)
					err_low_boot := int64(pt_low_boot) - exp_low_boot
					if err_low_boot > int64(q/2) {
						err_low_boot = -int64(q) + err_low_boot
					}
					limit_low_boot := float64(q) / float64(4*p)
					fmt.Printf("(L) Expected    : %020b\n(L) Actual      : %020b\n", exp_low_boot, int64(pt_low_boot))
					fmt.Printf("(L) MS Error    : %d/%d (%f%%)\n\n", err_low_boot, int64(limit_low_boot), math.Abs(float64(err_low_boot))/limit_low_boot*100.0)

					fmt.Printf("(L) %d -> %d\n\n", m_low, m_low_dec)
					fmt.Printf("=================================\n\n")
				}

				if m_low_dec == m_low_relu {
					err_calc := int64(pt_low_boot - (q / (2 * p) * m_low_relu))
					if err_calc > int64(q/2) {
						err_calc = err_calc - int64(q)
					}
					errs2 = append(errs2, float64(err_calc))
				}

				m_dec = m_dec*p + m_low_dec
			}

			if m >= uint64(P)/2 {
				// fmt.Printf("Bootstrap : %5d (>= %5d) -> %5d (%s)\n", m, uint64(P)/2, m_dec, Success[m_dec == 0])
				if m_dec == 0 {
					cnt++
				}
			} else {
				// fmt.Printf("Bootstrap : %5d (<  %5d) -> %5d (%s)\n", m, uint64(P)/2, m_dec, Success[m_dec == m])
				if m_dec == m {
					cnt++
				}
			}
		}

		em1, es1 := GetMeanStdev(errs1)
		fmt.Printf("\nError Stats : %f, %f\n", em1, es1)

		em2, es2 := GetMeanStdev(errs2)
		fmt.Printf("\nError Stats : %f, %f\n", em2, es2)

		// fmt.Printf("%v\n", norm_squares)
		nsm, _ := GetMeanStdev(norm_squares)
		fmt.Printf("\nNorm Square Mean : %f\n", nsm)

		fmt.Printf("Dec Success : %d/%d (%f%%)", cnt, tn, float64(cnt)/float64(tn)*100)
	})

	UNUSED(errs1, errs2, Success)
}

// Compare the speed of NTRU & RLWE blind rotation
func benchmarkNTRURLWE(b *testing.B, logp, logq, logB int, n uint64, Prob float64, tn, idx int) {
	prng, err := sampling.NewPRNG()
	require.NoError(b, err)

	// LWE params
	var p uint64 = 1 << logp // plaintext Modulus
	var q uint64 = 1 << logq // ciphertext Modulus

	var stdev_lwe = 3.2

	paramsLWE := lwe.NewParameters(n, q, p)

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(b, err)

	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: Prob}, false)
	require.NoError(b, err)
	skLWE := ternarySamplerLWE.ReadNew()

	skLWECoeffs := make([]*big.Int, ringQLWE.N())
	skLWEfloat := make([]float64, ringQLWE.N())
	for i := range skLWECoeffs {
		skLWECoeffs[i] = new(big.Int)
	}
	ringQLWE.PolyToBigintCentered(skLWE, 1, skLWECoeffs)
	var sk_ss float64 = 0
	for i := range skLWECoeffs {
		skLWEfloat[i] = float64(skLWECoeffs[i].Int64())
		if i < int(n) {
			sk_ss += skLWEfloat[i] * skLWEfloat[i]
		}
	}

	skvar, err := stats.StandardDeviation(skLWEfloat)
	if err != nil {
		panic(err)
	}
	fmt.Printf("lwe sk variance : %f\n", skvar)
	fmt.Printf("lwe sk ss : %f\n", sk_ss)

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)

	fmt.Printf("NTRU BR Context Generation\n")
	paramsNTRU, err := rlwe.NewParametersFromLiteral(testNTRUParamsLiteral[idx]) // N, Q
	require.NoError(b, err)

	// BR secret (f, f_inv) in NTT, Montgomery
	f, f_inv := ntru.NewKeyGenerator(paramsNTRU).GenSecretKeyPairNew()

	N := paramsNTRU.N()
	galEl := make([]uint64, 2*N-1)
	for i := 0; i < len(galEl); i++ {
		galEl[i] = 2*uint64(i) + 3
	}

	evkParamsNTRU := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsNTRU.MaxLevelQ()), LevelP: utils.Pointy(paramsNTRU.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	evkNTRU := rlwe.NewMemEvaluationKeySet(nil, ngsw.NewKeyGenerator(paramsNTRU).GenAutoKeysNew(galEl, f, f_inv, evkParamsNTRU)...)

	BRKNTRU := GenEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsNTRU, f_inv, evkParamsNTRU)
	evalBRNTRU := NewEvaluator(paramsNTRU, paramsLWE, evkNTRU)
	decNTRU := ntru.NewDecryptor(paramsNTRU, f)

	fmt.Printf("RLWE BR Context Generation")

	// RLWE params
	paramsRLWE, err := rlwe.NewParametersFromLiteral(testRLWEParamsLiteral[idx]) // N, Q
	require.NoError(b, err)

	skRLWE := rlwe.NewKeyGenerator(paramsRLWE).GenSecretKeyNew()

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsRLWE.MaxLevelQ()), LevelP: utils.Pointy(paramsRLWE.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	BRKRLWE, extrakey := GenRLWEEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsRLWE, skRLWE, evkParams)
	evalBRRLWE := NewRLWEEvaluator(paramsRLWE, paramsLWE)
	decRLWE := rlwe.NewDecryptor(paramsRLWE, skRLWE)

	// TestPolynomial Used for rotation
	F_ntru := InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsNTRU.RingQ())
	F_rlwe := InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsRLWE.RingQ())

	fmt.Printf("Setup finished, start Blind Rotation\n")

	// According to theory, NTRU should be about twice as fast as RLWE
	b.Run("NTRU Blind Rotation", func(b *testing.B) {
		var ntru_cnt = 0
		for ti := 0; ti < tn; ti++ {
			var m uint64 = uint64(rand.Intn(int(p))) // message to encrypt/decrypt
			ct := encryptorLWE.EncryptNew(encLWE.Encode(m, 0))

			ct_ntru_br, err := evalBRNTRU.BlindRotate(ct, &F_ntru, BRKNTRU) // ct_br is NTRU cipher
			require.NoError(b, err)

			pt_ntru_br := decNTRU.DecryptNew(ct_ntru_br)
			paramsNTRU.RingQ().INTT(pt_ntru_br.Value, pt_ntru_br.Value)

			m_br := (pt_ntru_br.Value.Coeffs[0][0] + uint64(float64(paramsNTRU.Q()[0])/float64(4*p))) % paramsNTRU.Q()[0]
			// fmt.Printf("NTRU Bootstrap : %d -> %d\n", m, int(float64(m_br*2*p)/float64(paramsNTRU.Q()[0])))
			if m == uint64(float64(m_br*2*p)/float64(paramsNTRU.Q()[0])) {
				ntru_cnt++
			}
		}
		fmt.Printf("NTRU BR Success : %d/%d (%f%%)\n", ntru_cnt, tn, float64(ntru_cnt)/float64(tn)*100.0)
	})

	b.Run("RLWE Blind Rotation", func(b *testing.B) {
		var rlwe_cnt = 0
		for ti := 0; ti < tn; ti++ {
			var m uint64 = uint64(rand.Intn(int(p))) // message to encrypt/decrypt
			ct := encryptorLWE.EncryptNew(encLWE.Encode(m, errGenLWE.GenErr()))

			// Do BlindRotation
			ct_br, err := evalBRRLWE.BlindRotate(ct, &F_rlwe, BRKRLWE, extrakey) // ct_br is NTRU cipher
			require.NoError(b, err)

			pt_rlwe_br := decRLWE.DecryptNew(ct_br)
			paramsRLWE.RingQ().INTT(pt_rlwe_br.Value, pt_rlwe_br.Value)

			m_br := (pt_rlwe_br.Value.Coeffs[0][0] + uint64(float64(paramsRLWE.Q()[0])/float64(4*p))) % paramsRLWE.Q()[0]
			// fmt.Printf("RLWE Bootstrap : %d -> %d\n", m, int(float64(m_br*2*p)/float64(paramsRLWE.Q()[0])))
			if m == uint64(float64(m_br*2*p)/float64(paramsRLWE.Q()[0])) {
				rlwe_cnt++
			}
		}
		fmt.Printf("RLWE BR Success : %d/%d (%f%%)\n", rlwe_cnt, tn, float64(rlwe_cnt)/float64(tn)*100.0)
	})
}

/*

func TestRLWEBR(t *testing.T) {
	prng, err := sampling.NewPRNG()
	require.NoError(t, err)

	// parameters of the TFHE sample
	var p uint64 = 1 << 3  // plaintext Modulus
	var q uint64 = 1 << 10 // ciphertext Modulus
	var n uint64 = 512     // TFHE dimension

	var logB = 2 // Decomposition Base used for Blind Rotation B = 2^logB
	// var BKS = 1 << 7 // Decomposition Base used for Key Switching

	var stdev_lwe = 3.2
	var bound_lwe = int64(5)
	// var stdev_lksk = 3.2

	var tn = 1

	/////////////////
	// LWE Context //
	/////////////////
	fmt.Printf("LWE Context Generation\n")
	paramsLWE := lwe.NewParameters(n, q, p) // n, q, p

	ringQLWE_n := 1 << int(math.Ceil(math.Log2(float64(n))))
	ringQLWE, err := ring.NewRing(ringQLWE_n, []uint64{12289})
	require.NoError(t, err)

	// generate skLWE
	ternarySamplerLWE, err := ring.NewSampler(prng, ringQLWE, ring.Ternary{P: 1.0 / 2.0}, false)
	require.NoError(t, err)
	skLWE := ternarySamplerLWE.ReadNew()

	// Error Generator
	errGenLWE := lwe.NewErrorGenerator(stdev_lwe, bound_lwe)

	// Encoder, Decoder for LWE
	encLWE := lwe.NewEncoder(paramsLWE)

	// Encryptor, Decryptor for LWE
	encryptorLWE := lwe.NewEncryptor(paramsLWE, ringQLWE, skLWE)
	// decryptorLWE := lwe.NewDecryptor(paramsLWE, ringQLWE, skLWE)

	//////////////////////////////////
	// BlindRotation Context (NTRU) //
	//////////////////////////////////
	fmt.Printf("BR Context Generation\n")
	paramsBR, err := rlwe.NewParametersFromLiteral(testRLWEParamsLiteral[0]) // N, Q
	require.NoError(t, err)
	ringQBR := paramsBR.RingQ()

	// BR secret (f, f_inv) in NTT, Montgomery
	skRLWE := rlwe.NewKeyGenerator(paramsBR).GenSecretKeyNew()

	evkParams := rlwe.EvaluationKeyParameters{LevelQ: utils.Pointy(paramsBR.MaxLevelQ()), LevelP: utils.Pointy(paramsBR.MaxLevelP()), BaseTwoDecomposition: utils.Pointy(logB)}
	BRK, extrakey := GenRLWEEvaluationKeyNew(ringQLWE, skLWE, paramsLWE, paramsBR, skRLWE, evkParams)

	// Encryptor, Decryptor for NTRU BR
	decryptorBR := rlwe.NewDecryptor(paramsBR, skRLWE)
	encryptorBR := rlwe.NewEncryptor(paramsBR, skRLWE)

	// encryptorBR := ntru.NewEncryptor(paramsBR, f_inv)
	// decryptorBR := ntru.NewDecryptor(paramsBR, f)

	// Evaluator for BR
	evalBR := NewRLWEEvaluator(paramsBR, paramsLWE)

	// TestPolynomial Used for rotation
	// F := InitTestPolynomialDebug(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ()) // InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ())
	F := InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ()) // InitTestPolynomial(paramsLWE.P(), paramsLWE.Q(), paramsBR.RingQ())

	///////////////////////////////////
	// Setup finished, testing start //
	///////////////////////////////////
	fmt.Printf("Setup finished, start bootstrapping\n")

	for ti := 0; ti < tn; ti++ {

		var m uint64 = uint64(rand.Intn(int(p))) // message to encrypt/decrypt
		ct := encryptorLWE.EncryptNew(encLWE.Encode(m, errGenLWE.GenErr()))
		// ct := encryptorLWE.EncryptNew(encLWE.Encode(m, 0))

		// Do BlindRotation
		ct_br, err := evalBR.BlindRotate(ct, &F, BRK, extrakey) // ct_br is NTRU cipher
		require.NoError(t, err)

		pt_br := decryptorBR.DecryptNew(ct_br)
		pt_br_poly := pt_br.Value

		ringQBR.INTT(pt_br_poly, pt_br_poly)

		m_br := (pt_br_poly.Coeffs[0][0] + uint64(float64(paramsBR.Q()[0])/float64(4*p))) % paramsBR.Q()[0]
		fmt.Printf("Bootstrap : %d -> %d\n", m, int(float64(m_br*2*p)/float64(paramsBR.Q()[0])))
	}

	UNUSED(ringQBR, tn, errGenLWE, encLWE, encryptorLWE, decryptorBR, evalBR, F, BRK, extrakey, encryptorBR)
}
*/

package hpntru

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/rgsw"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

const (
	// Parameter w of Algorithm 3 in https://eprint.iacr.org/2022/198
	windowSize = 5
)

type RLWEBlindRotationEvaluationKeySet interface {

	// GetBlindRotationKey should return RGSW(X^{s[i]})
	GetBlindRotationKey(i int) (brk *rgsw.Ciphertext, err error)

	// GetEvaluationKeySet should return an rlwe.EvaluationKeySet
	// providing access to all the required automorphism keys.
	GetEvaluationKeySet() (evk rlwe.EvaluationKeySet, err error)
}

// MemBlindRotationEvaluationKeySet is a basic in-memory implementation of the BlindRotationEvaluationKeySet interface.
type RLWEMemBlindRotationEvaluationKeySet struct {
	BlindRotationKeys []*rgsw.Ciphertext
	AutomorphismKeys  []*rlwe.GaloisKey
}

func (evk RLWEMemBlindRotationEvaluationKeySet) GetBlindRotationKey(i int) (*rgsw.Ciphertext, error) {
	return evk.BlindRotationKeys[i], nil
}

func (evk RLWEMemBlindRotationEvaluationKeySet) GetEvaluationKeySet() (rlwe.EvaluationKeySet, error) {
	return rlwe.NewMemEvaluationKeySet(nil, evk.AutomorphismKeys...), nil
}

// GenEvaluationKeyNew generates a new Blind Rotation evaluation key
func GenRLWEEvaluationKeyNew(ringQ *ring.Ring, skLWE ring.Poly, paramsLWE lwe.Parameters, paramsBR rlwe.Parameters, skRLWE *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (key RLWEMemBlindRotationEvaluationKeySet, extrakey *rgsw.Ciphertext) {

	skLWECopy := ringQ.NewPoly()
	skLWECopy.Copy(skLWE)
	sk := make([]*big.Int, skLWECopy.N())

	for i := range sk {
		sk[i] = new(big.Int)
	}
	ringQ.AtLevel(0).PolyToBigintCentered(skLWECopy, 1, sk)

	for i := int(paramsLWE.N()); i < skLWECopy.N(); i++ {
		sk[i] = new(big.Int).SetInt64(0)
	}

	encryptor := rgsw.NewEncryptor(paramsBR, skRLWE)
	levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(paramsBR, evkParams)

	// Generate upto n keys
	skiRGSW := make([]*rgsw.Ciphertext, paramsLWE.N())
	var sk_sum int = 0 // := 2 * paramsBR.N() * paramsLWEKG.N()
	ptXi := make(map[int]*rlwe.Plaintext)

	for i := 0; i < int(paramsLWE.N()); i++ {
		siInt := int(sk[i].Int64())
		sk_sum = sk_sum + siInt

		if _, ok := ptXi[siInt]; !ok {

			pt := &rlwe.Plaintext{}
			pt.MetaData = &rlwe.MetaData{}
			pt.IsNTT = true
			pt.IsMontgomery = true
			pt.Value = paramsBR.RingQ().NewMonomialXi(siInt)
			// fmt.Printf("%v\n", pt.Value.Coeffs[0])
			paramsBR.RingQ().NTT(pt.Value, pt.Value)
			paramsBR.RingQ().MForm(pt.Value, pt.Value)

			ptXi[siInt] = pt
		}

		skiRGSW[i] = rgsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)

		// Sanity check, this error should never happen unless this algorithm
		// has been improperly modified to provides invalid inputs.

		if err := encryptor.Encrypt(ptXi[siInt], skiRGSW[i]); err != nil {
			panic(err)
		}
	}

	// twoN := (2 * paramsBR.N())
	// sk_sum = ((sk_sum % twoN) + twoN) % twoN
	// fmt.Printf("Sum : %d\n", -sk_sum)

	pt := &rlwe.Plaintext{}
	pt.MetaData = &rlwe.MetaData{}
	pt.IsNTT = true
	pt.IsMontgomery = true
	pt.Value = paramsBR.RingQ().NewMonomialXi(-sk_sum)
	paramsBR.RingQ().NTT(pt.Value, pt.Value)
	paramsBR.RingQ().MForm(pt.Value, pt.Value)

	extrakey = rgsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)
	if err := encryptor.Encrypt(pt, extrakey); err != nil {
		panic(err)
	}

	kgen := rlwe.NewKeyGenerator(paramsBR)

	galEls := make([]uint64, windowSize)
	for i := 0; i < windowSize; i++ {
		galEls[i] = paramsBR.GaloisElement(i + 1)
		// fmt.Printf("%d\n", galEls[i])
	}

	galEls = append(galEls, paramsBR.RingQ().NthRoot()-ring.GaloisGen)

	gks := kgen.GenGaloisKeysNew(galEls, skRLWE, rlwe.EvaluationKeyParameters{
		LevelQ:               utils.Pointy(levelQ),
		LevelP:               utils.Pointy(levelP),
		BaseTwoDecomposition: utils.Pointy(BaseTwoDecomposition),
	})

	return RLWEMemBlindRotationEvaluationKeySet{BlindRotationKeys: skiRGSW, AutomorphismKeys: gks}, extrakey
}

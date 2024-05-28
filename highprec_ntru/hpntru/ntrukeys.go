package hpntru

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/lwe"
	"github.com/tuneinsight/lattigo/v5/core/ngsw"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
)

type BlindRotationEvaluationKeySet interface {

	// GetBlindRotationKey should return RGSW(X^{s[i]})
	GetBlindRotationKey(i int) (brk *ngsw.Ciphertext, err error)

	// GetEvaluationKeySet should return an rlwe.EvaluationKeySet
	// providing access to all the required automorphism keys.
	GetEvaluationKeySet() (evk rlwe.EvaluationKeySet, err error)
}

// MemBlindRotationEvaluationKeySet is a basic in-memory implementation of the BlindRotationEvaluationKeySet interface.
type MemBlindRotationEvaluationKeySet struct {
	BlindRotationKeys []*ngsw.Ciphertext
	AutomorphismKeys  []*rlwe.GaloisKey
}

func (evk MemBlindRotationEvaluationKeySet) GetBlindRotationKey(i int) (*ngsw.Ciphertext, error) {
	return evk.BlindRotationKeys[i], nil
}

func (evk MemBlindRotationEvaluationKeySet) GetEvaluationKeySet() (rlwe.EvaluationKeySet, error) {
	return rlwe.NewMemEvaluationKeySet(nil, evk.AutomorphismKeys...), nil
}

// GenEvaluationKeyNew generates a new Blind Rotation evaluation key
func GenEvaluationKeyNew(ringQ *ring.Ring, skLWE ring.Poly, paramsLWE lwe.Parameters, paramsBR rlwe.Parameters, f_inv *rlwe.SecretKey, evkParams ...rlwe.EvaluationKeyParameters) (key MemBlindRotationEvaluationKeySet) {

	skLWECopy := ringQ.NewPoly()
	skLWECopy.Copy(skLWE)
	sk := make([]*big.Int, skLWECopy.N())
	for i := range sk {
		sk[i] = new(big.Int)
	}
	ringQ.AtLevel(0).PolyToBigintCentered(skLWECopy, 1, sk)

	encryptor := ngsw.NewEncryptor(paramsBR, f_inv)
	levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(paramsBR, evkParams)

	// Generate upto n + 1 keys
	skiNGSW := make([]*ngsw.Ciphertext, paramsLWE.N()+2)
	var sk_sum int = 0 // := 2 * paramsBR.N() * paramsLWEKG.N()
	ptXi := make(map[int]*rlwe.Plaintext)

	for i := 0; i < int(paramsLWE.N()); i++ {
		siInt := int(sk[i].Int64())
		sk_sum -= siInt

		if i == 0 {
			pt := &rlwe.Plaintext{}
			pt.MetaData = &rlwe.MetaData{}
			pt.IsNTT = true
			// pt.IsMontgomery = true

			tmp := paramsBR.RingQ().NewPoly()
			tmp.Copy(f_inv.Value.Q)
			paramsBR.RingQ().INTT(tmp, tmp)
			paramsBR.RingQ().IMForm(tmp, tmp)

			pt.Value = MulByMonomial(paramsBR.RingQ(), tmp, siInt)
			paramsBR.RingQ().NTT(pt.Value, pt.Value)
			// paramsBR.RingQ().MForm(pt.Value, pt.Value)

			skiNGSW[i] = ngsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)

			// Sanity check, this error should never happen unless this algorithm
			// has been improperly modified to provides invalid inputs.

			if err := encryptor.Encrypt(pt, skiNGSW[i]); err != nil {
				panic(err)
			}
		} else {
			if _, ok := ptXi[siInt]; !ok {

				pt := &rlwe.Plaintext{}
				pt.MetaData = &rlwe.MetaData{}
				pt.IsNTT = true
				pt.IsMontgomery = true
				pt.Value = paramsBR.RingQ().NewMonomialXi(siInt)
				paramsBR.RingQ().NTT(pt.Value, pt.Value)
				paramsBR.RingQ().MForm(pt.Value, pt.Value)

				ptXi[siInt] = pt
			}

			skiNGSW[i] = ngsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)

			// Sanity check, this error should never happen unless this algorithm
			// has been improperly modified to provides invalid inputs.

			if err := encryptor.Encrypt(ptXi[siInt], skiNGSW[i]); err != nil {
				panic(err)
			}
		}
	}

	twoN := (2 * paramsBR.N())
	sk_sum = ((sk_sum % twoN) + twoN) % twoN

	pt := &rlwe.Plaintext{}
	pt.MetaData = &rlwe.MetaData{}
	pt.IsNTT = true
	pt.IsMontgomery = true
	pt.Value = paramsBR.RingQ().NewMonomialXi(sk_sum)
	paramsBR.RingQ().NTT(pt.Value, pt.Value)
	paramsBR.RingQ().MForm(pt.Value, pt.Value)

	skiNGSW[paramsLWE.N()] = ngsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)
	if err := encryptor.Encrypt(pt, skiNGSW[paramsLWE.N()]); err != nil {
		panic(err)
	}

	// final skiNGSW for X^si[0] (for encrypted TP)
	siInt := int(sk[0].Int64())
	if _, ok := ptXi[siInt]; !ok {

		pt := &rlwe.Plaintext{}
		pt.MetaData = &rlwe.MetaData{}
		pt.IsNTT = true
		pt.IsMontgomery = true
		pt.Value = paramsBR.RingQ().NewMonomialXi(siInt)
		paramsBR.RingQ().NTT(pt.Value, pt.Value)
		paramsBR.RingQ().MForm(pt.Value, pt.Value)

		ptXi[siInt] = pt
	}
	skiNGSW[paramsLWE.N()+1] = ngsw.NewCiphertext(paramsBR, levelQ, levelP, BaseTwoDecomposition)

	if err := encryptor.Encrypt(ptXi[siInt], skiNGSW[paramsLWE.N()+1]); err != nil {
		panic(err)
	}

	return MemBlindRotationEvaluationKeySet{BlindRotationKeys: skiNGSW, AutomorphismKeys: nil}

	/*


		pRLWE := *paramsRLWE.GetRLWEParameters()
		pLWE := *paramsLWE.GetRLWEParameters()

		skLWECopy := skLWE.CopyNew()
		pLWE.RingQ().AtLevel(0).INTT(skLWECopy.Value.Q, skLWECopy.Value.Q)
		pLWE.RingQ().AtLevel(0).IMForm(skLWECopy.Value.Q, skLWECopy.Value.Q)
		sk := make([]*big.Int, pLWE.N())
		for i := range sk {
			sk[i] = new(big.Int)
		}
		pLWE.RingQ().AtLevel(0).PolyToBigintCentered(skLWECopy.Value.Q, 1, sk)

		encryptor := rgsw.NewEncryptor(pRLWE, skRLWE)

		levelQ, levelP, BaseTwoDecomposition := rlwe.ResolveEvaluationKeyParameters(pRLWE, evkParams)

		skiRGSW := make([]*rgsw.Ciphertext, pLWE.N())

		ptXi := make(map[int]*rlwe.Plaintext)

		for i, si := range sk {

			siInt := int(si.Int64())

			if _, ok := ptXi[siInt]; !ok {

				pt := &rlwe.Plaintext{}
				pt.MetaData = &rlwe.MetaData{}
				pt.IsNTT = true
				pt.Value = pRLWE.RingQ().NewMonomialXi(siInt)
				pRLWE.RingQ().NTT(pt.Value, pt.Value)

				ptXi[siInt] = pt
			}

			skiRGSW[i] = rgsw.NewCiphertext(pRLWE, levelQ, levelP, BaseTwoDecomposition)

			// Sanity check, this error should never happen unless this algorithm
			// has been improperly modified to provides invalid inputs.
			if err := encryptor.Encrypt(ptXi[siInt], skiRGSW[i]); err != nil {
				panic(err)
			}
		}

		kgen := rlwe.NewKeyGenerator(pRLWE)

		galEls := make([]uint64, windowSize)
		for i := 0; i < windowSize; i++ {
			galEls[i] = pRLWE.GaloisElement(i + 1)
		}

		galEls = append(galEls, pRLWE.RingQ().NthRoot()-ring.GaloisGen)

		gks := kgen.GenGaloisKeysNew(galEls, skRLWE, rlwe.EvaluationKeyParameters{
			LevelQ:               utils.Pointy(levelQ),
			LevelP:               utils.Pointy(levelP),
			BaseTwoDecomposition: utils.Pointy(BaseTwoDecomposition),
		})

		return MemBlindRotationEvaluationKeySet{BlindRotationKeys: skiRGSW, AutomorphismKeys: gks}
	*/
}

package ngsw

import (
	"github.com/tuneinsight/lattigo/v5/core/ntru"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils"
)

type Encryptor struct {
	*ntru.Encryptor
	buffQP ringqp.Poly
}

func NewEncryptor(params rlwe.Parameters, key rlwe.EncryptionKey) *Encryptor {
	return &Encryptor{ntru.NewEncryptor(params, key), params.GetRLWEParameters().RingQP().NewPoly()}
}

func (enc Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) (err error) {
	var ntruvCt *Ciphertext
	var isNTRUV bool
	if ntruvCt, isNTRUV = ct.(*Ciphertext); !isNTRUV {
		return enc.Encryptor.Encrypt(pt, ct)
	}

	if err = enc.EncryptZero(ntruvCt); err != nil {
		return
	}

	params := enc.GetRLWEParameters()
	levelQ := ntruvCt.LevelQ()
	ringQ := params.RingQ().AtLevel(levelQ)

	if pt != nil {
		if !pt.IsNTT {
			ringQ.NTT(pt.Value, enc.buffQP.Q)

			if !pt.IsMontgomery {
				ringQ.MForm(enc.buffQP.Q, enc.buffQP.Q)
			}

		} else {
			if !pt.IsMontgomery {
				ringQ.MForm(pt.Value, enc.buffQP.Q)
			} else {
				enc.buffQP.Q.CopyLvl(levelQ, pt.Value)
			}
		}

		if err := AddPolyTimesGadgetVectorToGadgetCiphertext(
			enc.buffQP.Q,
			ntruvCt,
			*params.RingQP(),
			enc.buffQP.Q); err != nil {
			// Sanity check, this error should not happen.
			panic(err)
		}
	}

	return nil
}

func (enc Encryptor) EncryptZero(ct interface{}) (err error) {

	var ntruvCt *Ciphertext
	var isNTRUV bool
	if ntruvCt, isNTRUV = ct.(*Ciphertext); !isNTRUV {
		return enc.Encryptor.EncryptZero(ntruvCt)
	}

	BaseRNSDecompositionVectorSize := ntruvCt.Value.BaseRNSDecompositionVectorSize()
	BaseTwoDecompositionVectorSize := ntruvCt.Value.BaseTwoDecompositionVectorSize()

	metadata := &rlwe.MetaData{}
	metadata.IsMontgomery = true
	metadata.IsNTT = true

	for i := 0; i < BaseRNSDecompositionVectorSize; i++ {
		for j := 0; j < BaseTwoDecompositionVectorSize[i]; j++ {
			if err = enc.Encryptor.EncryptZero(rlwe.Element[ringqp.Poly]{MetaData: metadata, Value: []ringqp.Poly(ntruvCt.Value.Value[i][j])}); err != nil {
				return
			}
		}
	}

	return nil
}

// ShallowCopy creates a shallow copy of this Encryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc Encryptor) ShallowCopy() *Encryptor {
	return &Encryptor{Encryptor: enc.Encryptor.ShallowCopy(), buffQP: enc.GetRLWEParameters().RingQP().NewPoly()}
}

func AddPolyTimesGadgetVectorToGadgetCiphertext(pt ring.Poly, ct1 *Ciphertext, ringQP ringqp.Ring, buff ring.Poly) (err error) {
	ct := ct1.Value

	levelQ := ct.LevelQ()
	levelP := ct.LevelP()

	ringQ := ringQP.RingQ.AtLevel(levelQ)

	if levelP != -1 {
		ringQ.MulScalarBigint(pt, ringQP.RingP.AtLevel(levelP).Modulus(), buff) // P * pt
	} else {
		levelP = 0
		buff.CopyLvl(levelQ, pt) // 1 * pt
	}

	BaseRNSDecompositionVectorSize := len(ct.Value)

	BaseTwoDecompositionVectorSize := make([]int, len(ct.Value))
	for i := range BaseTwoDecompositionVectorSize {
		BaseTwoDecompositionVectorSize[i] = len(ct.Value[i])
	}

	N := ringQ.N()

	var index int
	for j := 0; j < utils.MaxSlice(BaseTwoDecompositionVectorSize); j++ {

		for i := 0; i < BaseRNSDecompositionVectorSize; i++ {

			if j < BaseTwoDecompositionVectorSize[i] {

				// e + (m * P * w^2j) * (q_star * q_tild) mod QP
				//
				// q_prod = prod(q[i*#Pi+j])
				// q_star = Q/qprod
				// q_tild = q_star^-1 mod q_prod
				//
				// Therefore : (pt * P * w^2j) * (q_star * q_tild) = pt*P*w^2j mod q[i*#Pi+j], else 0
				for k := 0; k < levelP+1; k++ {

					index = i*(levelP+1) + k

					// Handle cases where #pj does not divide #qi
					if index >= levelQ+1 {
						break
					}

					qi := ringQ.SubRings[index].Modulus
					p0tmp := buff.Coeffs[index]
					p1tmp := ct.Value[i][j][1].Q.Coeffs[index]
					for w := 0; w < N; w++ {
						p1tmp[w] = ring.CRed(p1tmp[w]+p0tmp[w], qi)
					}
				}
			}
		}

		// w^2j
		ringQ.MulScalar(buff, 1<<ct.BaseTwoDecomposition, buff)
	}

	return
}

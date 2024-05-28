package ntru

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

type Decryptor struct {
	params rlwe.Parameters
	ringQ  *ring.Ring
	buff   ring.Poly
	sk     *rlwe.SecretKey // will save the polynomial f
}

func NewDecryptor(params rlwe.ParameterProvider, sk *rlwe.SecretKey) *Decryptor {

	p := params.GetRLWEParameters()

	if sk.Value.Q.N() != p.N() {
		panic(fmt.Errorf("cannot NewDecryptor: secret_key ring degree does not match parameters ring degree"))
	}

	return &Decryptor{
		params: *p,
		ringQ:  p.RingQ(),
		buff:   p.RingQ().NewPoly(),
		sk:     sk,
	}
}

// GetRLWEParameters returns the underlying rlwe.Parameters.
func (d Decryptor) GetRLWEParameters() *rlwe.Parameters {
	return &d.params
}

// DecryptNew decrypts the Ciphertext and returns the result in a new Plaintext.
// Output pt MetaData will match the input ct MetaData.
func (d Decryptor) DecryptNew(ct *rlwe.Ciphertext) (pt *rlwe.Plaintext) {
	pt = rlwe.NewPlaintext(d.params, ct.Level())
	d.Decrypt(ct, pt)
	return
}

// ciphertext should be in Montgomery form
// NTT is resolved
func (d Decryptor) Decrypt(ct *rlwe.Ciphertext, pt *rlwe.Plaintext) {

	level := utils.Min(ct.Level(), pt.Level())

	ringQ := d.ringQ.AtLevel(level)

	pt.Resize(0, level)

	*pt.MetaData = *ct.MetaData

	if ct.IsNTT {
		pt.Value.CopyLvl(level, ct.Value[ct.Degree()])
	} else {
		ringQ.NTTLazy(ct.Value[ct.Degree()], pt.Value)
	}

	for i := ct.Degree(); i > 0; i-- {

		ringQ.MulCoeffsMontgomery(pt.Value, d.sk.Value.Q, pt.Value)

		if !ct.IsNTT {
			ringQ.NTTLazy(ct.Value[i-1], d.buff)
			ringQ.Add(pt.Value, d.buff, pt.Value)
		} else {
			ringQ.Add(pt.Value, ct.Value[i-1], pt.Value)
		}

		if i&7 == 7 {
			ringQ.Reduce(pt.Value, pt.Value)
		}
	}

	if (ct.Degree())&7 != 7 {
		ringQ.Reduce(pt.Value, pt.Value)
	}

	if !ct.IsNTT {
		ringQ.INTT(pt.Value, pt.Value)
	}
}

// ShallowCopy creates a shallow copy of Decryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Decryptor can be used concurrently.
func (d Decryptor) ShallowCopy() *Decryptor {
	return &Decryptor{
		params: d.params,
		ringQ:  d.ringQ,
		buff:   d.ringQ.NewPoly(),
		sk:     d.sk,
	}
}

// WithKey creates a shallow copy of Decryptor with a new decryption key, in which all the
// read-only data-structures are shared with the receiver and the temporary buffers
// are reallocated. The receiver and the returned Decryptor can be used concurrently.
func (d Decryptor) WithKey(sk *rlwe.SecretKey) *Decryptor {
	return &Decryptor{
		params: d.params,
		ringQ:  d.ringQ,
		buff:   d.ringQ.NewPoly(),
		sk:     sk,
	}
}

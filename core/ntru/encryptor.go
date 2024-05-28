package ntru

/*
ciphertext format will be same as rlwe,
but instead of (-as + m + e, a), the format (0, g*f_inv + m*f_inv) will be placed
encryptor should save the inverse key and the secret key generated from KeyGen

Decryption will take place same as RLWE, as (-as + m + e) + a * s = m + e
0 + (g*f_inv + m*f_inv) * f = m + g
*/

import (
	"fmt"
	"math/big"
	"reflect"

	"github.com/montanaflynn/stats"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

type Encryptor struct {
	params rlwe.Parameters
	*encryptorBuffers

	encKey         rlwe.EncryptionKey
	prng           sampling.PRNG
	xeSampler      ring.Sampler
	xsSampler      ring.Sampler
	tmpSampler     ring.Sampler
	basisextender  *ring.BasisExtender
	uniformSampler ringqp.UniformSampler
}

type encryptorBuffers struct {
	buffQ  [2]ring.Poly
	buffP  [3]ring.Poly
	buffQP ringqp.Poly
}

// GetRLWEParameters returns the underlying rlwe.Parameters.
func (enc Encryptor) GetRLWEParameters() *rlwe.Parameters {
	return &enc.params
}

// NewEncryptor creates a new Encryptor only accepts secret key
func NewEncryptor(params rlwe.ParameterProvider, key rlwe.EncryptionKey) *Encryptor {

	p := *params.GetRLWEParameters()

	enc := newEncryptor(p)
	var err error
	switch key := key.(type) {
	case *rlwe.SecretKey:
		err = enc.checkSk(key)
	case nil:
		return newEncryptor(p)
	default:
		// Sanity check
		panic(fmt.Errorf("key must be *rlwe.SecretKey but have %T", key))
	}

	if err != nil {
		// Sanity check, this error should not happen.
		panic(fmt.Errorf("key is not correct: %w", err))
	}

	enc.encKey = key
	return enc
}

func newEncryptor(params rlwe.Parameters) *Encryptor {

	prng, err := sampling.NewPRNG()
	if err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	// xeSampler, err := ring.NewSampler(prng, params.RingQ(), ring.Ternary{H: 4}, false)
	// fmt.Printf("%s\n", params.Xe())
	xeSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xe(), false)
	// xeSampler, err := ring.NewSampler(prng, params.RingQ(), , false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(fmt.Errorf("newEncryptor: %w", err))
	}

	xsSampler, err := ring.NewSampler(prng, params.RingQ(), params.Xs(), false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(fmt.Errorf("newEncryptor: %w", err))
	}

	tmpSampler, err := ring.NewSampler(prng, params.RingQ(), ring.Ternary{H: 1}, false)

	// Sanity check, this error should not happen.
	if err != nil {
		panic(fmt.Errorf("newEncryptor: %w", err))
	}

	return &Encryptor{
		params:           params,
		prng:             prng,
		xeSampler:        xeSampler,
		xsSampler:        xsSampler,
		tmpSampler:       tmpSampler,
		encryptorBuffers: newEncryptorBuffers(params),
		uniformSampler:   ringqp.NewUniformSampler(prng, *params.RingQP()),
		basisextender:    bc,
	}
}

func newEncryptorBuffers(params rlwe.Parameters) *encryptorBuffers {

	ringQ := params.RingQ()
	ringP := params.RingP()

	var buffP [3]ring.Poly
	if params.PCount() != 0 {
		buffP = [3]ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return &encryptorBuffers{
		buffQ:  [2]ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		buffP:  buffP,
		buffQP: params.RingQP().NewPoly(),
	}
}

// used for only type change from plaintext to ciphertext
func (enc Encryptor) EncryptPoly(pt *rlwe.Plaintext, ct interface{}) (err error) {
	if pt == nil {
		return enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *rlwe.Ciphertext:
			*ct.MetaData = *pt.MetaData
			level := utils.Min(pt.Level(), ct.Level())
			ct.Resize(ct.Degree(), level)
			enc.addPtToCt(level, pt, ct) // Difference with original
			return
		default:
			return fmt.Errorf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct))
		}
	}
}

func (enc Encryptor) Encrypt(pt *rlwe.Plaintext, ct interface{}) (err error) {
	if pt == nil {
		return enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *rlwe.Ciphertext:
			*ct.MetaData = *pt.MetaData
			level := utils.Min(pt.Level(), ct.Level())
			ct.Resize(ct.Degree(), level)
			if err = enc.EncryptZero(ct); err != nil {
				return fmt.Errorf("cannot Encrypt: %w", err)
			}
			enc.addPtDivfToCt(level, pt, ct) // Difference with original
			return
		default:
			return fmt.Errorf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct))
		}
	}
}

func (enc Encryptor) EncryptNew(pt *rlwe.Plaintext) (ct *rlwe.Ciphertext, err error) {
	ct = rlwe.NewCiphertext(enc.params, 1, pt.Level())
	return ct, enc.Encrypt(pt, ct)
}

// Only utilize one of the two polynomial terms of rlwe ciphertext structure
// encryptor only has secret key (to be exact, inverse of secret key)
func (enc Encryptor) EncryptZero(ct interface{}) (err error) {
	switch key := enc.encKey.(type) {
	case *rlwe.SecretKey:
		return enc.encryptZeroSk(key, ct)
	default:
		return fmt.Errorf("cannot encrypt: Encryptor has no encryption key")
	}
}

func (enc Encryptor) EncryptZeroNew(level int) (ct *rlwe.Ciphertext) {
	ct = rlwe.NewCiphertext(enc.params, 1, level)
	if err := enc.EncryptZero(ct); err != nil {
		// Sanity check, this error should not happen.
		panic(err)
	}
	return
}

///////////////////////////////////////////////// Parts to modify //////////////////////////////////////////////////////

func (enc Encryptor) encryptZeroSk(sk *rlwe.SecretKey, ct interface{}) (err error) {
	switch ct := ct.(type) {
	case *rlwe.Ciphertext:
		return enc.encryptZeroSkFromC1(sk, ct.Element)

	case rlwe.Element[ringqp.Poly]:
		return enc.encryptZeroSkFromC1QP(sk, ct)

	default:
		return fmt.Errorf("cannot EncryptZero: input ciphertext type %T is not supported", ct)
	}
}

// consider that the result is in montgomery domain
func (enc Encryptor) encryptZeroSkFromC1(sk *rlwe.SecretKey, ct rlwe.Element[ring.Poly]) (err error) {

	levelQ := ct.Level()
	ringQ := enc.params.RingQ().AtLevel(levelQ)

	c1 := ct.Value[1]
	enc.xeSampler.AtLevel(levelQ).Read(c1) // generating g

	errCoeffs := make([]*big.Int, enc.params.N())
	errfloat := make([]float64, enc.params.N())
	for i := range errCoeffs {
		errCoeffs[i] = new(big.Int)
	}
	ringQ.PolyToBigintCentered(c1, 1, errCoeffs)
	for i := range errCoeffs {
		errfloat[i] = float64(errCoeffs[i].Int64())
	}

	skvar, err := stats.StandardDeviation(errfloat)
	if err != nil {
		panic(err)
	}
	if true {
		fmt.Printf("NTRU err variance : %f\n", skvar)
	}

	ringQ.NTT(c1, c1)
	ringQ.MForm(c1, c1)
	ringQ.MulCoeffsMontgomery(c1, sk.Value.Q, c1) // inverse of secret key multiplied, c0 = NTT(g * f_inv)

	if !ct.IsNTT {
		fmt.Printf("Not in NTT C1!\n")
		ringQ.INTT(c1, c1)
	}

	return
}

// EncryptZeroSeeded generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler; if `sampler` is nil, then the internal sampler will be used.
// montgomery: returns the result in the Montgomery domain.
func (enc Encryptor) encryptZeroSkFromC1QP(sk *rlwe.SecretKey, ct rlwe.Element[ringqp.Poly]) (err error) {

	levelQ, levelP := ct.LevelQ(), ct.LevelP()
	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	c1 := ct.Value[1]

	// ct = (0, g)
	enc.xeSampler.AtLevel(levelQ).Read(c1.Q)
	errCoeffs := make([]*big.Int, enc.params.N())
	errfloat := make([]float64, enc.params.N())
	for i := range errCoeffs {
		errCoeffs[i] = new(big.Int)
	}
	ringQP.RingQ.PolyToBigintCentered(c1.Q, 1, errCoeffs)
	for i := range errCoeffs {
		errfloat[i] = float64(errCoeffs[i].Int64())
	}

	skvar, err := stats.StandardDeviation(errfloat)
	if err != nil {
		panic(err)
	}
	//
	if 1 == 0 {
		fmt.Printf("NGSW err variance : %f\n", skvar)
	}

	for w := 0; w < len(c1.Q.Coeffs[0]); w++ {
		c1.Q.Coeffs[0][w] = ring.CRed(c1.Q.Coeffs[0][w], ringQP.RingQ.Modulus().Uint64())
	}
	// c1.Q.Coeffs[0]
	// fmt.Printf("error : %v\n", c1.Q.Coeffs[0])

	if levelP != -1 {
		fmt.Print("Should not be here\n")
		ringQ := ringQP.RingQ
		ringP := ringQP.RingP
		be := ring.NewBasisExtender(ringQ, ringP)

		be.ModUpQtoP(levelQ, levelP, c1.Q, c1.P)
	}

	ringQP.NTT(c1, c1)
	// ct[1] is assumed to be sampled in of the Montgomery domain,
	// thus -as will also be in the Montgomery domain (s is by default), therefore 'e'
	// must be switched to the Montgomery domain.
	ringQP.MForm(c1, c1)

	// (0, g * f_inv)
	ringQP.MulCoeffsMontgomery(c1, sk.Value, c1)

	if !ct.IsNTT {
		fmt.Printf("Not in NTT C1QP!\n")
		ringQP.INTT(c1, c1)
	}

	return
}

// montgomery check는 어떻게 하지.... 어우...
func (enc Encryptor) addPtDivfToCt(level int, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {
	// sk := enc.encKey.(*rlwe.SecretKey)
	ringQ := enc.params.RingQ().AtLevel(level)
	switch sk := enc.encKey.(type) {
	case *rlwe.SecretKey:
		var buff ring.Poly
		if pt.IsNTT {
			if ct.IsNTT {
				buff = pt.Value
			} else {
				buff = enc.buffQ[0]
				ringQ.NTT(pt.Value, buff)
			}
		} else {
			if ct.IsNTT {
				buff = enc.buffQ[0]
				ringQ.INTT(pt.Value, buff)
			} else {
				buff = pt.Value
			}
		}

		ringQ.MulCoeffsMontgomery(buff, sk.Value.Q, buff)
		ringQ.Add(ct.Value[1], buff, ct.Value[1])

	default:
		return
	}
}

func (enc Encryptor) addPtToCt(level int, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {
	ringQ := enc.params.RingQ().AtLevel(level)
	var buff ring.Poly
	if pt.IsNTT {
		if ct.IsNTT {
			buff = pt.Value
		} else {
			buff = enc.buffQ[0]
			ringQ.NTT(pt.Value, buff)
		}
	} else {
		if ct.IsNTT {
			buff = enc.buffQ[0]
			ringQ.INTT(pt.Value, buff)
		} else {
			buff = pt.Value
		}
	}

	ringQ.Add(ct.Value[1], buff, ct.Value[1])
}

// checkPk checks that a given pk is correct for the parameters.
func (enc Encryptor) checkSk(sk *rlwe.SecretKey) (err error) {
	if sk.Value.Q.N() != enc.params.N() {
		return fmt.Errorf("sk ring degree does not match params ring degree")
	}
	return
}

func (enc Encryptor) ShallowCopy() *Encryptor {
	return NewEncryptor(enc.params, enc.encKey)
}

func (enc Encryptor) WithKey(key rlwe.EncryptionKey) *Encryptor {
	switch key := key.(type) {
	case *rlwe.SecretKey:
		if err := enc.checkSk(key); err != nil {
			// Sanity check, this error should not happen.
			panic(fmt.Errorf("cannot WithKey: %w", err))
		}
	case nil:
		return &enc
	default:
		// Sanity check, this error should not happen.
		panic(fmt.Errorf("invalid key type, want *rlwe.SecretKey, or nil but have %T", key))
	}
	enc.encKey = key
	return &enc
}

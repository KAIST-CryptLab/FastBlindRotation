package lwe

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
)

func TestEncryption(t *testing.T) {
	// parameters of the TFHE sample
	var p uint64 = 2    // plaintext Modulus
	var q uint64 = 1024 // ciphertext Modulus = 2^11
	var n uint64 = 512  // TFHE dimension
	var m uint64 = 1    // message to encrypt/decrypt

	var paramsLWE = NewParameters(n, q, p) // n, q, p
	var err error
	var params hefloat.Parameters

	if params, err = hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN: 8,
		LogQ: []int{50},
	}); err != nil {
		panic(err)
	}

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	skCopy := sk.CopyNew()
	params.RingQ().AtLevel(0).INTT(skCopy.Value.Q, skCopy.Value.Q)
	params.RingQ().AtLevel(0).IMForm(skCopy.Value.Q, skCopy.Value.Q)

	enc := NewEncryptor(paramsLWE, params.RingQ(), skCopy.Value.Q)
	dec := NewDecryptor(paramsLWE, params.RingQ(), skCopy.Value.Q)

	ct := enc.EncryptNew(m)
	m_dec := dec.DecryptNew(ct)

	fmt.Printf("%v\n", m_dec)
}

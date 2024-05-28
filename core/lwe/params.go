package lwe

type Parameters struct {
	n uint64 // for ciphertext dimension
	p uint64 // plaintext modulus
	q uint64 // ciphertext modulus
}

func NewParameters(n, q, p uint64) (params Parameters) {
	params = Parameters{
		n: n,
		q: q,
		p: p,
	}
	return params
}

func (param Parameters) N() uint64 {
	return param.n
}

func (param Parameters) P() uint64 {
	return param.p
}

func (param Parameters) Q() uint64 {
	return param.q
}

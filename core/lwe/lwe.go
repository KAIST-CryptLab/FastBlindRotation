package lwe

type Ciphertext struct {
	A []uint64
	B uint64 // b - a * s = Δm + e
}

func (ct *Ciphertext) GetA() []uint64 {
	return ct.A
}

func (ct *Ciphertext) GetB() uint64 {
	return ct.B
}

func NewCiphertext(n int) (ct *Ciphertext) {
	return &Ciphertext{
		A: make([]uint64, n),
		B: 0,
	}
}

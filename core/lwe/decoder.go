package lwe

type Decoder struct {
	lweParam Parameters
}

func NewDecoder(lweParam Parameters) *Decoder {
	return &Decoder{lweParam: lweParam}
}

func (dec *Decoder) Decode(m uint64) uint64 {
	q := dec.lweParam.Q()
	p := dec.lweParam.P()

	// dec_m := (m % q) * 2 * p / q
	dec_m := ((m + (q / (4 * p))) % q) * 2 * p / q
	return dec_m
}

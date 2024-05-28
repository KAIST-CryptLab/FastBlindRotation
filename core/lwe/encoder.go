package lwe

type Encoder struct {
	lweParam Parameters
}

func NewEncoder(lweParam Parameters) *Encoder {
	return &Encoder{lweParam: lweParam}
}

func (enc *Encoder) Encode(m uint64, e int64) uint64 {
	q := enc.lweParam.Q()
	p := enc.lweParam.P()

	err := ((e % int64(q)) + int64(q)) % int64(q)
	enc_m := ((q/(2*p))*m + uint64(err) + q/(4*p)) % q // in case error is negative, q/4p adjusts the value
	return enc_m
}

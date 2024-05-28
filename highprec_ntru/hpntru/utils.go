package hpntru

import (
	"github.com/montanaflynn/stats"
	"github.com/tuneinsight/lattigo/v5/ring"
)

func Decomp(n, B uint64, d int) []uint64 {
	var digits = make([]uint64, d)
	var num = n
	var idx = 0
	for num != 0 {
		digits[idx] = num % B
		num /= B
		idx++
	}
	for i, j := 0, len(digits)-1; i < j; i, j = i+1, j-1 {
		digits[i], digits[j] = digits[j], digits[i]
	}
	return digits
}

func GetMeanStdev(data []float64) (mean float64, stdev float64) {
	var err error
	if mean, err = stats.Mean(data); err != nil {
		panic(err)
	}
	if stdev, err = stats.StandardDeviation(data); err != nil {
		panic(err)
	}
	return
}

/*
// MulBySmallMonomialMod2N multiplies pol by x^n, with 0 <= n < N
func mulBySmallMonomialMod2N(mask uint64, pol ring.Poly, n int) {
	if n != 0 {
		N := len(pol.Coeffs[0])
		pol.Coeffs[0] = append(pol.Coeffs[0][N-n:], pol.Coeffs[0][:N-n]...)
		tmp := pol.Coeffs[0]
		for j := 0; j < n; j++ {
			tmp[j] = -tmp[j] & mask
		}
	}
}

func normalizeInv(x, a, b float64) (y float64) {
	return (x*(b-a) + b + a) / 2.0
}

func scaleUp(value float64, scale float64, Q uint64) (res uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int

	isNegative = false
	if value < 0 {
		isNegative = true
		xFlo = big.NewFloat(-scale * value)
	} else {
		xFlo = big.NewFloat(scale * value)
	}

	xFlo.Add(xFlo, big.NewFloat(0.5))

	xInt = new(big.Int)
	xFlo.Int(xInt)
	xInt.Mod(xInt, bignum.NewInt(Q))

	res = xInt.Uint64()

	if isNegative {
		res = Q - res
	}

	return
}
*/

func MulByMonomial(ringQ *ring.Ring, F ring.Poly, k int) (res ring.Poly) {
	Q := ringQ.ModuliChain()[:ringQ.Level()+1]
	N := ringQ.N()
	res = ringQ.NewPoly()

	ringQ.MultByMonomial(F, k, res)
	for j, qi := range Q {
		for i := 0; i < N; i++ {
			res.Coeffs[j][i] %= qi
		}
	}

	return
}

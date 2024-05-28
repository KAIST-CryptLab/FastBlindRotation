package ntru

import (
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
)

func GetDegree(N int, F ringqp.Poly) (deg int) {
	deg = -1

	for i := N - 1; i >= 0; i-- {
		if F.Q.Coeffs[0][i] != 0 {
			deg = i
			break
		}
	}

	return
}

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

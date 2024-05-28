package ntru

import (
	"fmt"
	"math/big"

	"github.com/montanaflynn/stats"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/ring/ringqp"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

type KeyGenerator struct {
	*Encryptor
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as EvaluationKeys.
func NewKeyGenerator(params rlwe.ParameterProvider) *KeyGenerator {
	return &KeyGenerator{
		Encryptor: NewEncryptor(params, nil),
	}
}

// generating secret keys and its inverse
func (kgen KeyGenerator) GenSecretKeyPairNew() (sk *rlwe.SecretKey, sk_inv *rlwe.SecretKey) {
	ringQP := kgen.params.RingQP()
	sk_pol, sk_inv_pol := kgen.GenInvertPolyPair(ringQP)

	sk = &rlwe.SecretKey{Value: sk_pol}
	sk_inv = &rlwe.SecretKey{Value: sk_inv_pol}
	return
}

// now ringP is not being used
func (kgen KeyGenerator) GenInvertPolyPair(ringQP *ringqp.Ring) (sk ringqp.Poly, sk_inv ringqp.Poly) {
	ringQ := ringQP.RingQ
	levelQ := ringQ.Level()

	Sampler := kgen.xsSampler

	Q := ringQ.ModuliChain()[:levelQ+1]
	N := ringQ.N()

	sk = ringQP.NewPoly()
	sk_inv = ringQP.NewPoly()
	sk.Q = Sampler.ReadNew()

	// Checks sk variance
	skCoeffs := make([]*big.Int, ringQ.N())
	var sk_ss float64 = 0
	skfloat := make([]float64, ringQ.N())
	for i := range skCoeffs {
		skCoeffs[i] = new(big.Int)
	}
	ringQ.PolyToBigintCentered(sk.Q, 1, skCoeffs)
	for i := range skCoeffs {
		skfloat[i] = float64(skCoeffs[i].Int64())
		sk_ss += skfloat[i] * skfloat[i]
	}

	skvar, err := stats.StandardDeviation(skfloat)
	if err != nil {
		panic(err)
	}
	fmt.Printf("sk variance : %f\n", skvar)
	fmt.Printf("sk ss : %f\n", sk_ss)

	// setup for algorithm
	f := ringQP.NewPoly()
	b := ringQP.NewPoly()
	f.Copy(sk)

	for j := range Q {
		b.Q.Coeffs[j][0] = 1
	}
	c := ringQP.NewPoly() // initialized as 0
	k := 0

	for f.Q.Coeffs[0][0] == 0 {
		f.Q = MulByMonomial(ringQ, f.Q, -1)
		k++
	}

	tmp := ringQP.NewPoly()
	g := ringQP.NewPoly()

	tmp.Copy(f)
	for j := range Q {
		tmp.Q.Coeffs[j][0] = 0
	}
	tmp.Q = MulByMonomial(ringQ, tmp.Q, -1)

	if f.Q.Coeffs[0][0] == 1 {
		ringQP.Sub(g, tmp, g)
		for j := range Q {
			g.Q.Coeffs[j][N-1] = 1
		}
		f, g = g, f
		b, c = c, b
		ringQP.Sub(b, c, b)
	} else {
		ringQP.Add(g, tmp, g)
		for j := range Q {
			g.Q.Coeffs[j][N-1] = 1
		}
		f, g = g, f
		b, c = c, b
		ringQP.Add(b, c, b)
	}
	c.Q = MulByMonomial(ringQ, c.Q, 1)
	k++

	// now we have f, g, b, c for loop
	for {
		for f.Q.Coeffs[0][0] == 0 {
			f.Q = MulByMonomial(ringQ, f.Q, -1)
			c.Q = MulByMonomial(ringQ, c.Q, 1)
			k++
		}

		if GetDegree(N, f) == 0 {
			for i, s := range ringQ.SubRings[:levelQ+1] {
				f0_inv := new(big.Int).SetUint64(f.Q.Coeffs[i][0])
				f0_inv.ModInverse(f0_inv, bignum.NewInt(s.Modulus))
				f0_inv.Sub(bignum.NewInt(s.Modulus), f0_inv)

				s.MulScalarMontgomery(b.Q.Coeffs[i], ring.MForm(f0_inv.Uint64(), s.Modulus, s.BRedConstant), b.Q.Coeffs[i])
			}

			b.Q = MulByMonomial(ringQ, b.Q, N-k)
			break
		}

		if GetDegree(N, f) < GetDegree(N, g) {
			f, g = g, f
			b, c = c, b
		}

		// multiply g0_inv * f0 to g

		tmp = ringQP.NewPoly()
		tmp2 := ringQP.NewPoly()

		for i, s := range ringQ.SubRings[:levelQ+1] {
			g0_inv := new(big.Int).SetUint64(g.Q.Coeffs[i][0])
			g0_inv.ModInverse(g0_inv, bignum.NewInt(s.Modulus))

			u := new(big.Int).SetUint64(f.Q.Coeffs[i][0])
			u.Mul(u, g0_inv)
			u.Mod(u, bignum.NewInt(s.Modulus))

			s.MulScalarMontgomery(g.Q.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp.Q.Coeffs[i])
			s.MulScalarMontgomery(c.Q.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp2.Q.Coeffs[i])
		}

		ringQP.Sub(f, tmp, f)
		ringQP.Sub(b, tmp2, b)
	}

	ringQP.NTT(sk, sk)
	ringQP.MForm(sk, sk)

	ringQP.NTT(b, b)
	ringQP.MForm(b, sk_inv)

	return
}

func GenInvertPolyPairQP(ringQP *ringqp.Ring) (sk ringqp.Poly, sk_inv ringqp.Poly) {
	ringQ := ringQP.RingQ
	ringP := ringQP.RingP
	be := ring.NewBasisExtender(ringQ, ringP)

	levelQ := ringQ.Level()
	levelP := ringP.Level()

	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}

	ternarySamplerMontgomeryQ, err := ring.NewSampler(prng, ringQ, ring.Ternary{H: 3}, false)
	if err != nil {
		panic(err)
	}

	Q := ringQ.ModuliChain()[:levelQ+1]
	P := ringP.ModuliChain()[:levelP+1]
	N := ringQ.N()

	sk = ringQP.NewPoly()
	sk_inv = ringQP.NewPoly()
	sk.Q = ternarySamplerMontgomeryQ.ReadNew()
	// ringQP.ExtendBasisSmallNormAndCenter(sk.Q, levelP, sk.Q, sk.P)
	be.ModUpQtoP(levelQ, levelP, sk.Q, sk.P)

	// setup for algorithm
	f := ringQP.NewPoly()
	b := ringQP.NewPoly()
	f.Copy(sk)

	for j := range Q {
		b.Q.Coeffs[j][0] = 1
	}
	for j := range P {
		b.P.Coeffs[j][0] = 1
	}
	c := ringQP.NewPoly() // initialized as 0
	k := 0

	for f.Q.Coeffs[0][0] == 0 {
		f.Q = MulByMonomial(ringQ, f.Q, -1)
		f.P = MulByMonomial(ringP, f.P, -1)
		k++
	}

	tmp := ringQP.NewPoly()
	g := ringQP.NewPoly()

	tmp.Copy(f)
	for j := range Q {
		tmp.Q.Coeffs[j][0] = 0
	}
	for j := range P {
		tmp.P.Coeffs[j][0] = 0
	}
	tmp.Q = MulByMonomial(ringQ, tmp.Q, -1)
	tmp.P = MulByMonomial(ringP, tmp.P, -1)

	if f.Q.Coeffs[0][0] == 1 {
		ringQP.Sub(g, tmp, g)
		for j := range Q {
			g.Q.Coeffs[j][N-1] = 1
		}
		for j := range P {
			g.P.Coeffs[j][N-1] = 1
		}
		f, g = g, f
		b, c = c, b
		ringQP.Sub(b, c, b)
	} else {
		ringQP.Add(g, tmp, g)
		for j := range Q {
			g.Q.Coeffs[j][N-1] = 1
		}
		for j := range P {
			g.P.Coeffs[j][N-1] = 1
		}
		f, g = g, f
		b, c = c, b
		ringQP.Add(b, c, b)
	}
	c.Q = MulByMonomial(ringQ, c.Q, 1)
	c.P = MulByMonomial(ringP, c.P, 1)
	k++

	// now we have f, g, b, c for loop
	for {
		for f.Q.Coeffs[0][0] == 0 {
			f.Q = MulByMonomial(ringQ, f.Q, -1)
			f.P = MulByMonomial(ringP, f.P, -1)
			c.Q = MulByMonomial(ringQ, c.Q, 1)
			c.P = MulByMonomial(ringP, c.P, 1)
			k++
		}

		if GetDegree(N, f) == 0 {
			for i, s := range ringQ.SubRings[:levelQ+1] {
				f0_inv := new(big.Int).SetUint64(f.Q.Coeffs[i][0])
				f0_inv.ModInverse(f0_inv, bignum.NewInt(s.Modulus))
				f0_inv.Sub(bignum.NewInt(s.Modulus), f0_inv)

				s.MulScalarMontgomery(b.Q.Coeffs[i], ring.MForm(f0_inv.Uint64(), s.Modulus, s.BRedConstant), b.Q.Coeffs[i])
			}
			for i, s := range ringP.SubRings[:levelP+1] {
				f0_inv := new(big.Int).SetUint64(f.P.Coeffs[i][0])
				f0_inv.ModInverse(f0_inv, bignum.NewInt(s.Modulus))
				f0_inv.Sub(bignum.NewInt(s.Modulus), f0_inv)

				s.MulScalarMontgomery(b.P.Coeffs[i], ring.MForm(f0_inv.Uint64(), s.Modulus, s.BRedConstant), b.P.Coeffs[i])
			}

			b.Q = MulByMonomial(ringQ, b.Q, N-k)
			b.P = MulByMonomial(ringP, b.P, N-k)
			break
		}

		if GetDegree(N, f) < GetDegree(N, g) {
			f, g = g, f
			b, c = c, b
		}

		// multiply g0_inv * f0 to g

		tmp = ringQP.NewPoly()
		tmp2 := ringQP.NewPoly()

		for i, s := range ringQ.SubRings[:levelQ+1] {
			g0_inv := new(big.Int).SetUint64(g.Q.Coeffs[i][0])
			g0_inv.ModInverse(g0_inv, bignum.NewInt(s.Modulus))

			u := new(big.Int).SetUint64(f.Q.Coeffs[i][0])
			u.Mul(u, g0_inv)
			u.Mod(u, bignum.NewInt(s.Modulus))

			s.MulScalarMontgomery(g.Q.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp.Q.Coeffs[i])
			s.MulScalarMontgomery(c.Q.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp2.Q.Coeffs[i])
		}

		for i, s := range ringP.SubRings[:levelP+1] {
			g0_inv := new(big.Int).SetUint64(g.P.Coeffs[i][0])
			g0_inv.ModInverse(g0_inv, bignum.NewInt(s.Modulus))

			u := new(big.Int).SetUint64(f.P.Coeffs[i][0])
			u.Mul(u, g0_inv)
			u.Mod(u, bignum.NewInt(s.Modulus))

			s.MulScalarMontgomery(g.P.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp.P.Coeffs[i])
			s.MulScalarMontgomery(c.P.Coeffs[i], ring.MForm(u.Uint64(), s.Modulus, s.BRedConstant), tmp2.P.Coeffs[i])
		}

		ringQP.Sub(f, tmp, f)
		ringQP.Sub(b, tmp2, b)
	}

	ringQP.NTT(sk, sk)
	ringQP.MForm(sk, sk)

	ringQP.NTT(b, b)
	ringQP.MForm(b, sk_inv)

	return
}

package rlwe

import (
	"bufio"
	"fmt"
	"io"

	"github.com/google/go-cmp/cmp"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils/buffer"
	"github.com/tuneinsight/lattigo/v4/utils/structs"
)

// SecretKey is a type for generic RLWE secret keys.
// The Value field stores the polynomial in NTT and Montgomery form.
type SecretKey struct {
	Value ringqp.Poly
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params ParametersInterface) *SecretKey {
	return &SecretKey{Value: params.RingQP().NewPoly()}
}

func (sk SecretKey) Equal(other *SecretKey) bool {
	return cmp.Equal(sk.Value, other.Value)
}

// LevelQ returns the level of the modulus Q of the target.
func (sk SecretKey) LevelQ() int {
	return sk.Value.Q.Level()
}

// LevelP returns the level of the modulus P of the target.
// Returns -1 if P is absent.
func (sk SecretKey) LevelP() int {
	return sk.Value.P.Level()
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk SecretKey) CopyNew() *SecretKey {
	return &SecretKey{sk.Value.CopyNew()}
}

// BinarySize returns the serialized size of the object in bytes.
func (sk SecretKey) BinarySize() (dataLen int) {
	return sk.Value.BinarySize()
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (sk SecretKey) WriteTo(w io.Writer) (n int64, err error) {
	return sk.Value.WriteTo(w)
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (sk *SecretKey) ReadFrom(r io.Reader) (n int64, err error) {
	return sk.Value.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (sk SecretKey) MarshalBinary() (p []byte, err error) {
	return sk.Value.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (sk *SecretKey) UnmarshalBinary(p []byte) (err error) {
	return sk.Value.UnmarshalBinary(p)
}

type tupleQP [2]ringqp.Poly

// NewPublicKey returns a new PublicKey with zero values.
func newTupleQP(params ParametersInterface) (pk tupleQP) {
	return [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}
}

// NewPublicKey returns a new PublicKey with zero values.
func newTupleQPAtLevel(params ParametersInterface, levelQ, levelP int) (pk tupleQP) {
	rqp := params.RingQP().AtLevel(levelQ, levelP)
	return [2]ringqp.Poly{rqp.NewPoly(), rqp.NewPoly()}
}

// CopyNew creates a deep copy of the target PublicKey and returns it.
func (p tupleQP) CopyNew() tupleQP {
	return [2]ringqp.Poly{p[0].CopyNew(), p[1].CopyNew()}
}

// Equal performs a deep equal.
func (p tupleQP) Equal(other *tupleQP) bool {
	return p[0].Equal(&other[0]) && p[1].Equal(&other[1])
}

func (p tupleQP) BinarySize() int {
	return structs.Vector[ringqp.Poly](p[:]).BinarySize()
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (p tupleQP) WriteTo(w io.Writer) (n int64, err error) {
	v := structs.Vector[ringqp.Poly](p[:])
	return v.WriteTo(w)
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (p *tupleQP) ReadFrom(r io.Reader) (n int64, err error) {
	v := structs.Vector[ringqp.Poly](p[:])
	n, err = v.ReadFrom(r)
	if len(v) != 2 {
		return n, fmt.Errorf("bad public key format")
	}
	return
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p tupleQP) MarshalBinary() ([]byte, error) {
	v := structs.Vector[ringqp.Poly](p[:])
	return v.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (p *tupleQP) UnmarshalBinary(b []byte) error {
	v := structs.Vector[ringqp.Poly](p[:])
	err := v.UnmarshalBinary(b)
	if len(v) != 2 {
		return fmt.Errorf("bad public key format")
	}
	return err
}

// PublicKey is a type for generic RLWE public keys.
// The Value field stores the polynomials in NTT and Montgomery form.
type PublicKey struct {
	Value tupleQP
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params ParametersInterface) (pk *PublicKey) {
	return &PublicKey{Value: newTupleQP(params)}
}

// CopyNew creates a deep copy of the target PublicKey and returns it.
func (p PublicKey) CopyNew() *PublicKey {
	return &PublicKey{Value: p.Value.CopyNew()}
}

// Equal performs a deep equal.
func (p PublicKey) Equal(other *PublicKey) bool {
	return p.Value.Equal(&other.Value)
}

func (p PublicKey) BinarySize() int {
	return p.Value.BinarySize()
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (p PublicKey) WriteTo(w io.Writer) (n int64, err error) {
	return p.Value.WriteTo(w)
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (p *PublicKey) ReadFrom(r io.Reader) (n int64, err error) {
	return p.Value.ReadFrom(r)
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (p PublicKey) MarshalBinary() ([]byte, error) {
	return p.Value.MarshalBinary()
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (p *PublicKey) UnmarshalBinary(b []byte) error {
	return p.Value.UnmarshalBinary(b)
}

// EvaluationKey is a public key indended to be used during the evaluation phase of a homomorphic circuit.
// It provides a one way public and non-interactive re-encryption from a ciphertext encrypted under `skIn`
// to a ciphertext encrypted under `skOut`.
//
// Such re-encryption is for example used for:
//
// - Homomorphic relinearization: re-encryption of a quadratic ciphertext (that requires (1, sk sk^2) to be decrypted)
// to a linear ciphertext (that required (1, sk) to be decrypted). In this case skIn = sk^2 an skOut = sk.
//
// - Homomorphic automorphisms: an automorphism in the ring Z[X]/(X^{N}+1) is defined as pi_k: X^{i} -> X^{i^k} with
// k coprime to 2N. Pi_sk is for exampled used during homomorphic slot rotations. Applying pi_k to a ciphertext encrypted
// under sk generates a new ciphertext encrypted under pi_k(sk), and an Evaluationkey skIn = pi_k(sk) to skOut = sk
// is used to bring it back to its original key.
type EvaluationKey struct {
	GadgetCiphertext
}

// NewEvaluationKey returns a new EvaluationKey with pre-allocated zero-value
func NewEvaluationKey(params ParametersInterface, levelQ, levelP int) *EvaluationKey {
	//evk := new(EvaluationKey)
	// drns := params.DecompRNS(levelQ, levelP)
	// dpw2 := params.DecompPw2(levelQ, levelP)
	// evk.Value = make(structs.Matrix[tupleQP], drns)
	// for i := range evk.Value {
	// 	evk.Value[i] = make([][2]ringqp.Poly, dpw2)
	// 	for j := range evk.Value[i] {
	// 		evk.Value[i][j] = NewPublicKey(params).Value
	// 	}
	// }
	return &EvaluationKey{GadgetCiphertext: *NewGadgetCiphertext(params, levelQ, levelP, params.DecompRNS(levelQ, levelP), params.DecompPw2(levelQ, levelP))}
}

// CopyNew creates a deep copy of the target EvaluationKey and returns it.
func (evk EvaluationKey) CopyNew() *EvaluationKey {
	return &EvaluationKey{GadgetCiphertext: *evk.GadgetCiphertext.CopyNew()}
}

// Equal performs a deep equal.
func (evk EvaluationKey) Equal(other *EvaluationKey) bool {
	return evk.GadgetCiphertext.Equal(&other.GadgetCiphertext)
}

// RelinearizationKey is type of evaluation key used for ciphertext multiplication compactness.
// The Relinearization key encrypts s^{2} under s and is used to homomorphically re-encrypt the
// degree 2 term of a ciphertext (the term that decrypt with s^{2}) into a degree 1 term
// (a term that decrypts with s).
type RelinearizationKey struct {
	EvaluationKey
}

// NewRelinearizationKey allocates a new RelinearizationKey with zero coefficients.
func NewRelinearizationKey(params ParametersInterface) *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP())}
}

// CopyNew creates a deep copy of the object and returns it.
func (rlk RelinearizationKey) CopyNew() *RelinearizationKey {
	return &RelinearizationKey{EvaluationKey: *rlk.EvaluationKey.CopyNew()}
}

// Equal performs a deep equal.
func (rlk RelinearizationKey) Equal(other *RelinearizationKey) bool {
	return rlk.EvaluationKey.Equal(&other.EvaluationKey)
}

// GaloisKey is a type of evaluation key used to evaluate automorphisms on ciphertext.
// An automorphism pi: X^{i} -> X^{i*GaloisElement} changes the key under which the
// ciphertext is encrypted from s to pi(s). Thus, the ciphertext must be re-encrypted
// from pi(s) to s to ensure correctness, which is done with the corresponding GaloisKey.
//
// Lattigo implements automorphismes differently than the usual way (which is to first
// apply the automorphism and then the evaluation key). Instead the order of operations
// is reversed, the GaloisKey for pi^{-1} is evaluated on the ciphertext, outputing a
// ciphertext encrypted under pi^{-1}(s), and then the automorphism pi is applied. This
// enables a more efficient evaluation, by only having to apply the automorphism on the
// final result (instead of having to apply it on the decomposed ciphertext).
type GaloisKey struct {
	GaloisElement uint64
	NthRoot       uint64
	EvaluationKey
}

// NewGaloisKey allocates a new GaloisKey with zero coefficients and GaloisElement set to zero.
func NewGaloisKey(params ParametersInterface) *GaloisKey {
	return &GaloisKey{EvaluationKey: *NewEvaluationKey(params, params.MaxLevelQ(), params.MaxLevelP()), NthRoot: params.RingQ().NthRoot()}
}

// Equal returns true if the two objects are equal.
func (gk GaloisKey) Equal(other *GaloisKey) bool {
	return gk.GaloisElement == other.GaloisElement && gk.NthRoot == other.NthRoot && cmp.Equal(gk.EvaluationKey, other.EvaluationKey)
}

// CopyNew creates a deep copy of the object and returns it
func (gk GaloisKey) CopyNew() *GaloisKey {
	return &GaloisKey{
		GaloisElement: gk.GaloisElement,
		NthRoot:       gk.NthRoot,
		EvaluationKey: *gk.EvaluationKey.CopyNew(),
	}
}

// BinarySize returns the serialized size of the object in bytes.
func (gk GaloisKey) BinarySize() (size int) {
	return gk.EvaluationKey.BinarySize() + 16
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (gk GaloisKey) WriteTo(w io.Writer) (n int64, err error) {
	switch w := w.(type) {
	case buffer.Writer:

		var inc int

		if inc, err = buffer.WriteUint64(w, gk.GaloisElement); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		if inc, err = buffer.WriteUint64(w, gk.NthRoot); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		var inc2 int64
		if inc2, err = gk.EvaluationKey.WriteTo(w); err != nil {
			return n + inc2, err
		}

		n += inc2

		return

	default:
		return gk.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (gk *GaloisKey) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:

		var inc int

		if inc, err = buffer.ReadUint64(r, &gk.GaloisElement); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		if inc, err = buffer.ReadUint64(r, &gk.NthRoot); err != nil {
			return n + int64(inc), err
		}

		n += int64(inc)

		var inc2 int64
		if inc2, err = gk.EvaluationKey.ReadFrom(r); err != nil {
			return n + inc2, err
		}

		n += inc2

		return
	default:
		return gk.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (gk GaloisKey) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(gk.BinarySize())
	_, err = gk.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (gk *GaloisKey) UnmarshalBinary(p []byte) (err error) {
	_, err = gk.ReadFrom(buffer.NewBuffer(p))
	return
}

// EvaluationKeySet is an interface implementing methods
// to load the RelinearizationKey and GaloisKeys in the Evaluator.
// Implementations of this interface must be safe for concurrent use.
type EvaluationKeySet interface {

	// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
	GetGaloisKey(galEl uint64) (evk *GaloisKey, err error)

	// GetGaloisKeysList returns the list of all the Galois elements
	// for which a Galois key exists in the object.
	GetGaloisKeysList() (galEls []uint64)

	// GetRelinearizationKey retrieves the RelinearizationKey.
	GetRelinearizationKey() (evk *RelinearizationKey, err error)
}

// MemEvaluationKeySet is a basic in-memory implementation of the EvaluationKeySet interface.
type MemEvaluationKeySet struct {
	Rlk *RelinearizationKey
	Gks structs.Map[uint64, GaloisKey]
}

// NewMemEvaluationKeySet returns a new EvaluationKeySet with the provided RelinearizationKey and GaloisKeys.
func NewMemEvaluationKeySet(relinKey *RelinearizationKey, galoisKeys ...*GaloisKey) (eks *MemEvaluationKeySet) {
	eks = &MemEvaluationKeySet{Gks: map[uint64]*GaloisKey{}}
	eks.Rlk = relinKey
	for _, k := range galoisKeys {
		eks.Gks[k.GaloisElement] = k
	}
	return eks
}

// GetGaloisKey retrieves the Galois key for the automorphism X^{i} -> X^{i*galEl}.
func (evk MemEvaluationKeySet) GetGaloisKey(galEl uint64) (gk *GaloisKey, err error) {
	var ok bool
	if gk, ok = evk.Gks[galEl]; !ok {
		return nil, fmt.Errorf("GaloiKey[%d] is nil", galEl)
	}

	return
}

// GetGaloisKeysList returns the list of all the Galois elements
// for which a Galois key exists in the object.
func (evk MemEvaluationKeySet) GetGaloisKeysList() (galEls []uint64) {

	if evk.Gks == nil {
		return []uint64{}
	}

	galEls = make([]uint64, len(evk.Gks))

	var i int
	for galEl := range evk.Gks {
		galEls[i] = galEl
		i++
	}

	return
}

// GetRelinearizationKey retrieves the RelinearizationKey.
func (evk MemEvaluationKeySet) GetRelinearizationKey() (rk *RelinearizationKey, err error) {
	if evk.Rlk != nil {
		return evk.Rlk, nil
	}

	return nil, fmt.Errorf("RelinearizationKey is nil")
}

func (evk MemEvaluationKeySet) BinarySize() (size int) {

	size++
	if evk.Rlk != nil {
		size += evk.Rlk.BinarySize()
	}

	size++
	if evk.Gks != nil {
		size += evk.Gks.BinarySize()
	}

	return
}

// WriteTo writes the object on an io.Writer. It implements the io.WriterTo
// interface, and will write exactly object.BinarySize() bytes on w.
//
// Unless w implements the buffer.Writer interface (see lattigo/utils/buffer/writer.go),
// it will be wrapped into a bufio.Writer. Since this requires allocations, it
// is preferable to pass a buffer.Writer directly:
//
//   - When writing multiple times to a io.Writer, it is preferable to first wrap the
//     io.Writer in a pre-allocated bufio.Writer.
//   - When writing to a pre-allocated var b []byte, it is preferable to pass
//     buffer.NewBuffer(b) as w (see lattigo/utils/buffer/buffer.go).
func (evk MemEvaluationKeySet) WriteTo(w io.Writer) (int64, error) {
	switch w := w.(type) {
	case buffer.Writer:

		var inc int
		var n, inc64 int64
		var err error

		if evk.Rlk != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.Rlk.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		if evk.Gks != nil {
			if inc, err = buffer.WriteUint8(w, 1); err != nil {
				return int64(inc), err
			}

			n += int64(inc)

			if inc64, err = evk.Gks.WriteTo(w); err != nil {
				return n + inc64, err
			}

			n += inc64

		} else {
			if inc, err = buffer.WriteUint8(w, 0); err != nil {
				return int64(inc), err
			}
			n += int64(inc)
		}

		return n, w.Flush()

	default:
		return evk.WriteTo(bufio.NewWriter(w))
	}
}

// ReadFrom reads on the object from an io.Writer. It implements the
// io.ReaderFrom interface.
//
// Unless r implements the buffer.Reader interface (see see lattigo/utils/buffer/reader.go),
// it will be wrapped into a bufio.Reader. Since this requires allocation, it
// is preferable to pass a buffer.Reader directly:
//
//   - When reading multiple values from a io.Reader, it is preferable to first
//     first wrap io.Reader in a pre-allocated bufio.Reader.
//   - When reading from a var b []byte, it is preferable to pass a buffer.NewBuffer(b)
//     as w (see lattigo/utils/buffer/buffer.go).
func (evk *MemEvaluationKeySet) ReadFrom(r io.Reader) (n int64, err error) {
	switch r := r.(type) {
	case buffer.Reader:
		var inc int
		var n, inc64 int64
		var err error

		var hasKey uint8

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		if hasKey == 1 {

			if evk.Rlk == nil {
				evk.Rlk = new(RelinearizationKey)
			}

			if inc64, err = evk.Rlk.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		if inc, err = buffer.ReadUint8(r, &hasKey); err != nil {
			return int64(inc), err
		}

		n += int64(inc)

		if hasKey == 1 {

			if evk.Gks == nil {
				evk.Gks = structs.Map[uint64, GaloisKey]{}
			}

			if inc64, err = evk.Gks.ReadFrom(r); err != nil {
				return n + inc64, err
			}

			n += inc64
		}

		return n, nil

	default:
		return evk.ReadFrom(bufio.NewReader(r))
	}
}

// MarshalBinary encodes the object into a binary form on a newly allocated slice of bytes.
func (evk MemEvaluationKeySet) MarshalBinary() (p []byte, err error) {
	buf := buffer.NewBufferSize(evk.BinarySize())
	_, err = evk.WriteTo(buf)
	return buf.Bytes(), err
}

// UnmarshalBinary decodes a slice of bytes generated by
// MarshalBinary or WriteTo on the object.
func (evk *MemEvaluationKeySet) UnmarshalBinary(p []byte) (err error) {
	_, err = evk.ReadFrom(buffer.NewBuffer(p))
	return
}

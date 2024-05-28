package hpntru

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
)

var (
	Q10                   uint64 = 995329
	Q11                   uint64 = 44421121 // 536903681
	testNTRUParamsLiteral        = []rlwe.ParametersLiteral{
		{
			LogN:    10,            // should be same as q for LWE
			Q:       []uint64{Q10}, // ~2^19.9
			NTTFlag: true,
			Xe:      ring.Ternary{H: 50}, // ring.Ternary{P: 1. / 2.},
			Xs:      ring.Ternary{H: 512},
		},
		{
			LogN:    10,            // should be same as q for LWE
			Q:       []uint64{Q10}, // ~2^19.9
			NTTFlag: true,
			Xe:      ring.Ternary{H: 50}, // Xe:      ring.Ternary{H: 20},
			Xs:      ring.Ternary{H: 512},
		},
		{
			LogN:    11,            // should be same as q for LWE
			Q:       []uint64{Q11}, // ~2^29
			NTTFlag: true,
			Xe:      ring.Ternary{P: 1. / 2.},
			Xs:      ring.Ternary{H: 1024},
		},
	}

	testRLWEParamsLiteral = []rlwe.ParametersLiteral{
		{
			LogN:    10,                  // should be same as q for LWE
			Q:       []uint64{0x7fff801}, //
			NTTFlag: true,
			Xe:      ring.Ternary{P: 1. / 2.},
		},
		{
			LogN:    10,                  // should be same as q for LWE
			Q:       []uint64{0x7fff801}, //
			NTTFlag: true,
			Xe:      ring.Ternary{H: 5},
		},
		{
			LogN:    11,                   // should be same as q for LWE
			Q:       []uint64{1073750017}, //
			NTTFlag: true,
			Xe:      ring.Ternary{P: 1. / 2.},
			Xs:      ring.Ternary{H: 1024},
		},
	}

	testMVBRParamsLiteral = []rlwe.ParametersLiteral{
		{
			LogN:    10,            // should be same as q for LWE
			Q:       []uint64{Q10}, // ~2^19.9
			NTTFlag: true,
			Xe:      ring.Ternary{H: 5},
		},
		{
			LogN:    10,            // should be same as q for LWE
			Q:       []uint64{Q10}, // ~2^19.9
			NTTFlag: true,
			Xe:      ring.Ternary{H: 1},
		},
		{
			LogN:    11,            // should be same as q for LWE
			Q:       []uint64{Q11}, // ~2^29
			NTTFlag: true,
			Xe:      ring.Ternary{H: 32},
			Xs:      ring.Ternary{H: 1024},
		},
	}

	// testPlaintextModulus = []uint64{0x101, 0xffc001}
	// testParams = []heint.ParametersLiteral{testInsecure}
)

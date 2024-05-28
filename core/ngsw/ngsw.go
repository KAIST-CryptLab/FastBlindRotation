package ngsw

import (
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
)

var (
	testParamsLiteral = []hefloat.ParametersLiteral{
		{
			LogN:            4,
			LogQ:            []int{50, 40, 40},
			LogP:            []int{60},
			LogDefaultScale: 40,
			RingType:        ring.Standard,
		},
		{
			LogN:            14,
			LogQ:            []int{50, 40, 40, 40, 40, 40, 40, 40},
			LogP:            []int{60},
			LogDefaultScale: 40,
			RingType:        ring.Standard,
		},
	}
)

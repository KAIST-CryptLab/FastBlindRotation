package ngsw

import (
	"container/list"
	"math"
	"math/bits"
)

type Queue struct {
	v *list.List
}

func NewQueue() *Queue {
	return &Queue{list.New()}
}

func (q *Queue) Push(v interface{}) {
	q.v.PushBack(v)
}

func (q *Queue) Pop() interface{} {
	front := q.v.Front()
	if front == nil {
		return nil
	}
	return q.v.Remove(front)
}

func (q *Queue) Len() int {
	return q.v.Len()
}

// used for generating bit-reversal permutation
func BitRevGen(n uint) []uint {
	m := uint(math.Logb(float64(n)))
	r := make([]uint, n)
	for i := uint(0); i < n; i++ {
		r[i] = bits.Reverse(i) >> (bits.UintSize - m)
	}
	return r
}

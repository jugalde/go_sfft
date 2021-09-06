package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Filter ...
type Filter struct {
	time  []complex128
	sizet int
	freq  []complex128
}

func shift(x []complex128, r int) []complex128 {
	return append(x[r+1:], x[:r+1]...)
}

// Chebyshev ..
func Chebyshev(m, x float64) float64 {
	if math.Abs(x) <= 1 {
		return math.Cos(m * math.Acos(x))
	}
	inner := complex(m, 0) * cmplx.Acosh(complex(x, 0))
	return real(cmplx.Cosh(inner))
}
func makeDolphChebyshevT(lobeFrac, tolerance float64) ([]complex128, int) {
	w := int((1 / math.Pi) * (1 / lobeFrac) * math.Acosh(1./tolerance))
	if w%2 == 0 {
		w--
	}
	x := make([]complex128, w)
	t0 := math.Cosh(math.Acosh(1/tolerance) / (float64)(w-1))

	for i := 0; i < w; i++ {
		x[i] = complex(Chebyshev((float64)(w-1), t0*math.Cos((float64)(math.Pi)*(float64)(i)/(float64)(w)))*tolerance, 0)
	}
	fftwDft(x, false)
	x = shift(x, w/2)

	for i := 0; i < w; i++ {
		x[i] = complex(real(x[i]), 0)
	}
	return x, w
}

func makeMultipleT(x []complex128, w, n, b int) (Filter, error) {
	if b > n || w > n {
		return Filter{}, fmt.Errorf("b and w must be less than n")
	}
	g := make([]complex128, n)
	h := make([]complex128, n)
	copy(g, x[w/2:])
	copy(g[n-w/2:], x[:w/2])

	fftwDft(g, false)
	s := complex(0, 0)

	for i := 0; i < b; i++ {
		s += g[i]
	}

	var max float64 = 0
	offset := b / 2
	for i := 0; i < n; i++ {
		h[(i+n+offset)%n] = s
		max = math.Max(max, cmplx.Abs(s))
		s = s + (g[(i+b)%n] - g[i])
	}

	for i := 0; i < n; i++ {
		h[i] /= complex(max, 0)
	}

	offsetC := complex(1, 0)
	step := cmplx.Exp(complex((float64)(-2*math.Pi)*(float64)(w/2)/(float64)(n), 0) * complex(0, 1))

	for i := 0; i < n; i++ {
		h[i] *= offsetC
		offsetC *= step
	}

	hCopy := make([]complex128, n)
	copy(hCopy, h)
	fftwDft(hCopy, true)
	copy(x, hCopy[:w])

	for i := -0; i < w; i++ {
		x[i] /= complex((float64)(n), 0)

	}

	return Filter{
		time:  x,
		sizet: w,
		freq:  h,
	}, nil
}

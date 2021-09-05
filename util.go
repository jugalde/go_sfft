package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"strconv"

	"github.com/cpmech/gosl/fun/fftw"
)

type listOfFloats [][]float64

// Len ...
func (s listOfFloats) Len() int {
	return len(s)
}

// Swap ...
func (s listOfFloats) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

// Less ...
func (s listOfFloats) Less(i, j int) bool {
	if s[i][0] == s[j][0] {
		return s[i][1] < s[j][1]
	}
	return s[i][0] < s[j][0]
}

type floats []float64

// Len ...
func (s floats) Len() int {
	return len(s)
}

// Swap ...
func (s floats) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

// Less ...
func (s floats) Less(i, j int) bool {
	return s[i] < s[j]
}

type listOfPairs [][]int

// Len ...
func (s listOfPairs) Len() int {
	return len(s)
}

// Swap ...
func (s listOfPairs) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

// Less ...
func (s listOfPairs) Less(i, j int) bool {
	if s[i][0] == s[j][0] {
		return s[i][1] < s[j][1]
	}
	return s[i][0] < s[j][0]
}

// n must be less than len a
func memsetComplex(a []complex128, v complex128, n int) {
	if len(a) == 0 {
		return
	}
	a[0] = v
	for bp := 1; bp < n; bp *= 2 {
		copy(a[bp:], a[:bp])
	}
}

// n must be less than len a
func memsetInt(a []int, v, n int) {
	if len(a) == 0 {
		return
	}
	a[0] = v
	for bp := 1; bp < n; bp *= 2 {
		copy(a[bp:], a[:bp])
	}
}

func timesmod(x, a, n int) int {
	return (x * a) % n
}

func radix(byt, size int, A, TEMP []int) {

	COUNT := make([]int, 256)
	byt = byt << 3
	for i := 0; i < size; i++ {
		COUNT[((A[i])>>(byt))&0xFF]++
	}

	for i := 0; i < 256; i++ {
		COUNT[i] += COUNT[i-1]
	}
	for i := size - 1; i >= 0; i-- {
		TEMP[COUNT[(A[i]>>(byt))&0xFF]-1] = A[i]
		COUNT[(A[i]>>(byt))&0xFF]--
	}
}

func radixSort(A []int, size int) {
	TEMP := make([]int, size)
	for i := 0; i < strconv.IntSize; i++ {
		radix(i, size, A, TEMP)
		radix(i+1, size, TEMP, A)
	}

}

// GCD via Euclidean algorithm
func GCD(a, b int) int {
	for b != 0 {
		t := b
		b = a % b
		a = t
	}
	return a
}

// crappy inversion code I stole from elsewhere
// Undefined if gcd(a, n) > 1
func modInverse(a, n int) int {
	i := n
	v := 0
	d := 1
	for a > 0 {
		t := i / a
		x := a
		a = i % x
		i = x
		x = d
		d = v - t*x
		v = x
	}
	v %= n
	if v < 0 {
		v = (v + n) % n
	}
	return v
}

func binomialCdf(prob float64, n, needed int) float64 {
	var ans float64 = 0
	var choose float64 = 1
	for i := n; i >= needed; i-- {
		ans += choose * math.Pow(prob, (float64)(i)) * math.Pow((float64)(1-prob), (float64)(n-i))
		choose = choose * (float64)(i) / (float64)(n-i+1)
	}
	return ans
}

func floorToPow2(x float64) int {
	var ans uint
	for ans = 1; (float64)(ans) <= x; ans <<= 1 {
		continue
	}
	return (int)(ans / 2)
}

// AWGN ...
func AWGN(x []complex128, n int, stdNoise float64) float64 {

	if stdNoise == 0 {
		return 1000000000.
	}

	gn := complex(0, 0)

	var sigPower float64 = 0
	var noisePower float64 = 0
	var snr float64
	var u, v float64

	for h := 0; h < n; h++ {
		sigPower += cmplx.Abs(x[h]) * cmplx.Abs(x[h])

		u = rand.Float64()
		v = rand.Float64()
		gn = complex(stdNoise*math.Sqrt((float64)(-2)*math.Log(u)), 0) * cmplx.Exp(complex(2*math.Pi*v, 0)*complex(0, 1))

		noisePower += -2 * math.Log(u)

		x[h] += gn
	}

	noisePower = noisePower * stdNoise * stdNoise
	snr = sigPower / noisePower

	return snr
}

func fftwDft(out []complex128, backwards bool) {
	plan := fftw.NewPlan1d(out, backwards, false)
	plan.Execute()
	plan.Free()
}

func nthElementImmutable(input []float64, n, num int) float64 {
	x := make([]float64, n)
	copy(x, input)
	//nth.Element(floats(x), num)
	sort.Float64s(x)
	return x[num]
}

func findLargestIndices(output []int, num int, samples []float64, n int) error {
	if n < num+1 {
		return fmt.Errorf("n must be greater than num")
	}

	cutoff := nthElementImmutable(samples, n, num)
	count := 0

	for i := 0; i < n; i++ {
		if samples[i] < cutoff {
			output[count] = i
			count++
		}
	}
	if count < num {
		for i := 0; i < n; i++ {
			if samples[i] == cutoff {
				output[count] = i
				count++
				if count >= num {
					break
				}
			}
		}
	}
	sort.Ints(output)
	if count != num {
		return fmt.Errorf("count should equal num")
	}
	return nil
}

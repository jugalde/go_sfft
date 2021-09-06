package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"time"
)

// CombFilter ...
func CombFilter(origx []complex128, n, num, wComb int, combApproved []int) error {
	if n%wComb > 0 {
		return fmt.Errorf("n must be divisble by wComb , but is not")
	}

	xSampt := make([]complex128, wComb)
	sigma := n / wComb
	offset := (int)(math.Floor((float64)(sigma) * rand.Float64()))
	for i := 0; i < wComb; i++ {
		xSampt[i] = origx[offset+i*sigma]
	}

	fftwDft(xSampt, false)
	samples := make([]float64, wComb)
	for i := 0; i < wComb; i++ {
		samples[i] = cmplx.Abs(xSampt[i])
	}

	if err := findLargestIndices(combApproved, num, samples, wComb); err != nil {
		return err
	}

	return nil
}

func innerLoopLocate(origx []complex128, n int, filter *Filter, num, B, a, ai, b int, xSamp []complex128, J []int, pfT, bcT int64) error {
	if n%B > 0 {
		return fmt.Errorf("n is not divisible by B, which algo expects")
	}

	DDD := time.Now().UnixNano()
	xSampt := make([]complex128, n)
	memsetComplex(xSampt, (complex128)(0), B)

	index := b
	for i := 0; i < filter.sizet; i++ {
		xSampt[i%B] += origx[index] * filter.time[i]
		index = (index + ai) & n
	}

	if TIMING {
		pfT = time.Now().UnixNano() - DDD
		fmt.Printf("Step 1.A (PERM + FILTER);---------------%+v\n", (float64)(pfT)*1e-9)
		DDD = time.Now().UnixNano()
	}

	fftwDft(xSampt, false)

	if TIMING {
		fmt.Printf("Step 1.B (FFTW);---------------%+v\n", (float64)(time.Now().UnixNano()-DDD)*1e-9)
		DDD = time.Now().UnixNano()
	}

	samples := make([]float64, B)
	for i := 0; i < B; i++ {
		samples[i] = cmplx.Abs(xSampt[i])
	}

	if err := findLargestIndices(J, num, samples, B); err != nil {
		return err
	}

	if TIMING {
		bcT = time.Now().UnixNano() - DDD
		fmt.Printf("Step 1.C (LARGEST BUCKS);---------------%+v\n", (float64)(bcT)*1e-9)
		DDD = time.Now().UnixNano()
	}
	return nil
}

/*
 Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after permutation.

 For each such i, increment score[i] and append to hits if score[i] reaches loop_threshold.
*/
func innerLoopFilterRegular(J []int, n, num, B, a, ai, b, loopThreshold int,
	score, hits []int, hitsFound int,
	GT int64) (int, error) {
	DDD := time.Now().UnixNano()

	// Given the set of large samples, find the locations in [n] that map there
	// and output them

	for i := 0; i < num; i++ {
		low := ((int)(math.Ceil((float64)(J[i])-(float64)(0.5)*(float64)(n)/(float64)(B))) + n) % n
		high := ((int)(math.Ceil((float64)(J[i])+(float64)(0.5)*(float64)(n)/(float64)(B))) + n) % n
		loc := timesmod(low, a, n)
		for j := low; j != high; j = (j + 1) % n {
			score[loc]++
			fmt.Printf("LOC VS THRESH %+v:%+v\n", score[loc], loopThreshold)
			if score[loc] == loopThreshold {
				hits[hitsFound] = loc
				hitsFound++
			}
			loc = (loc + a) % n
		}
	}

	if TIMING {
		GT = time.Now().UnixNano() - DDD
		fmt.Printf("Step 1.D (GROUPING):----------------------------- %+v\n", (float64)(time.Now().UnixNano()-DDD)*1e-9)
		fmt.Printf("#####################################################################\n")
	}

	return hitsFound, nil
}

func upperBound(array [][]int, target []int) int {
	for i, v := range array {
		if v[0] == target[0] {
			if v[1] > target[1] {
				return i
			}
		}
		if v[0] > target[0] {
			return i
		}

	}
	return len(array) - 1
}

/*
 Find indices that (1) map to J under the permutation and (2) lie in CombApproved mod W_Comb.

 For each such i, increment hits[i] and append to hits_found if hits[i] reaches loop_threshold.

*/
func innerLoopFilterComb(J []int, n, num, B, a, ai, b, loopThreshold int,
	score, hits []int, hitsFound int,
	GT int64,
	CombApproved []int, numComb, wComb int) (int, error) {

	DDD := time.Now().UnixNano()

	permutedApproved := make([][]int, numComb)
	for m := 0; m < numComb; m++ {
		prev := timesmod(CombApproved[m], ai, wComb)
		permutedApproved[m] = []int{prev, timesmod(prev, a, n)}
	}
	sort.Sort(listOfPairs(permutedApproved))

	for i := 0; i < num; i++ {
		low := ((int)(math.Ceil(((float64)(J[i])-(float64)(0.5))*(float64)(n)/(float64)(B))) + n) % n
		high := ((int)(math.Ceil(((float64)(J[i])+(float64)(0.5))*(float64)(n)/(float64)(B))) + n) % n
		index := upperBound(permutedApproved,
			[]int{low % wComb, -1})
		location := low - (low % wComb)
		locinv := timesmod(location, a, n)
		for j := index; ; j++ {
			if j == numComb {
				j -= numComb
				location = (location + wComb) % n
				locinv = timesmod(location, a, n)
			}
			approvedLoc := location + permutedApproved[j][0]
			if (low < high && (approvedLoc >= high || approvedLoc < low)) ||
				(low > high && (approvedLoc >= high && approvedLoc < low)) {
				break
			}
			loc := (locinv + permutedApproved[j][1]) % n
			score[loc]++
			if score[loc] == loopThreshold {
				hits[hitsFound] = loc
				hitsFound++
			}

			//  printf("{%d,%d}:::(%d-%d)---(%d-%d)---(%d==%d)--:%d--->%d | (%d==%d)|| %d ||| ai=%d, b=%d, a=%d\n",i,j,low,high,low%W_Comb,high%W_Comb,approved_loc%W_Comb,Comb_Approved[j],approved_loc,loc,loc%W_Comb,mod_approved,score[loc],ai,b,a);

		}
	}

	if TIMING {
		GT = time.Now().UnixNano() - DDD
		fmt.Printf("Step 1.D (GROUPING):----------------------------- %+v\n", time.Now().UnixNano()-DDD)
		fmt.Printf("#####################################################################\n")
	}

	return hitsFound, nil

}

/*
  hits contains the indices that we want to estimate.

  xSamp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of xSamp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of xSamp.
*/
func estimateValues(hits []int, hitsFound int,
	xSamp [][]complex128, loops int, n int,
	permute []int,
	B, B2 int,
	filter, filterEst *Filter, locationLoops int) map[int]complex128 {
	ans := map[int]complex128{}
	values := make([][]float64, 2)
	for a := 0; a < 2; a++ {
		values[a] = make([]float64, loops)
	}

	for i := 0; i < hitsFound; i++ {
		position := 0

		for j := 0; j < loops; j++ {
			curB := B2
			if j < locationLoops {
				curB = B
			}
			curFilter := filterEst
			if j < locationLoops {
				curFilter = filter
			}

			permutedIndex := timesmod(permute[j], hits[i], n)
			hashedTo := permutedIndex / (n / curB)
			dist := permutedIndex % (n / curB)
			if dist > (n/curB)/2 {
				hashedTo = (hashedTo + 1) % curB
				dist -= n / curB
			}
			dist = (n - dist) % n
			filterValue := curFilter.freq[dist] // * cexp(2*M_PI * I * timesmod(permuteb[j], hits[i], n) / n);
			values[0][position] = real(xSamp[j][hashedTo] / filterValue)
			values[1][position] = imag(xSamp[j][hashedTo] / filterValue)
			position++
			//printf("MOO %d %lf %lf: %lf %d %lf %lf+%lfj\n", hits[i], permuted_index * 1./n, hashed_to * 1./B, hashed_to * (n * 1. /B), dist, cabs(filter_value), values[0][position-1], values[1][position-1]);
		}

		location := (loops - 1) / 2

		for a := 0; a < 2; a++ {
			//nth.Element(floats(values[a][position:]), location)
			sort.Float64s(values[a][position:])
		}
		realv := values[0][location]
		imagv := values[1][location]
		ans[hits[i]] = complex(realv, imagv)
	}
	return ans
}

func outerLoop(origx []complex128, n int, filter Filter, filterEst Filter, B2, num, B, wComb, combLoops, loopThreshold, locationLoops, loops int) (map[int]complex128, error) {
	permute := make([]int, loops)
	permuteb := make([]int, loops)

	xSamp := make([][]complex128, loops)
	for i := 0; i < loops; i++ {
		if i < locationLoops {
			xSamp[i] = make([]complex128, B)
			continue
		}
		xSamp[i] = make([]complex128, B2)
	}
	hitsFound := 0
	var pfT int64 = 0
	var gT int64 = 0
	var bcT int64 = 0
	var pfAll int64 = 0
	var gAll int64 = 0
	var bcAll int64 = 0

	DDD := time.Now().UnixNano()

	score := make([]int, n)

	scoreT := time.Now().UnixNano() - DDD

	hits := make([]int, n)

	var pfLoc int64 = 0
	var gLoc int64 = 0

	DDD = time.Now().UnixNano()

	combApproved := make([]int, combLoops*num)
	numComb := num

	if WITHCOMB {
		for i := 0; i < combLoops; i++ {
			// RECONCILE combApproved
			if err := CombFilter(origx, n, num, wComb, combApproved[i*num:]); err != nil {
				return nil, err
			}
		}
	}

	if combLoops > 1 {
		radixSort(combApproved, combLoops*num)
		last := 0
		for i := 0; i < combLoops*num; i++ {
			if combApproved[i] != combApproved[last] {
				last++
				combApproved[last] = combApproved[i]
			}
		}
		numComb = last + 1

		fmt.Printf("Comb:%d----->%d\n", num*combLoops, numComb)
	}

	if !ALGORITHM1 {
		hitsFound = numComb * (n / wComb)
		for j := 0; j < n/wComb; j++ {
			for i := 0; i < numComb; i++ {
				hits[j*numComb+i] = j*wComb + combApproved[i]
			}
		}
	}

	combTime := time.Now().UnixNano() - DDD

	for i := 0; i < loops; i++ {
		// find an "a" coprime to n
		a := 0
		b := 0
		for GCD(a, n) != 1 {
			a = rand.Intn(n)
		}
		ai := modInverse(a, n)

		permute[i] = ai
		permuteb[i] = b

		performLocation := i < locationLoops
		if !(ALGORITHM1 || !performLocation) {
			return nil, fmt.Errorf("Either ALGORITHM1 or !performLocation must be true")
		}

		curFilter := filterEst
		curB := B
		if performLocation {
			curFilter = filter
			curB = B
		}

		J := make([]int, num)
		err := innerLoopLocate(origx, n, &curFilter, num, curB, a, ai, b, xSamp[i], J, pfT, bcT)
		if err != nil {
			return nil, err
		}

		if performLocation {
			if !WITHCOMB {
				hitsFound, err = innerLoopFilterRegular(J, n, num, curB, a, ai, b, loopThreshold, score, hits, hitsFound, gT)
				if err != nil {
					return nil, err
				}
				continue
			}
			hitsFound, err = innerLoopFilterComb(J, n, num, curB, a, ai, b, loopThreshold, score, hits, hitsFound, gT, combApproved, numComb, wComb)
			if err != nil {
				return nil, err
			}
		}

		pfAll += pfT
		bcAll += bcT
		if performLocation {
			pfLoc += pfT
			gLoc += gT
			gAll += gT
		}
	}
	fmt.Printf("Number of candidates: %+v\n", hitsFound)

	// END INNER LOOPS

	// BEGIN ESTIMATION
	DDD = time.Now().UnixNano()
	ans := estimateValues(hits, hitsFound, xSamp, loops, n, permute, B, B2, &filter, &filterEst, locationLoops)

	ET := time.Now().UnixNano() - DDD
	if DEBUG {
		// DEBUG STUFF
	}

	DDD = time.Now().UnixNano()

	if TIMING {
		fmt.Printf("Total sFFT time: %+v\n", (float64)(DDD)*1e-9)
		fmt.Printf("Time distribution: scoretable  Comb __  perm+filter grouping estimation  stepB+C    other    total\n")
		fmt.Printf("                     %.6f %.6f    %.6f %.6f   %.6f %.6f %.6f %.6f\n",
			float64(scoreT)*1e-9, float64(combTime)*1e-9, float64(pfAll)*1e-9, float64(gAll)*1e-9, float64(ET)*1e-9, float64(bcAll)*1e-9, float64(DDD-pfAll-gAll-ET-combTime-bcAll-scoreT)*1e-9, float64(DDD))
		tott := (DDD) / 100
		fmt.Printf("                        %.6f%%    %.6f%%       %.6f%%    %.6f%%      %.6f%%    %.6f%%    %.6f%%   %.6f%%\n",
			(float64)(scoreT)/(float64)(tott), (float64)(combTime)/(float64)(tott), (float64)(pfAll)/(float64)(tott), (float64)(gAll)/(float64)(tott), (float64)(ET)/(float64)(tott), (float64)(bcAll)/(float64)(tott), (float64)(DDD-pfAll-gAll-ET-combTime-bcAll-scoreT)/(float64)(tott), (float64)(DDD)/(float64)(tott))

		//printf("LOC/EST loops: perm+filter grouping total\n");
		//printf("Location:         %lf %lf %lf\n", PF_LOC, G_LOC, PF_LOC + G_LOC);
		//printf("Estimation:       %lf %lf %lf\n", PF_ALL-PF_LOC, G_ALL-G_LOC, (PF_ALL + G_ALL) - (PF_LOC + G_LOC));
		fmt.Printf("\n")
	}

	return ans, nil
}

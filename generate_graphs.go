package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
)

/*
func runGraphExperiment(x []complex128, n int,
	lobefracLoc, toleranceLoc float64, bLoc,
	BLoc, BThresh, loopsLoc, loopsThresh int,
	lobefracEst, toleranceEst float64, bEst int,
	BEst, loopsEst, wComb, combLoops,
	repetitions int, fftwOpt bool,
	largeFreq []int, k int, xF []complex128,
	sffttime, fftwTime, errorPerEntry float64) error {
	if WITHCOMB {
		if BThresh >= wComb {
			return fmt.Errorf("BThresh needs to be less than wComb")
		}
	}
	if ALGORITHM1 {
		if BThresh >= BLoc {
			return fmt.Errorf("BThresh needs to be less than BLoc")
		}
		if loopsThresh > loopsLoc {
			return fmt.Errorf("loopsLoc needs to be less than loopsThresh")
		}
	}

	filterT, wLoc := makeDolphChebyshevT((float64)(lobefracLoc), (float64)(toleranceLoc))
	filter, err := makeMultipleT(filterT, wLoc, n, bLoc)
	if err != nil {
		return fmt.Errorf("Error hit on makeMultipleT: %+v", err)
	}

	filtertEst, wEst := makeDolphChebyshevT((float64)(lobefracEst), (float64)(toleranceEst))
	filterEst, err := makeMultipleT(filtertEst, wEst, n, bEst)
	if err != nil {
		return fmt.Errorf("Error hit on makeMultipleT: %+v", err)
	}
	fmt.Printf(" Window size: Location Filter : %+v; Estimation Filter : %+v;\n", wLoc, wEst)

	var filterNoise float64 = 0
	var filterNoiseEst float64 = 0

	for i := 0; i < 10; i++ {
		filterNoise = math.Max(filterNoise,
			math.Max(cmplx.Abs(filter.freq[n/2+i]),
				cmplx.Abs(filter.freq[n/2-i])))
		filterNoiseEst = math.Max(filterNoiseEst,
			math.Max(cmplx.Abs(filterEst.freq[n/2+i]),
				cmplx.Abs(filterEst.freq[n/2-i])))
	}


ans := map[int]complex128{}

		ans, err = outerLoop(x, n, filter, filterEst, BEst, BThresh, BLoc, wComb,
			combLoops, loopsThresh, loopsLoc, loopsLoc+loopsEst)
		if err != nil {
			fmt.Printf("Error on outerLoop: %+v\n", err)
			return
		}

	numCandidates := len(ans)
	candidates := make([][]float64, numCandidates)
	xFLarge := make([]complex128, n)
	ansLarge := make([]complex128, n)

	counter := 0

	for key, value := range ans {
		candidates[counter] = []float64{cmplx.Abs(value), (float64)(key)}
		counter++
	}

	//Enter ALL large frequences as zero
	for i := 0; i < k; i++ {
		xFLarge[largeFreq[i]] = xF[largeFreq[i]]
	}

	return nil
}*/

func graphExperiment() {
	nPtr := flag.Bool("N", false, "N")
	kPtr := flag.Bool("K", false, "K")
	sPtr := flag.Bool("S", false, "S")
	repetitionsPtr := flag.Int("R", 10, "repetitions")
	optPtr := flag.Bool("O", false, "fftwOpt")
	vPtr := flag.Bool("v", false, "verbose")
	withCombPtr := flag.Bool("W", false, "withComb")

	n := 4 * 128 * 8192
	k := 100
	var snr float64 = 100
	//var stdNoise float64 = 0
	graphType := 1
	FFTWOPT := false

	// Parse Flags
	flag.Parse()
	repetitions := *repetitionsPtr
	if *withCombPtr {
		WITHCOMB = true
	}

	if *optPtr {
		FFTWOPT = true
	}

	if *vPtr {
		VERBOSE = true
	}
	if *nPtr {
		graphType = 1
	}
	if *kPtr {
		graphType = 2
	}

	if *sPtr {
		graphType = 3
	}

	length := 1
	nVec := []int{8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216}
	kVec := []int{50, 100, 200, 500, 1000, 2000, 2500, 4000}
	snrVec := []float64{-20, -10, -7, -3, 0, 3, 7, 10, 20, 30, 40, 50, 60, 120}

	if graphType == 1 {
		//length = 12
		length = 7
	} else if graphType == 2 {
		//length = 8
		length = 7
	} else {
		//length = 14
		length = 7
	}

	sfftTime := make([][]float64, length)
	fftwTime := make([][]float64, length)
	sfftError := make([][]float64, length)

	for pp := 0; pp < length; pp++ {
		n = 4194304
		k = 50
		snr = snrVec[pp]
		BcstLoc, BcstEst, combCst, locLoops, estLoops, thresholdLoops, combLoops, toleranceLoc, toleranceEst := getExpermientVsNParameters(n, WITHCOMB)
		if graphType == 1 {
			n = nVec[pp]
			k = 50
			BcstLoc, BcstEst, combCst, locLoops, estLoops, thresholdLoops, combLoops, toleranceLoc, toleranceEst = getExpermientVsNParameters(n, WITHCOMB)
		} else if graphType == 2 {
			n = 4194304
			k = kVec[pp]
			BcstLoc, BcstEst, combCst, locLoops, estLoops, thresholdLoops, combLoops, toleranceLoc, toleranceEst = getExpermientVsKParameters(k, WITHCOMB)
		}

		// SET FILTER PARAMETERS

		if !(ALGORITHM1 || WITHCOMB) {
			fmt.Printf("ALOGRITHM1 or WITHCOMB must be true\n")
			return
		}
		BBLoc := (BcstLoc * math.Sqrt((float64)(n*k)/(math.Log2((float64)(n)))))
		BBEst := (BcstEst * math.Sqrt((float64)(n*k)/(math.Log2((float64)(n)))))

		lobefracLoc := 0.5 / (BBLoc)
		lobefracEst := 0.5 / (BBEst)

		bLoc := int(1.2 * 1.1 * ((float64)(n) / BBLoc))
		bEst := int(1.4 * 1.1 * ((float64)(n) / BBEst))

		BLoc := floorToPow2(BBLoc)
		BThresh := 2 * k
		BEst := floorToPow2(BBEst)

		wComb := floorToPow2(combCst * (float64)(n) / (float64)(BLoc))

		if graphType != 3 {
			fmt.Printf("Running SFFT and FFTW %d times for N=%d and K=%d\n", repetitions, n, k)
		} else {
			fmt.Printf("Running SFFT and FFTW %d times for N=%d and K=%d and SNR=%f dB\n", repetitions, n, k, snr)
		}

		var avgSfftTime float64 = 0
		var avgFftwTime float64 = 0
		var avgSfftError float64 = 0
		var itSfft float64
		var itFftw float64
		var itError float64

		for rr := 0; rr < repetitions; rr++ {

			x := make([]complex128, n)
			xF := make([]complex128, n)
			largeFreq := make([]int, k)

			// Randomized the None Zero Bins and Generate Time Domain Data
			for i := 0; i < k; i++ {
				largeFreq[i] = (int)(math.Floor(rand.Float64() * (float64)(n)))
				xF[largeFreq[i]] = 1.0 // Will ADD Random Phase and Amplitude Later.
			}

			copy(x, xF)
			fftwDft(x, true)

			if graphType == 3 {
				//stdNoise = math.Sqrt((float64)(k) / (2. * math.Pow(10, snr/10)))
				//	snrAchieved := AWGN(x, n, stdNoise)

				copy(xF, x)
				fftwDft(xF, true)
				for i := 0; i < n; i++ {
					xF[i] /= complex((float64)(n), 0)
				}

			}

			itSfft, itFftw, itError = runExperiment(x, n,
				lobefracLoc, toleranceLoc, bLoc,
				BLoc, BThresh, locLoops, thresholdLoops,
				lobefracEst, toleranceEst, bEst,
				BEst, estLoops, wComb, combLoops,
				repetitions, FFTWOPT, largeFreq, k, xF)

			avgSfftTime += itSfft
			avgFftwTime += itFftw
			avgSfftError += itError

		}
		if graphType == 1 {
			sfftTime[pp] = []float64{(float64)(n), avgSfftTime / (float64)(repetitions)}
			fftwTime[pp] = []float64{(float64)(n), avgFftwTime / (float64)(repetitions)}
		} else if graphType == 2 {
			sfftTime[pp] = []float64{(float64)(k), avgSfftTime / (float64)(repetitions)}
			fftwTime[pp] = []float64{(float64)(k), avgFftwTime / (float64)(repetitions)}
		} else {
			sfftError[pp] = []float64{snr, avgSfftError / (float64)(repetitions)}
		}
	}

	//fname := ""
	//persist := false
	//debug := true

	//	p, err := gnuplot.NewPlotter(fname, persist, debug)
	//	if err != nil {
	//		errString := fmt.Sprintf("** err: %v\n", err)
	//		panic(errString)
	//	}

	//plotTitle := "SFFT 1.0\nFFTW"
	//if WITHCOMB {
	//	plotTitle = "SFFT 2.0\nFFTW"
	//}
	if graphType == 3 {
		//x, y := arrOfArrToTwoArrs(sfftError)
		//err = p.PlotXY(x, y, plotTitle)
		fmt.Printf("ERROR IS %+v\n", sfftError)
	} else {
		//x1, y1 := arrOfArrToTwoArrs(sfftTime)
		//x2, y2 := arrOfArrToTwoArrs(fftwTime)
		//	err = p.PlotXY(append(x1, x2...), append(y1, y2...), plotTitle)
		fmt.Printf("SFFT IS %+v\n", sfftTime)
		fmt.Printf("FFTW IS %+v\n", fftwTime)
	}

	//	if err != nil {
	//		errString := fmt.Sprintf("** err: %v\n", err)
	//		panic(errString)
	//	}
	//
	//	defer func() {
	//		err := p.Close()
	//		if err != nil {
	//			errString := fmt.Sprintf("** err: %v\n", err)
	//			panic(errString)
	//		}
	//	}()

}

func arrOfArrToTwoArrs(arrOfArrs [][]float64) ([]float64, []float64) {
	x := []float64{}
	y := []float64{}
	for _, val := range arrOfArrs {
		x = append(x, val[0])
		y = append(y, val[1])
	}
	return x, y
}

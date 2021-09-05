package main

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"

	"github.com/cpmech/gosl/fun/fftw"
)

// TIMING ...
var TIMING = true

// ALGORITHM1 ...
var ALGORITHM1 = true

// WITHCOMB ...
var WITHCOMB = true

// DEBUG ...
var DEBUG = false

// FFTWOPT ...
var FFTWOPT = false

// FFTWFORWARD ...
var FFTWFORWARD = false

// FFTWMEASURE ...
var FFTWMEASURE = false

// FFTWESTIMATE ...
var FFTWESTIMATE = false

// VERBOSE ...
var VERBOSE = false

func runExperiment(x []complex128, n int,
	lobefracLoc, toleranceLoc float64, bLoc,
	BLoc, BThresh, loopsLoc, loopsThresh int,
	lobefracEst, toleranceEst float64, bEst int,
	BEst, loopsEst, wComb, combLoops,
	repetitions int, fftwOpt bool,
	largeFreq []int, k int, xF []complex128) {

	fmt.Printf("sFFT filter parameters for: n=%+v, k=%+v.\n", n, k)
	fmt.Printf("******************************************************************************\n")
	if WITHCOMB {
		fmt.Printf(" Comb Filter: loops: %+v mod: %+v/%+v\n", combLoops, BThresh, wComb)
	} else {
		fmt.Printf(" Comb Filter: none\n")
	}
	if ALGORITHM1 {
		fmt.Printf(" Location Filter: (numlobes=%+v, tol=%+v, b=%+v) B: %+v/%+v loops: %+v/%+v\n", (float64)(0.5)/(float64)(lobefracLoc), toleranceLoc, bLoc, BThresh, BLoc, loopsThresh, loopsLoc)
	} else {
		fmt.Printf(" Location Filter: none\n")
	}
	fmt.Printf(" Estimation Filter: (numlobes=%+v, tol=%+v, b=%+v) B: %+v loops: %+v\n", (float64)(0.5)/(float64)(lobefracEst), toleranceEst, bEst, BEst, loopsEst)
	fmt.Printf("\n")

	if WITHCOMB && BThresh >= wComb {
		fmt.Printf("If using comb, BThresh (%+v) must be less than wComb (%+v)\n", BThresh, wComb)
		return
	}
	if ALGORITHM1 {
		if BThresh >= BLoc {
			fmt.Printf("If using algo1, BThresh must be less than Bloc\n")
		}
		if loopsThresh > loopsLoc {
			fmt.Printf("If using algo1, loopsThresh must be less than loopsLoc\n")
		}
	}

	filterT, wLoc := makeDolphChebyshevT((float64)(lobefracLoc), (float64)(toleranceLoc))
	filter, err := makeMultipleT(filterT, wLoc, n, bLoc)
	if err != nil {
		fmt.Printf("Error hit on makeMultipleT: %+v\n", err)
		return
	}

	filtertEst, wEst := makeDolphChebyshevT((float64)(lobefracEst), (float64)(toleranceEst))
	filterEst, err := makeMultipleT(filtertEst, wEst, n, bEst)
	if err != nil {
		fmt.Printf("Error hit on makeMultipleT: %+v\n", err)
		return
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
	fmt.Printf(" Noise in filter: Location Filter : %+v; Estimation Filter %+v\n", filterNoise, filterNoiseEst)
	fmt.Printf("******************************************************************************\n\n")

	if DEBUG {
		// DO DEBUG STUFF
	}

	fmt.Printf("sFFT Results\n")
	fmt.Printf("******************************************************************************\n")

	ans := map[int]complex128{}

	for i := 0; i < repetitions; i++ {
		//reset_timer();
		ans, err = outerLoop(x, n, filter, filterEst, BEst, BThresh, BLoc, wComb,
			combLoops, loopsThresh, loopsLoc, loopsLoc+loopsEst)
		if err != nil {
			fmt.Printf("Error on outerLoop: %+v\n", err)
			return
		}
		fmt.Printf("????ANS %+v:%+v\n\n", ans, repetitions)
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

	//nth.Element(listOfFloats(candidates[0:numCandidates]), numCandidates-k)
	sort.Sort(listOfFloats(candidates[0:numCandidates]))
	fmt.Printf("??? %+v\n", candidates)
	for i := 0; i < k; i++ {
		key := candidates[numCandidates-k+i][1]
		ansLarge[(int)(key)] = ans[(int)(key)]
	}

	largeFound := 0
	FOUND := 0
	for i := 0; i < k; i++ {
		if _, ok := ans[largeFreq[i]]; ok {
			FOUND++
		}
		if ansLarge[(largeFreq[i])] != 0 {
			largeFound++
		}
	}

	//Estimate error as the difference of the K largest entries of x_f and ans)
	var ERROR float64 = 0
	for i := 0; i < n; i++ {
		ERROR += cmplx.Abs(ansLarge[i] - xFLarge[i])
	}

	//fmt.Printf("---------------CONSIDER K LARGEST ONLY ---------------------------------\n");
	///fmt.Printf("ERROR:\n")
	//fmt.Printf("K=%d; MISSED (estimation, result) = (%d, %d); L1 ERROR= %+v  (%+v per large frequency)\n", k, k-FOUND, k-large_found, ERROR, ERROR/k)
	fmt.Printf("******************************************************************************\n\n")

	fmt.Printf("FFTW Results\n")
	fmt.Printf("******************************************************************************\n")
	xtmp := make([]complex128, n)
	copy(xtmp, x)
	//reset_timer();

	var p *fftw.Plan1d

	if FFTWOPT {
		p = fftw.NewPlan1d(xtmp, FFTWFORWARD, FFTWMEASURE)
	} else {
		p = fftw.NewPlan1d(xtmp, FFTWFORWARD, FFTWESTIMATE)
	}

	//  fmt.Printf("Time to create FFTW plan: %+V\n", get_time());
	// reset_timer();
	for i := 0; i < repetitions; i++ {
		p.Execute()
	}
	//fmt.Printf("Time to run FFTW : %lf\n", get_time());
	// fftw_destroy_plan(p);
	p.Free()
	fmt.Printf("******************************************************************************\n\n")

}

func evaluateRuntime(n int, lobefrac, tolerance,
	lobefrac2, tolerance2 float64,
	num, B, B2, locationLoops,
	estLoops, loopThreshold, wComb,
	combLoops int) float64 {
	w := int((float64)(1/math.Pi) * (float64)(1) / (float64)(lobefrac) * math.Acosh((float64)(1.)/(float64)(tolerance)))
	w2 := int((float64)(1/math.Pi) * (float64)(1) / (float64)(lobefrac2) * math.Acosh((float64)(1.)/(float64)(tolerance2)))
	fmt.Printf("W ARE %+v:%+v\n", w, w2)

	if WITHCOMB {
		if num >= wComb {
			return -1
		}
	}
	if ALGORITHM1 {
		if num >= B {
			return -1
		}
		if loopThreshold > locationLoops {
			return -1
		}
	}
	if w > n || w2 > n {
		return -1
	}

	loops := locationLoops + estLoops
	mFrac := 1
	if WITHCOMB {
		mFrac = num * 1. / wComb
	}
	var projectedHits float64 = (float64)(n * mFrac)
	if ALGORITHM1 {
		projectedHits = binomialCdf((float64)(num)*((float64)(1.)/(float64)(B)-(float64)(1.)/(float64)(n)), locationLoops, loopThreshold)*(float64)(n)*(float64)(mFrac) + (float64)(num)/(float64)(2)
	}
	//XXX B2 for some, B for some
	var kEst float64 = (float64)(num) / (float64)(2)

	projectedNoiseOnK := (float64)(2) * binomialCdf(kEst*((float64)(1.)/(float64)(B2)-(float64)(1.)/(float64)(n))/(float64)(2), loops, (loops+1)/2) * kEst
	projectedErrorRate := (float64)(2)*binomialCdf(kEst*((float64)(1.)/(float64)(B2)-(float64)(1.)/(float64)(n))/(float64)(4), loops, (loops+1)/2)*(float64)(n*mFrac) + projectedNoiseOnK
	//double projected_error_rate = binomial_cdf((num/2) * (1. / B2 - 1./n), est_loops, (est_loops+1)/2) * (projected_hits - num/2);
	fmt.Printf("Projected error rate: %+v (%+v per large frequency)\n", projectedErrorRate, projectedNoiseOnK)

	pagesToSet := num * (n / B) * mFrac * locationLoops * 1024
	willArrayMemset := pagesToSet > n

	var constScorearray float64 = 1.8
	if n < (1 << 21) {
		constScorearray = 0.3
	}
	constPermfilt := 38.0
	constCombtime := 90.0
	var constEstimationtime float64 = 150
	if WITHCOMB {
		constEstimationtime = 140
	}
	var constGrouping float64 = 23
	var constGroupSort float64 = 0
	if WITHCOMB {
		constGroupSort = 30
	}
	constBplusctime := 41

	var timeScorearray float64 = 0
	willArrayMemsetNum := 1
	if willArrayMemset {
		willArrayMemsetNum = 0
		timeScorearray = constScorearray * (float64)(n)
	}

	timeComb := constCombtime * (float64)(wComb*combLoops)

	timePermfilt := constPermfilt * (float64)(w*1.*locationLoops+w2*1.*estLoops)
	timeGrouping := ((float64)(locationLoops)*(constGrouping*(float64)(num*(n/B)*mFrac)+constGroupSort*(float64)(num)*math.Log(float64(num))) + constScorearray*(float64)(willArrayMemsetNum*pagesToSet))
	timeEstimation := constEstimationtime * projectedHits * (float64)(loops)
	var timeBplusc float64 = (float64)(constBplusctime * (locationLoops*B + estLoops*B2))
	timeTotal := timeScorearray + timeComb + timePermfilt + timeGrouping + timeEstimation + timeBplusc

	return timeTotal
}

func main() {
	nPtr := flag.Int("N", 4*128*8192, "N")
	kPtr := flag.Int("K", 100, "K")
	repetitionsPtr := flag.Int("R", 1, "repetitions")

	BcstLocPtr := flag.Float64("B", 2., "BcstLoc")
	BcstEstPtr := flag.Float64("E", 0.2, "BcstEst")
	BcstEst := *BcstEstPtr
	combCstPtr := flag.Float64("M", 16., "combCst")

	locLoopsPtr := flag.Int("l", 3, "locLoops")

	estLoopsPtr := flag.Int("L", 12, "estLoops")
	thresholdLoopsPtr := flag.Int("r", 2, "thresholdLoops")
	combLoopsPtr := flag.Int("m", 1, "combLoops")
	snrPtr := flag.Float64("S", 1000000000, "snr")
	algPtr := flag.Bool("A", false, "NotAlgortihm1")

	optPtr := flag.Bool("O", false, "fftwOpt")

	vPtr := flag.Bool("v", false, "verbose")

	FFTWOPT := false
	toleranceLocPtr := flag.Float64("t", 1.e-6, "toleranceLoc")
	toleranceEstPtr := flag.Float64("e", 1.e-8, "toleranceEst")
	simulatePtr := flag.Bool("s", false, "simulate")

	// Parse Flags
	flag.Parse()

	n := *nPtr
	k := *kPtr
	repetitions := *repetitionsPtr
	BcstLoc := *BcstLocPtr

	combCst := *combCstPtr
	if combCst > 0 {
		WITHCOMB = true
	}

	locLoops := *locLoopsPtr
	if locLoops == 0 {
		ALGORITHM1 = false
	}

	estLoops := *estLoopsPtr
	thresholdLoops := *thresholdLoopsPtr
	combLoops := *combLoopsPtr
	snr := *snrPtr
	stdNoise := math.Sqrt((float64)(k) / (2. * snr))
	if *algPtr {
		ALGORITHM1 = false
		locLoops = 0
	}
	if *optPtr {
		FFTWOPT = true
	}

	if *vPtr {
		VERBOSE = true
	}
	toleranceLoc := *toleranceLocPtr

	toleranceEst := *toleranceEstPtr
	simulate := *simulatePtr

	n = floorToPow2((float64)(n))
	if !(ALGORITHM1 || WITHCOMB) {
		fmt.Printf("ALOGRITHM1 or WITHCOMB must be true\n")
		return
	}

	BBLoc := (BcstLoc * math.Sqrt((float64)(n*k)/(math.Log2((float64)(n)))))
	BBEst := (BcstEst * math.Sqrt((float64)(n*k)/(math.Log2((float64)(n)))))

	lobefracLoc := 0.5 / (BBLoc)
	lobefracEst := 0.5 / (BBEst)

	bLoc := int(1.2 * 1.1 * ((float64)(n) / BBLoc))
	//b_loc = 1;
	//real_t BB2 = (unsigned) (Bcst2*sqrt((double)n*k/(log2(n))));
	bEst := int(1.4 * 1.1 * ((float64)(n) / BBEst))

	BLoc := floorToPow2(BBLoc)
	BThresh := 2 * k
	BEst := floorToPow2(BBEst)

	wComb := floorToPow2(combCst * (float64)(n) / (float64)(BLoc))

	fmt.Printf("\n\nRUNNING EXPERIMENT: n=%+v, k=%+v.\n", n, k)

	fmt.Printf("\n\nSimulation:\n")
	fmt.Printf("******************************************************************************\n")
	fmt.Printf("Expected running time: %+v\n",
		evaluateRuntime(n, lobefracLoc, toleranceLoc, lobefracEst, toleranceEst,
			BThresh, BLoc, BEst, locLoops, estLoops, thresholdLoops, wComb, combLoops)*1e-9)
	fmt.Printf("******************************************************************************\n")
	if simulate {
		return
	}

	x := make([]complex128, n)

	//srand(17);
	// srand48( time(NULL) ^ (getpid() * 171717));

	//Randomized the None Zero Bins

	xF := make([]complex128, n)
	copy(xF, x)

	largeFreq := make([]int, k)

	// fmt.Printf("LARGE BINS:");

	//  for(int i = 0; i < repetitions; i++){

	for i := 0; i < k; i++ {
		largeFreq[i] = (int)(math.Floor(rand.Float64() * (float64)(n)))
		xF[largeFreq[i]] = 1.0 // Will ADD Random Phase and Amplitude Later.
		//      fmt.Printf("%d, ",LARGE_FREQ[i]);
	}
	fmt.Printf("\n")

	fftwDft(xF, true)

	//ADDED NOISE
	snrAchieved := AWGN(x, n, stdNoise)
	if stdNoise != 0 {
		fmt.Printf("SNR = %+v / %+v dB \n\n", snrAchieved, 10*math.Log10(snrAchieved))
	}

	copy(x, xF)
	fftwDft(x, true)
	for i := 0; i < n; i++ {
		xF[i] /= complex((float64)(n), 0)
	}

	runExperiment(x, n,
		lobefracLoc, toleranceLoc, bLoc,
		BLoc, BThresh, locLoops, thresholdLoops,
		lobefracEst, toleranceEst, bEst,
		BEst, estLoops, wComb, combLoops,
		repetitions, FFTWOPT, largeFreq, k, xF)

	return
}

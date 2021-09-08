package main

func getExpermientVsNParameters(N int, WITHCOMB bool) (bcstLoc, bcstEst, combCst float64, locLoops, estLoops, thresholdLoops, combLoops int, toleranceLoc, toleranceEst float64) {
	bcstLoc = 1
	bcstEst = 1
	combCst = 2
	combLoops = 1
	locLoops = 3
	estLoops = 16
	thresholdLoops = 3
	toleranceLoc = 1.e-8
	toleranceEst = 1.e-8

	if WITHCOMB {

		switch N {
		case 8192:
			bcstLoc = 2
			bcstEst = 2
			combCst = 32
			combLoops = 8
			estLoops = 16
			locLoops = 7
			thresholdLoops = 6
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 16384:
			bcstLoc = 4
			bcstEst = 4
			combCst = 32
			combLoops = 8
			estLoops = 10
			locLoops = 6
			thresholdLoops = 5
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 32768:
			bcstLoc = 4
			bcstEst = 2
			combCst = 64
			combLoops = 4
			estLoops = 8
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 65536:
			bcstLoc = 4
			bcstEst = 2
			combCst = 128
			combLoops = 6
			estLoops = 10
			locLoops = 4
			thresholdLoops = 2
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 131072:
			bcstLoc = 1
			bcstEst = 1
			combCst = 8
			combLoops = 2
			estLoops = 12
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 262144:
			bcstLoc = 1
			bcstEst = 1
			combCst = 8
			combLoops = 2
			estLoops = 14
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 524288:
			bcstLoc = 0.5
			bcstEst = 0.5
			combCst = 8
			combLoops = 1
			estLoops = 10
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 1048576:
			bcstLoc = 0.5
			bcstEst = 0.5
			combCst = 8
			combLoops = 2
			estLoops = 12
			locLoops = 4
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 2097152:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 8
			combLoops = 1
			estLoops = 10
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 4194304:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 8
			combLoops = 1
			estLoops = 8
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 8388608:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 8
			combLoops = 1
			estLoops = 8
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 16777216:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 16
			combLoops = 1
			estLoops = 8
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		}
	} else {
		switch N {
		case 8192:
			bcstLoc = 2
			bcstEst = 2
			combCst = 1
			combLoops = 1
			estLoops = 16
			locLoops = 7
			thresholdLoops = 6
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 16384:
			bcstLoc = 4
			bcstEst = 4
			combCst = 1
			combLoops = 1
			estLoops = 10
			locLoops = 6
			thresholdLoops = 5
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 32768:
			bcstLoc = 4
			bcstEst = 2
			combCst = 1
			combLoops = 1
			estLoops = 8
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 65536:
			bcstLoc = 4
			bcstEst = 2
			combCst = 1
			combLoops = 1
			estLoops = 8
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-8
			toleranceEst = 1e-8

			break

		case 131072:
			bcstLoc = 2
			bcstEst = 1
			combCst = 1
			combLoops = 1
			estLoops = 10
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 262144:
			bcstLoc = 2
			bcstEst = 0.5
			combCst = 1
			combLoops = 1
			estLoops = 14
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 524288:
			bcstLoc = 1
			bcstEst = 0.5
			combCst = 1
			combLoops = 1
			estLoops = 12
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 1048576:
			bcstLoc = 2
			bcstEst = 0.5
			combCst = 1
			combLoops = 1
			estLoops = 12
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 2097152:
			bcstLoc = 2
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 15
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 4194304:
			bcstLoc = 4
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 10
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 8388608:
			bcstLoc = 2
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 8
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		case 16777216:
			bcstLoc = 4
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 8
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1e-8

			break

		}
	}

	return

}

func getExpermientVsKParameters(K int, WITHCOMB bool) (bcstLoc, bcstEst, combCst float64, locLoops, estLoops, thresholdLoops, combLoops int, toleranceLoc, toleranceEst float64) {
	bcstLoc = 1
	bcstEst = 1
	combCst = 2
	combLoops = 1
	locLoops = 3
	estLoops = 16
	thresholdLoops = 3
	toleranceLoc = 1.e-8
	toleranceEst = 1.e-8

	if WITHCOMB {

		switch K {
		case 50:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 16
			combLoops = 1
			estLoops = 10
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1.0e-8

			break

		case 100:
			bcstLoc = 0.5
			bcstEst = 0.2
			combCst = 16
			combLoops = 1
			estLoops = 12
			locLoops = 4
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1.0e-8

			break

		case 200:
			bcstLoc = 0.5
			bcstEst = 0.5
			combCst = 32
			combLoops = 1
			estLoops = 8
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 0.5e-8

			break

		case 500:
			bcstLoc = 0.5
			bcstEst = 0.5
			combCst = 64
			combLoops = 1
			estLoops = 10
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 0.5e-8

			break

		case 1000:
			bcstLoc = 1
			bcstEst = 1
			combCst = 128
			combLoops = 3
			estLoops = 12
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 0.5e-8

			break

		case 2000:
			bcstLoc = 1
			bcstEst = 1
			combCst = 512
			combLoops = 3
			estLoops = 16
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-7
			toleranceEst = 0.2e-8

			break

		case 2500:
			bcstLoc = 1
			bcstEst = 1
			combCst = 512
			combLoops = 3
			estLoops = 16
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-7
			toleranceEst = 0.2e-8

			break

		case 4000:
			bcstLoc = 1
			bcstEst = 2
			combCst = 512
			combLoops = 3
			estLoops = 14
			locLoops = 8
			thresholdLoops = 7
			toleranceLoc = 1e-8
			toleranceEst = 0.5e-8

			break

		}
	} else {
		switch K {
		case 50:
			bcstLoc = 4
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 10
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1.0e-8

			break

		case 100:
			bcstLoc = 2
			bcstEst = 0.2
			combCst = 1
			combLoops = 1
			estLoops = 12
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 1.0e-8

			break

		case 200:
			bcstLoc = 4
			bcstEst = 0.5
			combCst = 1
			combLoops = 1
			estLoops = 10
			locLoops = 3
			thresholdLoops = 2
			toleranceLoc = 1e-6
			toleranceEst = 0.5e-8

			break

		case 500:
			bcstLoc = 2
			bcstEst = 1
			combCst = 1
			combLoops = 1
			estLoops = 12
			locLoops = 4
			thresholdLoops = 3
			toleranceLoc = 1e-6
			toleranceEst = 0.5e-8

			break

		case 1000:
			bcstLoc = 2
			bcstEst = 1
			combCst = 1
			combLoops = 1
			estLoops = 12
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-6
			toleranceEst = 1.0e-8

			break

		case 2000:
			bcstLoc = 2
			bcstEst = 1
			combCst = 1
			combLoops = 1
			estLoops = 16
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-7
			toleranceEst = 0.5e-8

			break

		case 2500:
			bcstLoc = 2
			bcstEst = 1
			combCst = 1
			combLoops = 1
			estLoops = 16
			locLoops = 5
			thresholdLoops = 4
			toleranceLoc = 1e-7
			toleranceEst = 0.5e-8

			break

		case 4000:
			bcstLoc = 2
			bcstEst = 2
			combCst = 1
			combLoops = 1
			estLoops = 14
			locLoops = 6
			thresholdLoops = 5
			toleranceLoc = 1e-8
			toleranceEst = 1.0e-8

			break

		}
	}

	return

}

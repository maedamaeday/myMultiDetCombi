# calculate detection performance with multiple detector results combined
-----

## contents
 - myMultiDetCombi.py : class definition to calculate combined performance, to be imported by the main script
 - findBestDetCombi.py : main script.
 - draw2oo3.C : A [ROOT](https://cern.root.ch) macro to show performance improvement by 2-out-of-3, 3-out-of-4 and 2-out-of-4 logics.

## how to run
### calculate combined performance
Run with a list of false positive (FP) rates. The false negative (FN) rate is fixed to 0.01, which can be modified by ```-n``` option. Then, performance (combined false positive rates, false negative rates and F-values for each condition) with M-out-of-N (MooN) method and likelihood method. Here is an example:
```
python findBestDetCombi.py -p 0.01 0.02 0.1 0.3 0.5
input false positive rates : [0.01, 0.02, 0.1, 0.3, 0.5]
input false negative rate : 0.01 (common for all the detectors)
****M-out-of-N logic****
(4, 5) : (0.0004840000000000001, 0.0009801495999999998, 0.9992678636138077)
(3, 4) : (0.0009620000000000002, 0.00059203, 0.99922295075394)
(2, 3) : (0.0031600000000000005, 0.00029800000000000003, 0.9982689486922889)
best condition : N = 5, M = 4, FP = 0.0004840000000000001, FN = 0.0009801495999999998
****likelihood methods****
(0, 5) : (0.000581, 0.0006890599, 0.9993649671289094)
(0, 4) : (0.0009620000000000002, 0.00059203, 0.99922295075394)
(0, 3) : (0.0031600000000000005, 0.00029800000000000003, 0.9982689486922889)
(0, 2) : (0.01, 0.01, 0.99)
best condition : N = 5, FP = 0.000581, FN = 0.0006890599
```
where each tuple in the right of colon(:) means (M,N). M=0 in case of the likelihood method. Less N values than the number of elements in the given FP list are cases where detector(s) with the worst FPs are removed.
Tuples in the left of a colon show combined detector performance. Three values are FP rate, FN rate and F-value (harmonic average of 1-FP and 1-FN), respectively. The "best" performance is determined based on the F-value.

### calculate likelihood
Use ```findBestDetCombi.py``` with ```-o``` option followed by a list of observation for each detector. Observation can be given by 1 (positive detection) or 0. A list of FPs is also needed. The order of observation result must be the same with that of the FP list. Here is an example
```
python findBestDetCombi.py -p 0.01 0.02 0.1 0.3 0.5 -o 0 0 1 0 1
(13.835611229671274, 3.3826602606637444)
```
where the first value is a negative log likelihood with positive assumption and the second value is for negative assumption. Smaller value is more likely.

### show combined performance M-out-of-N logic (where N=3 or 4)
You need to install [ROOT](https://cern.root.ch), which is a C++ library developed mainly for particle physics data analysis. Then, you can produce a plot by running the following command:
```
root -l draw2oo3.C+
```
The graph shows ratios of combined FP (or FN) rate to the original ones.

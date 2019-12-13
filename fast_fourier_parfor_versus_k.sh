#!/bin/bash
# The above shows that it is a bash file

KS="2 4 6 8 16 32 64"

for K in $KS; do 
	export HPCE_FFT_LOOP_K=${K}
	bin/time_fourier_transform hpce.satish.fast_fourier_transform_parfor 0 30 "${HPCE_FFT_LOOP_K}, " > dumploop_${HPCE_FFT_LOOP_K}.csv
done 
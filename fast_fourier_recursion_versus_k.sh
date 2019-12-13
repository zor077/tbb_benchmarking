#!/bin/bash
# The above shows that it is a bash file

KS="2 4 8 16 32 64"

for K in $KS; do 
	export HPCE_FFT_RECURSION_K=${K}
	bin/time_fourier_transform hpce.satish.fast_fourier_transform_taskgroup 0 30 "${HPCE_FFT_RECURSION_K}, " > dump_rec_${HPCE_FFT_RECURSION_K}.csv
done 
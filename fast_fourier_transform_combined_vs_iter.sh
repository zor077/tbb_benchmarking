#!/bin/bash
# The above shows that it is a bash file

export HPCE_FFT_RECURSION_K=1
export HPCE_FFT_LOOP_K=1

bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dump_RK1_LK1.csv

LKS="2 4 8 16 32 64 128 256"

export HPCE_FFT_RECURSION_K=1
for LK in $LKS; do
	export HPCE_FFT_LOOP_K=${LK}
	bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dumpRK1_${HPCE_FFT_LOOP_K}.csv
done 

RKS="2 4 8 16 32 64 128 256"

export HPCE_FFT_LOOP_K=1
for RK in $RKS; do
	export HPCE_FFT_RECURSION_K=${RK}
	bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dump_${HPCE_FFT_RECURSION_K}_LK1.csv
done 

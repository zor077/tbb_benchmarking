#!/bin/bash
# The above shows that it is a bash file

export HPCE_FFT_RECURSION_K=1
export HPCE_FFT_LOOP_K=1

bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dumpK${HPCE_FFT_RECURSION_K}_LK${HPCE_FFT_LOOP_K}.csv

export HPCE_FFT_RECURSION_K=1
export HPCE_FFT_LOOP_K=128

bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dump_RK${HPCE_FFT_RECURSION_K}_LK${HPCE_FFT_LOOP_K}.csv

export HPCE_FFT_RECURSION_K=256
export HPCE_FFT_LOOP_K=1

bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dumpRK${HPCE_FFT_RECURSION_K}_LK${HPCE_FFT_LOOP_K}.csv

export HPCE_FFT_RECURSION_K=256
export HPCE_FFT_LOOP_K=128

bin/time_fourier_transform hpce.satish.fast_fourier_transform_combined 0 30 "${HPCE_FFT_RECURSION_K}, ${HPCE_FFT_LOOP_K}, " > dumpRK${HPCE_FFT_RECURSION_K}_LK${HPCE_FFT_LOOP_K}.csv

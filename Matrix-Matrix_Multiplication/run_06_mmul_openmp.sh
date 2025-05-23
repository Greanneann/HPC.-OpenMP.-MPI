#! /bin/bash

make clean
make openmp || exit 1
PROGRAM="./main_openmp"

# Iteration range
N_SEQ=$(seq 4 10)

N_=""
for N in $N_SEQ; do
	N_+=" $((2**N))"
	#N_+=" $((512*N))"
done

# Variants
VARIANT_="22 32 24 34 35"
# VARIANT_="22 32 35"
# VARIANT_="35"
# VARIANT_="35 34"

source benchmark_base.sh

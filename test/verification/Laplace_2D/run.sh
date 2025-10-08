#!/usr/bin/env bash

export FI_PROVIDER=tcp mpirun

for DIM in 16 32 64 128 256 512
do
    mkdir -p "${DIM}"
    
    cd "${DIM}"
    
    rm -rf plt* chk*
    mpirun -np 8 ../vidyut2d.llvm.MPI.ex ../inputs2d max_step=10 amr.n_cell="${DIM}" "${DIM}" 1
    
    cd ..
done

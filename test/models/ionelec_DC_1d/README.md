# Linear advection test

This case tests linear advection with the 5th order WENO scheme and 
verifies its accuracy. 

### Build instructions

make sure $AMREX_HOME is set to your clone of amrex
`$ export AMREX_HOME=/path/to/amrex`

To build a serial executable with gcc do
`$ make -j COMP=gnu`

To build a serial executable with clang++ do
`$ make -j COMP=llvm`

To build a parallel executable with gcc do
`$ make -j COMP=gnu USE_MPI=TRUE`

To build a parallel executable with gcc, mpi and cuda
`$ make -j COMP=gnu USE_CUDA=TRUE USE_MPI=TRUE`

### Run instructions

Use the `plasjob` script and then the `postprocess.sh` script to run 
multiple 1d cases and verify order of convergence.

Run with inputs2d to see how advection of a Gaussian feature is 
captured with AMR. `$ mpirun -n 1 ./*.ex inputs2d`

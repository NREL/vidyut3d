# Linear advection test

This case solves a coupled diffusion and poisson problem 
with two species

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

Runs can be done with any of the input files:
`inputs_x`,`inputs_y`,`inputs_z` to test if numerics 
are implemented correctly along each direction, e.g.
`$ mpirun -n 1 ./*.ex inputs_x`

Testing can be done also with the cell-mask feature 
using `inputs_ib` file

to view solution and compare with analytic solution, use 
python script:`python verify_spec.py plt00001`

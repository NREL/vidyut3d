# Gaseous Electronics Conference Radio-frequency cell

This case simulates GEC RF cell at 100 mTorr with argon 
plasma chemistry. This run takes about 4 mins to do
a single RF cycle with 64 processors. To get to steady state, 
about 500 cycles is required. This case uses the
cell masking feature in the 2D axisymmetric mode.

### Build instructions

make sure $AMREX_HOME is set to your clone of amrex
`$ export AMREX_HOME=/path/to/amrex`

To build a serial executable with gcc do
`$ make -j COMP=gnu`

To build a serial executable with clang++ do
`$ make -j COMP=llvm`

To build a parallel executable with gcc do
`$ make -j COMP=gnu USE_MPI=TRUE`

### Run instructions

Run with inputs2d 
`$ mpirun -n 64 ./*.ex inputs2d`

# Simple ion-electron plasma

This case simulates a DC discharge between two parallel plates.
It captures the cathode sheath formation and associated electron,ion and 
electron temperature distributions with a simple ion-electron plasma with 
just one ionization reaction.

### Build instructions

make sure $AMREX_HOME is set to your clone of amrex
`$ export AMREX_HOME=/path/to/amrex`

If you are copying this case folder elsewhere then
make sure $VIDYUT_DIR is set to your clone of vidyut
`$ export VIDYUT_DIR=/path/to/vidyut`

To build a serial executable with gcc do
`$ make -j COMP=gnu`

To build a serial executable with clang++ do
`$ make -j COMP=llvm`

To build a parallel executable with gcc do
`$ make -j COMP=gnu USE_MPI=TRUE`

To build a parallel executable with gcc, mpi and cuda
`$ make -j COMP=gnu USE_CUDA=TRUE USE_MPI=TRUE`

### Run instructions
This case can be built with `DIM=1`, `DIM=2`, or `DIM=3`, 1d being
the quickest to run.
`$ mpirun -n 4 ./*.ex inputs1d`

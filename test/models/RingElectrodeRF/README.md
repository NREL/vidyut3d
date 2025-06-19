# RF discharge with ring electrodes

This case simulates He discharge in a cylindrical tube
with ring electrodes. This case includes the thin
dielectric layer boundary condition along with surface charge solves.

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

`$ mpirun -n 16 ./*.ex inputs2d`

<img src="https://github.com/user-attachments/assets/4f860eea-af13-4b20-b438-922bd87a8057" width=600>

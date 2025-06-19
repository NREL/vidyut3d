# RF discharge with two offset electrodes

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
To run do
`$ mpirun -n 8 ./*.ex inputs2d`
It takes about 3 mins to do a single RF cycle.
See below, a picture of electron temperature, during the 3rd cycle.
<img src="https://github.com/user-attachments/assets/aa9530e9-d401-438d-ba8c-f7ba89dc80b6" width=900>

# Ar Streamers at atmospheric pressure

This case simulates propagation and interaction of 5 cathode directed streamers
initiated by seed charges between two parallel plates.

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

This is a 3D run with 4 AMR levels.
`mpirun -n 32 ./*.ex inputs3d`
Over about 12 ns, the streamer in the middle propagates faster
and the ones in the periphery tend to merge with the middle one.

<img src="https://github.com/user-attachments/assets/e2f4414f-7434-4255-8596-6b505d91633d" width=500>



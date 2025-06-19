# Ar-H2 Streamers at atmospheric pressure

This case simulates propagation and interaction of 14 cathode directed streamers
initiated by seed charges between two parallel plates with Ar-H2 plasma
chemistry. Some interesting things like formation of Argonium ions and 
H2 vibrational excited states can be seen in this case.

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

This is quite a large 3D run. With 4 AMR levels, the 
number of cells can approach 300 million.
This run finishes in approximately 3 hours with 200 NVIDIA-H100 GPUs
`mpirun -n 200 ./<gpu-ex> inputs3d`

<img src="https://github.com/user-attachments/assets/1acbc367-c4e0-478c-80dd-b81caf5da95d" width=500>


# Diffusion-Poisson test

This case solves a coupled diffusion and Poisson problem 
with two species:


$$\frac{ds_1}{dt}+\frac{d (v s_1)}{dx}=D_1\frac{d^2 s_1}{d x^2} \quad s_1(0)=0 \quad s_1(1)=1$$
$$\frac{ds_2}{dt}=\frac{d^2 s_2}{d x^2}+1.0 \quad s_2(0)=s_2(1)=0$$
$$\frac{d2\phi}{dx^2}=0 \quad \phi(0)=1 \quad \phi(1)=0$$
$$D_1=0.1 \quad v=E=-\frac{d\phi}{dx}$$

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
This case can be built in 1D,2D or 3D. Set the DIM variable 
accordingly when building the executable, 
e.g. `make -j COMP=gnu USE_MPI=TRUE DIM=1`

Runs can be done with any of the input files, inputs_x`,`inputs_y`,`inputs_z` to test if numerics 
are implemented correctly along each direction, e.g.
`$ mpirun -n 1 ./*.ex inputs_x` 
Note that at least a 
2D build is required to use `inputs_y` and a 3D build for `inputs_z``

Testing can be done also with the cell-mask feature 
using `inputs_ib` file

to view solution and compare with analytic solution, use 
python script:`python verify_spec.py plt00001`

<img src="https://github.com/user-attachments/assets/6c6fc5f7-7f02-4ec3-956c-a3fb0ae9f3ad" width="600">

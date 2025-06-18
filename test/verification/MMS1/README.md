# Method-of-Manufactured solutions (MMS) test 1

This case solves a coupled system that includes the species-Poisson numerical 
coupling in the plasma fluid model

$$\frac{dn}{dt}+\frac{d (\mu n E)}{dx}=\frac{d}{dx}\left(D\frac{dn}{dx}\right) + \left(\frac{5 x^4}{3} - \frac{x}{6} -2\right) \quad n(0)=0 \quad n(1)=1$$
$$\frac{d^2\phi}{dx^2}=n \quad \phi(0)=1 \quad \phi(1)=0$$
$$\mu=-1 \quad D=1 \quad E=-\frac{d\phi}{dx}$$

The exact solution to this system is:

$$n=x^2 \quad \phi=(x^4-x)/12$$

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

Runs can be done with any of the input files, `inputs_x`,`inputs_y`,`inputs_z`  
to test if numerics are implemented correctly along each direction, e.g.
`$ mpirun -n 1 ./*.ex inputs_x` 
Note that at least a 
2D build is required to use `inputs_y` and a 3D build for `inputs_z`.

Testing can be done also with the cell-mask feature 
using `inputs_ib_x/y/z` files.

to view solution and compare with analytic solution, use 
python script:`python verify_spec.py plt00001`

To calculate L2 norm error and also get the solution, use 
python script: `python get_L2norm_error.py plt00001 0` (use 0,1,2 for 
x,y or z direction for the second argument here).

<img src="https://github.com/user-attachments/assets/94a02372-ac25-4f4d-b3bf-aca7ac93c2c5" width="600">

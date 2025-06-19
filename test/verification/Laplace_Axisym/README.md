# Laplace solve with axisymmetric coordinate system

This case tests the Laplace equation solution in axisymmetric coordinate system.

$$\nabla^2\phi=0 \quad \frac{d^2\phi}{dx^2}+\frac{1}{r}\frac{d\phi}{dr}=0$$
$$\phi(r=R_{min})=\phi_1 \quad \phi(r=R_{max})=\phi_2$$

The exact solution for this equation is

$$\phi(r)=\frac{1}{\log\left(\frac{R_{max}}{R_{min}}\right)}\left(\phi_2 \log\left(\frac{r}{R_{min}}\right) + \phi_1 \log\left(\frac{R_{max}}{r}\right)\right)$$

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

Use the `run.sh` script to run an array of cases with Dirichlet and Neumann boundary conditions.
You will need the `fextract` executable from amrex, which can be built from within 
https://github.com/AMReX-Codes/amrex/tree/development/Tools/Plotfile

Alternatively, you can just do `mpirun -n 1 ./*.ex inputs2d`

<img src="https://github.com/user-attachments/assets/efd361a2-f2e1-4da5-be2f-9066374d037d" width=500>

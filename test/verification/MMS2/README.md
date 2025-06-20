# Method-of-Manufactured solutions (MMS) test 2

This case solves all equations in the plasma model with MMS source terms:

$$\frac{dn_e}{dt}+\frac{d \Gamma_e}{dx}=k_i n_e + S_e$$
$$\frac{dn_i}{dt}+\frac{d \Gamma_i}{dx}=k_i n_e + S_i$$
$$\frac{d^2\phi}{dx^2}=\frac{e(n_e-n_i)}{\epsilon_0} \quad E=-\frac{d\phi}{dx}$$
$$\frac{dE_e}{dt}+\frac{d \Gamma_E}{dx}=-e \Gamma_e E - k_i n_e E_i - \frac{3}{2} n_e k_B T_e \nu \frac{2 m_e}{m_h} + S_\epsilon$$
$$\Gamma_e= \mu_e n_e E - D_e \frac{dn_e}{dx}$$
$$\Gamma_i= \mu_i n_i E - D_i \frac{dn_i}{dx}$$
$$\Gamma_\epsilon= \frac{5}{3}\left( \mu_e E_\epsilon E - D_e \frac{dE_\epsilon}{dx}\right)$$

These values are assumed for various reaction and transport parameters where e and $$m_p$$ are 
electronic charge and proton mass, respectively:

$$k_i=5.0~\exp\left(-\frac{2.0}{k_BTe}\right)$$
$$\mu_e=-1.0 \quad D_e=1.0 \quad \mu_i=0.5 \quad D_i=0.5$$
$$E_i=\frac{4}{e} \quad \nu=10000.0 \quad m_h=4 m_p$$

We assume the exact solution to this system is:

$$n_e=\frac{x^2}{\alpha} + n_0$$
$$n_i=\frac{x^2}{2.0 \alpha}+n_0$$
$$\phi=\frac{1}{24}\left(x^4-x\right)$$
$$E_\epsilon=\frac{x^2}{\alpha}+n_0$$
$$alpha=\frac{e}{\epsilon_0}$$
$$n_0=10^6$$

We then compute the sources ($$S_e,S_i,S_\epsilon$$) 
for each of the equations by substituting the exact solution.

All boundary conditions in this case are Dirichlet type with the 
boundary values directly computed from the exact solution.

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
using `inputs_ib_x` file.

To calculate L2 norm error and also get the solution, use 
python script: `python get_L2norm_error.py plt00001 0` (use 0,1,2 for 
x,y or z direction for the second argument here).

<img src="https://github.com/user-attachments/assets/44ba5010-2e21-4d99-9696-3a38c5b8ffdc" width="600">

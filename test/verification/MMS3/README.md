# Method-of-Manufactured solutions (MMS) test 3

This case tests the time accuracy of our schemes. 
So, solves are just time-dependent ordinary differential equations 
at every computational cell.

$$\frac{dn_e}{dt} = k_i n_e$$
$$\frac{dn_i}{dt} = k_i n_e$$
$$\frac{dE_\epsilon}{dt} = k_\epsilon E_\epsilon$$

$$k_i=-5.0 \quad k_\epsilon=-10.0$$

Initially,

$$n_e=n_i=E_\epsilon=n_0 \left(1+\sin\left(\frac{\pi x}{L}\right)\right)$$

Exact solution after time t is 

$$n_e=n_i=n_0 \left(1+\sin\left(\frac{\pi x}{L}\right)\right) \exp(-5.0t)$$
$$E_\epsilon=n_0 \left(1+\sin\left(\frac{\pi x}{L}\right)\right) \exp(-10.0t)$$

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

Use the `run1.sh`, `run2.sh`, or `run4.sh` to run this case
with 1,2 and 4 timestep correctors. The implicit RK scheme in the code
will give 1st order accuracy with 1 step while 2nd order accuracy with 2 and 4 steps.

Use the python script get_L2norm_error.py as `python get_L2norm_error.py <pltfile> <direction> <finaltime> <dt>`
to get solution error in time and plots of the solution.


<img src="https://github.com/user-attachments/assets/593b1db5-131b-49fa-882b-ab256f69c1dd" width="600">


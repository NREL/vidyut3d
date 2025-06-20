# Linear advection test

This case tests linear advection with the 5th order WENO scheme and 
verifies its accuracy. 

### Build instructions

make sure $AMREX_HOME is set to your clone of amrex as
`$ export AMREX_HOME=/path/to/amrex`

If you are copying this case folder elsewhere then
make sure $VIDYUT_DIR is set to your clone of vidyut as
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

This case can be run in 2d as well as 1d. The 1d case is to verify 
5th order convergence with the WENO5 scheme. The default build is 2d.
If you want to build in 1d do `$ make -j COMP=gnu USE_MPI=TRUE DIM=1`

Run with inputs2d to see how advection of a Gaussian feature is 
captured with AMR. 
Using WENO5 just do `$ mpirun -n 1 ./*.ex inputs2d`
Using 1st order upwind do `$ mpirun -n 1 ./*.ex inputs2d vidyut.hyp_order=1`

<img src="https://github.com/user-attachments/assets/88434d4f-0ae5-4405-9a46-a13d8fb37f73" width="500">

Use the `plasjob` script and then the `postprocess.sh` script to run 
multiple 1d cases and verify order of convergence. The `pic_adv_weno.gp` 
script will plot convergence rates using gnuplot, do `gnuplot pic_adv_weno.gp`

<img src="https://github.com/user-attachments/assets/c3f27f30-edce-4bd3-b560-953f0b4e2e7f" width="400">




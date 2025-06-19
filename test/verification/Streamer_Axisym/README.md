# Axisymmetric streamer propagation

This case simulates air streamer propagation and is a 
replication of case 1 in the paper by
`Bagheri, Behnaz, et al. "Comparison of six simulation codes for positive streamers in air, 
Plasma Sources Science and Technology 27.9 (2018): 095002` 
see https://iopscience.iop.org/article/10.1088/1361-6595/aad768/meta

A seed of positive ions are initialized in a Gaussian kernel between 
two electrode boundaries separated by 1.25 cm with a powered electrode (top) 
at 18.75 kV and the bottom electrode grounded. 
This case is solved in the axisymmetric mode.
Local field approximation is used and Electron energy solves are turned off.

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

### Run instructions

Run with inputs2d that has a total of 3 AMR levels
`$ mpirun -n 64 ./*.ex inputs2d`

Run with inputs2d_amr_opt that has a coarser base grid 
but with a total 5 AMR levels
`$ mpirun -n 64 ./*.ex inputs2d_amr_opt`

<img src="https://github.com/user-attachments/assets/d0e73fa4-6d5c-44a2-a3db-483edcc0684a" width=600>



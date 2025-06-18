# Parallel plate He capacitive discharge

This is a 1d case that simulates He capacitive discharge

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

Run with inputs1d 
`$ mpirun -n 16 ./*.ex inputs1d`

 Run with inputs1d_ib to test
 cell masking
`$ mpirun -n 16 ./*.ex inputs1d_ib`

Use the track_residual.py script to see if 
quantities are approaching steady state.
`python track_residual.py "plt?????" "HEp"

Use the make_lineplot.py script to extract
quantities along an axial line: 
`python make_lineplot.py plt00010`

Use the compare.gp gnuplot script to get a plot 
of ion density compared against literature data
`gnuplot compare.gp`

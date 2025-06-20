# Parallel plate He capacitive discharge

This is a 1d case that simulates He capacitive discharge between two parallel plates. 
This case has been used as benchmark for PIC and fluid simulations. see

Turner, M. M., Derzsi, A., Donko, Z., Eremin, D., Kelly, S. J., Lafleur, T., & Mussenbrock, T. (2013). 
Simulation benchmarks for low-pressure plasmas: Capacitive discharges. Physics of Plasmas, 20(1).
https://pubs.aip.org/aip/pop/article/20/1/013507/1017414

We also compare our solution with other plasma fluid codes, 
for e.g. SOMAFOAM (https://www.sciencedirect.com/science/article/pii/S0010465521000229)
and mps1d (https://iopscience.iop.org/article/10.1088/1361-6463/ad7ecb/meta)

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

These runs will take about an hour to run 1000 cycles over
which the ion density profiles become steady.

Use the track_residual.py script to see if 
quantities are approaching steady state.
`python track_residual.py "plt?????" "HEp"

Use the make_lineplot.py script to extract
quantities along an axial line: 
`python make_lineplot.py plt00010`

Use the compare.gp gnuplot script to get a plot 
of ion density compared against literature data
`gnuplot compare.gp`
<img src="https://github.com/user-attachments/assets/29dda91f-09ce-4562-ac51-0ffcea583f42" width=500>



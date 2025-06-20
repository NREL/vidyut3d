# Ar discharge in a 3 electrode reactor

This case simulates a three electrode reactor similar to the one used
for thin film deposition. The two powered electrodes on the bottom 
boundary are "targets" from which metal is sputtered off by the plasma 
while the top boundary has a larger substrate electrode where deposition takes place.
In this case, we do not include any sputtering reactions. We include 
plasma reactions and secondary electron emission.

Reactor geometry is similar to the one in
Fébba, Davi M., et al. "Autonomous sputter synthesis of thin film nitrides with composition controlled by Bayesian optimization of optical plasma emission." APL Materials 11.7 (2023).

https://pubs.aip.org/aip/apm/article/11/7/071119/2903572

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

This is a 3D run with 2 levels of AMR and a cell count 
reaching approximately 3.4 million. It takes about
15 mins to do 1 RF cycle with 256 processors.
`mpirun -n 256 ./*.ex inputs3d`

Change the `amr.plot_int` variable to print outputs more often.

<img src="https://github.com/user-attachments/assets/03aeea0d-b38d-4d1c-9f46-7fcb668b5f6b" width=500>

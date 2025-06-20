# Vidyut3d: A fluid solver for simulating low-temperature plasmas

Vidyut (meaning lightning or electricity in Sanskrit) is a massively-parallel plasma-fluid solver 
for low-temperature plasmas (LTPs) that supports both local field (LFA) and local mean energy (LMEA) approximations, 
as well as complex gas and surface-phase chemistry. The solver supports 1D, 2D, 2D-axisymmetric and 3D domains, and uses adaptive mesh 
refinement capabilities to increase the grid resolution around complex structures (e.g. streamer heads and sheaths) 
while maintaining a tractable problem size. Vidyut specializes in simulating various types of gas-phase discharges, 
as well as plasma/surface interactions and surface chemistry (e.g. for plasma-mediated catalysis applications). 
The solver also supports hybrid CPU/GPU parallelization strategies, and has demonstrated excellent scaling on 
various high performance computing architectures for problem sizes consisting of O(100 M) control volumes.

## Models and Features

- multi-species two temperature model with self-consistent Poisson and electron energy equation solves
- LFA and LMEA models for solving the plasma-fluid equations with a drift-diffusion approximation
- Support for complex gas and surface-phase chemistry
- Second order semi-implicit scheme that handles drift and reactive source terms explicitly, and diffusive sources implicitly 
- Parallelization via OpenMPI/MPICH and GPU Acceleration with CUDA (NVidia) and HIP (AMD)
  
## Build requirements

* A clone of vidyut is required as `git clone https://github.com/NREL/vidyut3d`

* Vidyut heavily relies on AMReX, a performance portable adaptive Cartesian grid management and
  solver library. A clone of AMReX is required, preferably in the same location as the vidyut folder from
  previous step - `git clone https://github.com/AMReX-Codes/amrex` 

* Compilers: gcc and an MPI library (openMPI/MPICH) for CPU builds. cuda-11.0 is also required for GPU builds

## Code structure

* Like many AMReX based codes, the main numerical schemes are implemented in the `vidyut3d/Source` folder
* Any simulation with `vidyut` needs its own executable. It is built by linking case specific headers and c++ files with
  the source files. For an example case, please go to the `vidyut3d/test` folder and try out any of the
  `verification` or `model` cases. Each case here has an associated README file that is most readable through github on a web browser.
  

* Each example/run case must include
   - `ProbParm.H` : header with user controlled problem specific parameters, e.g., voltage at an electrode
   - `Prob.H` : header with functions that set initial conditions and problem specific parameters in ProbParm.H
   - `UserFunctions.H`: User specific functions for setting boundary conditions
   - `UserSources.H` : User specific sources for each equation
   - `Chemistry.H/cpp` : Plasma chemistry description including mobility,diffusion coefficients, and production Rates from chemical kinetics
   - `inputs` : includes run time parameters, e.g., final time, number of AMR levels
      
* Vidyut relies on AMReX's multi-level-multi-grid linear solvers and therefore uses all of the
  various boundary condition handles such as Dirichlet, Homogenous/inhomogenous Neumann and Robin. The setting of these conditions
  is done in `UserFunctions.H` within each case. More information about the implementation can be found here -
  https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#boundary-conditions
  
## Visualization instructions

* The outputs for a case are in the form of AMReX plotfiles
* These plot files can be open using AMReX grid reader in ParaView (see https://amrex-codes.github.io/amrex/docs_html/Visualization.html#paraview)
* Alternatively visit can be used. see https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html
* There are several python scripts in the `test` folder (e.g., /test/verification/He_RF_1d
/make_lineplot.py) that utilize yt (https://yt-project.org/). This provides a way to obtain data from plot files
  in the form of numpy arrays for analysis.
* AMReX's plotfile tools (https://github.com/AMReX-Codes/amrex/tree/development/Tools/Plotfile) are also great for getting average values or data along a line, among others.
* Amrvis (https://amrex-codes.github.io/amrex/docs_html/Visualization.html) is a great tool for quick visualizations 

## Acknowledgments

This work was authored by the National Renewable Energy Laboratory (NREL) under software record SWR-24-101, operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by funding from DOE Laboratory Directed Research and Development (LDRD) and DOE Basic Energy Sciences (BES). The research was performed using computational resources sponsored by the Department of Energy's Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory.

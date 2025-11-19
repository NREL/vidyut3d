#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <UnivConstants.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real Ar_frac = 0.3;
    amrex::Real CO2_frac = 0.6;
    amrex::Real H2_frac = 0.1;
    amrex::Real pres = 101325.0;
    amrex::Real temp = 300.0;
    amrex::Real seedion_frac = 1e-8;
    amrex::Real Vns = 800.0;
    amrex::Real pulse_freq = 1e4;
    amrex::Real pulse_width = 3e-8;
};
#endif

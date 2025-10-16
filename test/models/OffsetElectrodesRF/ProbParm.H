#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real pres = 133.0;
    amrex::Real temp = 300.0;
    amrex::Real len1 = 0.01;
    amrex::Real len2 = 0.01;
    amrex::Real rfvoltage = 120.0;
    amrex::Real dielthick = 0.0001;
    amrex::Real Kdiel = 5.0;
};
#endif

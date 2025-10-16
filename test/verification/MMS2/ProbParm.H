#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    int axdir = 0;
    amrex::Real solnfrac = 0.5;
    amrex::Real n0 = 1e6;
};
#endif

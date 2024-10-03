#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <UnivConstants.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real V=1.0;
    amrex::Real f=13.56e6;
    amrex::Real bg_gasden=1.0;
};
#endif

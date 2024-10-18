#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real r1 = 0.0; // source term for species 1
    amrex::Real r2 = 1.0; // source term for species 2
    int axdir=0;
};
#endif

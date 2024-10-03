#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    int one_dim=0;
    amrex::Real blobrad = 0.1; 
    amrex::Real velx = 1.0;
    amrex::Real vely = 0.0;
    amrex::Real velz = 0.0;
};
#endif

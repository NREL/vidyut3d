#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    int nseeds=1;
    amrex::Real* loc_x=nullptr;
    amrex::Real* loc_y=nullptr;
    amrex::Real* loc_z=nullptr;

    amrex::Real rad=0.1;
    amrex::Real ampl=0.1;
    amrex::Real bg=1.0;

    amrex::Real V1=1.0;
    amrex::Real V2=0.0;

    int twod=1;
};
#endif

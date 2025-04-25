#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <UnivConstants.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real init_E=1e15;
    amrex::Real init_I=1e15;
    amrex::Real init_etemp=0.5*eV;
    amrex::Real V1=1.0;
    amrex::Real V2=0.0;

    amrex::Real ns=1.0;
    amrex::Real phis=1.0;
};
#endif

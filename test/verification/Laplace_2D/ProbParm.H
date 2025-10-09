#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real inner_radius=0.1;
    amrex::Real outer_radius=0.2;
    amrex::Real inner_voltage=10.0;
    amrex::Real outer_voltage=20.0;
    int inner_dirichlet = 1;
    int outer_dirichlet = 1;
};
#endif

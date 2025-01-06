#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real left_voltage=10.0;
    amrex::Real right_voltage=20.0;
    int left_dirc=1;
    int right_dirc=1;
    int left_neumann=0;
    int right_neumann=0;
};
#endif

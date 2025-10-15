#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real pres = 133.0;
    amrex::Real temp = 300.0;
    amrex::Real y1 = 0.015;
    amrex::Real l1 = 0.01;
    amrex::Real y2 = 0.035;
    amrex::Real l2 = 0.01;
    amrex::Real rfv1 = 120.0;
    amrex::Real rfv2 = 0.0;
    amrex::Real freq1 = 13.56e6;
    amrex::Real freq2 = 0.0;
    amrex::Real dielthick = 0.0001;
    amrex::Real Kdiel = 5.0;
};
#endif

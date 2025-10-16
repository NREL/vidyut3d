#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <UnivConstants.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real Efield_arc = 1e3;
    amrex::Real wallpot = 0.0;
    amrex::Real bulkplasmatemp = 10000;
    amrex::Real bulkelectemp = 10000;
    amrex::Real walltemp = 3000;
    amrex::Real Arp_arc = 1e22;
    amrex::Real Hp_arc = 1e21;
    amrex::Real H_arc = 1e21;
    amrex::Real spec_arc = 1e15;
    amrex::Real E_arc = 1e21;
    amrex::Real gama = 0.02;
    amrex::Real h2frac = 0.05;
};
#endif

#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real pres = 133.0;
    amrex::Real temp = 300.0;
    amrex::Real subs_minx = 1e5;
    amrex::Real subs_miny = 1e5;
    amrex::Real subs_maxx = 1e5;
    amrex::Real subs_maxy = 1e5;
    amrex::Real trode1_minx = 1e5;
    amrex::Real trode1_miny = 1e5;
    amrex::Real trode1_maxx = 1e5;
    amrex::Real trode1_maxy = 1e5;
    amrex::Real trode2_minx = 1e5;
    amrex::Real trode2_miny = 1e5;
    amrex::Real trode2_maxx = 1e5;
    amrex::Real trode2_maxy = 1e5;
    amrex::Real v1 = 300;
    amrex::Real v2 = 300;
    amrex::Real ph1 = 0.0;
    amrex::Real ph2 = 0.0;
    amrex::Real freq1 = 13.56e6;
    amrex::Real freq2 = 13.56e6;
    amrex::Real ionden_init = 1e14;

    int circle_electrodes = 0;
    amrex::Real trode1_rad = 0.025;
    amrex::Real trode2_rad = 0.025;
    amrex::Real subs_rad = 0.05;
    amrex::Real trode1_cent = -0.025;
    amrex::Real trode2_cent = 0.025;

    // secondary emission coefficients
    amrex::Real gama_trode1 = 0.05;
    amrex::Real gama_trode2 = 0.05;
    amrex::Real gama_subs = 0.05;
};
#endif

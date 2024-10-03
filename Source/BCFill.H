#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <Chemistry.H>

using namespace amrex;

struct AmrCoreFill
{
    AMREX_GPU_DEVICE
    void operator()(
        const IntVect& iv,
        Array4<Real> const& data,
        const int dcomp,
        const int numcomp,
        GeometryData const& geom,
        const Real time,
        const BCRec* bcr,
        const int bcomp,
        const int orig_comp) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
        // NOTE: we dont need ext_dir as we are doing 
        // linear solves for all equations
        // we are spoofing fillpatch to foextrap
        // we may lose order of accuracy at boundaries for
        // hyperbolic flux, oh well!

        /*const int* domlo = geom.Domain().loVect();
          const int* domhi = geom.Domain().hiVect();
          const amrex::Real* prob_lo = geom.ProbLo();
          const amrex::Real* prob_hi = geom.ProbHi();
          const amrex::Real* dx = geom.CellSize();
          const amrex::Real x[AMREX_SPACEDIM] = {
          AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0], 
          prob_lo[1] + (iv[1] + 0.5) * dx[1], prob_lo[2] + (iv[2] + 0.5) * dx[2])};

          const int* bc = bcr->data();

          amrex::Real s_int[NVAR] = {0.0};
          amrex::Real s_ext[NVAR] = {0.0};

        // xlo and xhi
        int idir = 0;
        if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir]))
        {
        amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
        for (int n = 0; n < NVAR; n++)
        {
        s_int[n] = data(loc, n);
        }
        externalbc(x, s_int, s_ext, idir, +1, time, geom);
        for (int n = 0; n < NVAR; n++)
        {
        data(iv, n) = s_ext[n];
        }
        }*/
    }
};

namespace {
    static AmrCoreFill amrcore_fill_func;
}

#endif

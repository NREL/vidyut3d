#ifndef _REACTIONS_H_
#define _REACTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>
#include <VarDefines.H>

using namespace amrex;
namespace user_sources
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void add_user_react_sources(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& rxnsrc,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time,
                ProbParm const& prob_parm,
                amrex::Real Tg,amrex::Real Pg)
        {
            rxnsrc(i,j,k,S1_ID)+=prob_parm.r1;
            rxnsrc(i,j,k,S2_ID)+=prob_parm.r2;
        }

    AMREX_GPU_DEVICE AMREX_INLINE
        void add_user_potential_sources(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& source,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time,
                ProbParm const& prob_parm)
        {
          //add any user specific sources
        }
}
#endif

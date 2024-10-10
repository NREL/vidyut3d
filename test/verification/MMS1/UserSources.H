#ifndef _REACTIONS_H_
#define _REACTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace user_sources
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void  add_user_react_sources(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& rxnsrc,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time,
                              ProbParm const& prob_parm,
                              amrex::Real Tg,amrex::Real Pg)
    {
        int bigdir=0;
        amrex::Real src=0.0;
        amrex::Real xlen=prob_hi[0]-prob_lo[0];
#if AMREX_SPACEDIM > 1
        amrex::Real ylen=prob_hi[1]-prob_lo[1];
#if AMREX_SPACEDIM == 3
        amrex::Real zlen=prob_hi[2]-prob_lo[2];
#endif
#endif

        amrex::Real maxlen=xlen;
#if AMREX_SPACEDIM > 1
        maxlen = std::max(maxlen,ylen);
#if AMREX_SPACEDIM == 3
        maxlen = std::max(maxlen,zlen);
#endif
#endif

        if(maxlen==xlen) bigdir=0;
#if AMREX_SPACEDIM > 1
        if(maxlen==ylen) bigdir=1;
#if AMREX_SPACEDIM == 3
        if(maxlen==zlen) bigdir=2;
#endif
#endif

        IntVect iv{AMREX_D_DECL(i, j, k)};
        amrex::Real x=prob_lo[bigdir]+(iv[bigdir]+0.5)*dx[bigdir];
        rxnsrc(i,j,k,AR_ID) += 5.0*std::pow(x,4.0)/3.0-x/6.0-2.0;
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
        source(i,j,k)+=phi(i,j,k,AR_ID);
    }
}
#endif

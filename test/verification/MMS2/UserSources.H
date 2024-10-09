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
        amrex::Real n0=prob_parm.n0;
        amrex::Real alpha=ECHARGE/EPS0;
        amrex::Real Ee_exact=(x*x)/alpha+n0;
        amrex::Real ne_exact=(x*x)/alpha+n0;
        amrex::Real Te_exact=(2.0/3.0)*(Ee_exact/ne_exact)/K_B;
        amrex::Real Ta=2.0/K_B;
        amrex::Real Efield_exact=-1/24.0*(4*x*x*x-1.0);
        amrex::Real mue=-1.0;
        amrex::Real De=1.0;
        amrex::Real dnedx=2.0*x/alpha;
        amrex::Real gama_e=mue*ne_exact*Efield_exact-De*dnedx;

        amrex::Real izrate=5.0*std::exp(-Ta/Te_exact)*ne_exact*phi(i,j,k,HE_ID);
        amrex::Real dGamaedx=(1.0/alpha)*(5.0*std::pow(x,4.0)/6.0-x/12.0-2.0)+n0*x*x/2.0;
        amrex::Real dGamaidx=(1.0/(2.0*alpha))*(-5.0*std::pow(x,4.0)/12.0+x/24.0-1.0)-n0*x*x/4.0;
        amrex::Real dGamaEedx=(5.0/3.0)*dGamaedx;
        
        rxnsrc(i,j,k,E_ID) += (dGamaedx-izrate);
        rxnsrc(i,j,k,HEp_ID) += (dGamaidx-izrate);
        
        amrex::Real jheat=-ECHARGE*gama_e*Efield_exact;
        amrex::Real inel_loss=4.0*izrate;
        amrex::Real nu=10000.0;
        amrex::Real el_loss= 3.0/2.0 * K_B * ne_exact * (Te_exact-Tg) * nu * (2.0*ME/(4.0*M_AMU));
        rxnsrc(i,j,k,EEN_ID) += (dGamaEedx-jheat+inel_loss+el_loss);
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
        //source(i,j,k)+=phi(i,j,k,AR_ID);
    }
}
#endif

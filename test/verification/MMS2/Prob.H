#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

AMREX_INLINE
void amrex_probinit (ProbParm & h_prob_parm, ProbParm & d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    pp.query("axial_dir", h_prob_parm.axdir );
    pp.query("solnfrac", h_prob_parm.solnfrac );

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif

}

AMREX_GPU_DEVICE
AMREX_INLINE
void initdomaindata(Box const& bx, Array4<Real> const& phi, 
        GeometryData const& geomdata,
        ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();
    
#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                for(int c=0;c<NVAR;c++)
                {
                   phi(i,j,k,c)=1e-12;
                }
                Real alpha=ECHARGE/EPS0;
                IntVect iv(AMREX_D_DECL(i,j,k));
                Real frac=prob_parm->solnfrac;
                Real lenx=prob_hi[prob_parm->axdir]-prob_lo[prob_parm->axdir];
                Real x = prob_lo[prob_parm->axdir] 
                + (0.5+iv[prob_parm->axdir]) * dx[prob_parm->axdir];

                phi(i,j,k,CMASK_ID)=1.0;
                if(x*(x-1.0)>0.0)
                {
                   phi(i,j,k,CMASK_ID)=0.0;
                }

                //little twiddle so the electron density is not 0
                phi(i,j,k,E_ID)    = frac*(x*x*x/alpha+prob_parm->n0)+1e-5*prob_parm->n0;
                phi(i,j,k,HE_ID)   = 1.0; 
                phi(i,j,k,HEp_ID)  = frac*(x*x*x/(2.0*alpha)+prob_parm->n0);
                
                phi(i,j,k,POT_ID) = frac*(1.0/24.0)*(x*x*x*x*x-x);
                phi(i,j,k,EFX_ID) = 0.0;
                phi(i,j,k,EFY_ID) = 0.0;
                phi(i,j,k,EFZ_ID) = 0.0;
                phi(i,j,k,EEN_ID) = frac*(x*x*x/alpha+prob_parm->n0);
                phi(i,j,k,ETEMP_ID) = 1.0;
            }
        }
    }
}
#endif

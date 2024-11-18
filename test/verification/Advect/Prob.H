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
    pp.query("one_dim",h_prob_parm.one_dim);
    pp.query("blobrad", h_prob_parm.blobrad );
    pp.query("velx", h_prob_parm.velx );
    pp.query("vely", h_prob_parm.vely );
    pp.query("velz", h_prob_parm.velz );

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

    amrex::Real midx=0.5*(prob_lo[0]+prob_hi[0]);
    amrex::Real midy=0.5*(prob_lo[1]+prob_hi[1]);
    
#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real z = prob_lo[2] + (0.5+k) * dx[2];
            Real y = prob_lo[1] + (0.5+j) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real x = prob_lo[0] + (0.5+i) * dx[0];
                
                if(prob_parm->one_dim)
                {
                    //using same conditions as
                    //Ghosh, Debojyoti, and James D. Baeder. "Weighted non-linear 
                    //compact schemes for the direct numerical simulation of 
                    //compressible, turbulent flows." 
                    //Journal of Scientific Computing 61 (2014): 61-89.
                    phi(i,j,k,AR_ID)=sin(PI*x - sin(PI*x)/PI);
                }
                else
                {
                    phi(i,j,k,AR_ID)  = exp(-0.25*((x-midx)*(x-midx)+(y-midy)*(y-midy))
                                            /prob_parm->blobrad/prob_parm->blobrad);
                }

                phi(i,j,k,E_ID) = 0.0;
                phi(i,j,k,POT_ID) = 0.0;
                phi(i,j,k,EFX_ID) = 0.0;
                phi(i,j,k,EFY_ID) = 0.0;
                phi(i,j,k,EFZ_ID) = 0.0;
                phi(i,j,k,EEN_ID) = 0.0;
                phi(i,j,k,ETEMP_ID) = 0.5*eV;
            }
        }
    }
}
#endif

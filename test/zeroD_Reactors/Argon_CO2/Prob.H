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
    pp.query("init_E", h_prob_parm.init_E );
    pp.query("init_I", h_prob_parm.init_I );
    pp.query("init_etemp", h_prob_parm.init_etemp );
    pp.query("V1", h_prob_parm.V1 );
    pp.query("V2", h_prob_parm.V2 );
    pp.query("nscale",h_prob_parm.ns);
    pp.query("potscale",h_prob_parm.phis);

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
    
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real z = prob_lo[2] + (0.5+k) * dx[2];
                Real y = prob_lo[1] + (0.5+j) * dx[1];
                Real x = prob_lo[0] + (0.5+i) * dx[0];

                // No masking of reactor cells
                phi(i,j,k,CMASK_ID)=1.0;

                phi(i,j,k,0)  = 1e9;         // electron num dens.
                phi(i,j,k,3)  = 0.9*2.446e25;    // Ar num dens.
                phi(i,j,k,4)  = 0.1*2.446e25;    // CO2 num dens.
                // phi(i,j,k,4)  = 2.446e25;    // CO2 num dens.
                phi(i,j,k,ETEMP_ID) = 300;
                phi(i,j,k,EEN_ID) = 1.5*phi(i,j,k,0)*K_B*phi(i,j,k,ETEMP_ID);
                phi(i,j,k,EFX_ID) = 0.0;
                phi(i,j,k,EFY_ID) = 0.0;
                phi(i,j,k,EFZ_ID) = 0.0;
                phi(i,j,k,POT_ID) = 0.0;
            }
        }
    }
}
#endif

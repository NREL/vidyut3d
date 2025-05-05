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
    pp.query("left_voltage", h_prob_parm.left_voltage );
    pp.query("right_voltage", h_prob_parm.right_voltage );
    pp.query("left_dirc", h_prob_parm.left_dirc );
    pp.query("right_dirc", h_prob_parm.right_dirc );
    pp.query("left_neumann", h_prob_parm.left_neumann );
    pp.query("right_neumann", h_prob_parm.right_neumann );

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
                for(int sp=0;sp<NVAR;sp++)
                {
                   phi(i,j,k,sp)=0.0;
                }
                phi(i,j,k,CMASK_ID)=1;
                
                phi(i,j,k,E_ID)  = 1e12;
                phi(i,j,k,S1_ID)  = 1e12;
                phi(i,j,k,S2_ID)  = 1e12;
                phi(i,j,k,ETEMP_ID) = 0.5*eV;
            }
        }
    }
}
#endif

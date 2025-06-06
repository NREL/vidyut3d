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
    pp.query("y1", h_prob_parm.y1 );
    pp.query("l1", h_prob_parm.l1 );
    pp.query("y2", h_prob_parm.y2 );
    pp.query("l2", h_prob_parm.l2 );
    pp.query("rfv1", h_prob_parm.rfv1 );
    pp.query("rfv2", h_prob_parm.rfv2 );
    pp.query("dielectic_thickness", h_prob_parm.dielthick );
    pp.query("dielectic_constant", h_prob_parm.Kdiel );
    
    amrex::ParmParse pp_solver("vidyut");
    pp_solver.query("gas_pressure",h_prob_parm.pres); 
    pp_solver.query("gas_temperature",h_prob_parm.temp); 

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
                
                for(int c=0;c<NVAR;c++)
                {
                   phi(i,j,k,c)=0.0;
                }
                phi(i,j,k,CMASK_ID)=1.0;

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,sp)=1e12;
                }
               
                phi(i,j,k,HE_ID)   = prob_parm->pres/K_B/prob_parm->temp;  
                phi(i,j,k,HEp_ID)  = 1e14;
                phi(i,j,k,E_ID)    = 1e14;
                phi(i,j,k,ETEMP_ID) = 0.5*eV;
                phi(i,j,k,EEN_ID) = 1.5*phi(i,j,k,E_ID)*K_B*phi(i,j,k,ETEMP_ID);
            }
        }
    }
}
#endif

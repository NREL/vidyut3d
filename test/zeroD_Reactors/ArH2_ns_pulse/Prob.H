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
    pp.query("Efield_arc",h_prob_parm.Efield_arc);
    pp.query("wallpot",h_prob_parm.wallpot);
    pp.query("bulkplasmatemp",h_prob_parm.bulkplasmatemp);
    pp.query("bulkelectemp",h_prob_parm.bulkelectemp);
    pp.query("walltemp",h_prob_parm.walltemp);
    pp.query("Arp_arc",h_prob_parm.Arp_arc);
    pp.query("Hp_arc",h_prob_parm.Hp_arc);
    pp.query("H_arc",h_prob_parm.H_arc);
    pp.query("spec_arc",h_prob_parm.spec_arc);
    pp.query("gama",h_prob_parm.gama);
    pp.query("h2frac",h_prob_parm.h2frac);

    //skips electrons in the next loop anyway
    //so setting E_ID is ok here
    amrex::Real arcspecies[NUM_SPECIES]={h_prob_parm.spec_arc};
    arcspecies[ARp_ID]=h_prob_parm.Arp_arc;
    arcspecies[Hp_ID]=h_prob_parm.Hp_arc;
    arcspecies[H_ID]=h_prob_parm.H_arc;

    int chrgspec[NUM_SPECIES]={0};
    CKCHRG(chrgspec);
    
    amrex::Real net_ion_charge=0.0;
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        if(sp!=E_ID)
        {
            net_ion_charge += arcspecies[sp]*chrgspec[sp];
        }
    }
    h_prob_parm.E_arc=net_ion_charge;

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

    amrex::Real L=prob_hi[0]-prob_lo[0];
    amrex::Real deltaTemp=prob_parm->bulkplasmatemp-prob_parm->walltemp;


    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {

            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real z = prob_lo[2] + (0.5+k) * dx[2];
                Real y = prob_lo[1] + (0.5+j) * dx[1];
                Real x = prob_lo[0] + (0.5+i) * dx[0];

                int chrgspec[NUM_SPECIES]={0};
                CKCHRG(chrgspec);

                Real gastemp=prob_parm->bulkplasmatemp-deltaTemp/L*(x-prob_lo[0]);

                for(int c=0;c<NVAR;c++)
                {
                    phi(i,j,k,c)=0.0;
                }
                phi(i,j,k,CMASK_ID)=1.0;
                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,sp)=prob_parm->spec_arc;
                }
                
                phi(i,j,k,AR_ID)=(1.0-prob_parm->h2frac)*P_NTP/K_B/gastemp;
                phi(i,j,k,H2_ID)=prob_parm->h2frac*P_NTP/K_B/gastemp;

                //phi(i,j,k,ARp_ID)=prob_parm->Arp_arc;
                //phi(i,j,k,Hp_ID)=prob_parm->Hp_arc;
                //phi(i,j,k,H_ID)=prob_parm->H_arc;

                amrex::Real net_ion_charge=0.0;
                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    if(sp!=E_ID)
                    {
                        net_ion_charge += phi(i,j,k,sp)*chrgspec[sp];
                    }
                }

                phi(i,j,k,E_ID)    = net_ion_charge;
                phi(i,j,k,ETEMP_ID) = prob_parm->bulkelectemp;
                phi(i,j,k,EEN_ID) = 1.5*phi(i,j,k,E_ID)*K_B*phi(i,j,k,ETEMP_ID);
            }
        }
    }
}
#endif

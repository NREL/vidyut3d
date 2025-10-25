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
void amrex_probinit(ProbParm& h_prob_parm, ProbParm& d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    pp.query("Ar_frac", h_prob_parm.Ar_frac);
    pp.query("CO2_frac", h_prob_parm.CO2_frac);
    pp.query("seedion_frac", h_prob_parm.seedion_frac);
    h_prob_parm.H2_frac=1.0-h_prob_parm.Ar_frac-h_prob_parm.CO2_frac;

    amrex::ParmParse ppv("vidyut");
    pp.query("gas_pressure", h_prob_parm.pres);
    pp.query("gas_temperature", h_prob_parm.temp);


#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif
}

AMREX_GPU_DEVICE
AMREX_INLINE
void initdomaindata(
    Box const& bx,
    Array4<Real> const& phi,
    GeometryData const& geomdata,
    ProbParm* prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx = geomdata.CellSize();
    
    int chrgspec[NUM_SPECIES] = {0};
    CKCHRG(chrgspec);

    for (int k = lo.z; k <= hi.z; ++k)
    {
        for (int j = lo.y; j <= hi.y; ++j)
        {

            for (int i = lo.x; i <= hi.x; ++i)
            {
                Real z = prob_lo[2] + (0.5 + k) * dx[2];
                Real y = prob_lo[1] + (0.5 + j) * dx[1];
                Real x = prob_lo[0] + (0.5 + i) * dx[0];

                Real numden=prob_parm->pres/K_B/prob_parm->temp;

                // No masking of reactor cells
                phi(i, j, k, CMASK_ID) = 1.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                {
                    phi(i, j, k, sp) = prob_parm->seedion_frac*0.001*numden;
                }

                phi(i, j, k, ARp_ID)  = prob_parm->seedion_frac*numden/3.0;
                phi(i, j, k, H2p_ID)  = prob_parm->seedion_frac*numden/3.0;
                phi(i, j, k, CO2p_ID) = prob_parm->seedion_frac*numden/3.0;

                amrex::Real net_ion_charge = 0.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                {
                    if (sp != E_ID)
                    {
                        net_ion_charge += phi(i, j, k, sp) * chrgspec[sp];
                    }
                }
                phi(i, j, k, E_ID) = net_ion_charge;

                phi(i, j, k, AR_ID) = prob_parm->Ar_frac*numden;
                phi(i, j, k, H2_ID) = prob_parm->H2_frac*numden; 
                phi(i, j, k, CO2_ID) = prob_parm->CO2_frac*numden;

                phi(i, j, k, ETEMP_ID) = 0.5*eV;
                phi(i, j, k, EEN_ID) =
                1.5 * phi(i, j, k, 0) * K_B * phi(i, j, k, ETEMP_ID);
                phi(i, j, k, EFX_ID) = 0.0;
                phi(i, j, k, EFY_ID) = 0.0;
                phi(i, j, k, EFZ_ID) = 0.0;
                phi(i, j, k, POT_ID) = 0.0;
            }
        }
    }
}
#endif

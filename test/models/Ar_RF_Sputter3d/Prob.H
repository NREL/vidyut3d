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
    pp.query("v1", h_prob_parm.v1);
    pp.query("v2", h_prob_parm.v2);
    pp.query("ph1", h_prob_parm.ph1);
    pp.query("ph2", h_prob_parm.ph2);
    pp.query("freq1", h_prob_parm.freq1);
    pp.query("freq2", h_prob_parm.freq2);
    pp.query("ionden_init", h_prob_parm.ionden_init);
    pp.query("gama_trode1", h_prob_parm.gama_trode1);
    pp.query("gama_trode2", h_prob_parm.gama_trode2);
    pp.query("gama_subs", h_prob_parm.gama_subs);

    pp.query("trode1_minx", h_prob_parm.trode1_minx);
    pp.query("trode1_miny", h_prob_parm.trode1_miny);
    pp.query("trode1_maxx", h_prob_parm.trode1_maxx);
    pp.query("trode1_maxy", h_prob_parm.trode1_maxy);

    pp.query("trode2_minx", h_prob_parm.trode2_minx);
    pp.query("trode2_miny", h_prob_parm.trode2_miny);
    pp.query("trode2_maxx", h_prob_parm.trode2_maxx);
    pp.query("trode2_maxy", h_prob_parm.trode2_maxy);

    pp.query("subs_minx", h_prob_parm.subs_minx);
    pp.query("subs_miny", h_prob_parm.subs_miny);
    pp.query("subs_maxx", h_prob_parm.subs_maxx);
    pp.query("subs_maxy", h_prob_parm.subs_maxy);

    amrex::ParmParse pp_solver("vidyut");
    pp_solver.query("gas_pressure", h_prob_parm.pres);
    pp_solver.query("gas_temperature", h_prob_parm.temp);

    pp.query("trode1_rad", h_prob_parm.trode1_rad);
    pp.query("trode2_rad", h_prob_parm.trode2_rad);
    pp.query("subs_rad", h_prob_parm.subs_rad);
    pp.query("trode1_cent", h_prob_parm.trode1_cent);
    pp.query("trode2_cent", h_prob_parm.trode2_cent);
    pp.query("circle_electrodes", h_prob_parm.circle_electrodes);

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

    for (int k = lo.z; k <= hi.z; ++k)
    {
        for (int j = lo.y; j <= hi.y; ++j)
        {

            for (int i = lo.x; i <= hi.x; ++i)
            {
                Real z = prob_lo[2] + (0.5 + k) * dx[2];
                Real y = prob_lo[1] + (0.5 + j) * dx[1];
                Real x = prob_lo[0] + (0.5 + i) * dx[0];

                for (int n = 0; n < NVAR; n++)
                {
                    phi(i, j, k, n) = 0.0;
                }
                phi(i, j, k, CMASK_ID) = 1.0;

                for (int sp = 0; sp < NUM_SPECIES; sp++)
                {
                    phi(i, j, k, sp) = 1e12;
                }

                phi(i, j, k, AR_ID) = prob_parm->pres / K_B / prob_parm->temp;
                phi(i, j, k, ARp_ID) = prob_parm->ionden_init;
                phi(i, j, k, E_ID) =
                    phi(i, j, k, ARp_ID) + phi(i, j, k, AR2p_ID);
                phi(i, j, k, ETEMP_ID) = 0.5 * eV;
                phi(i, j, k, EEN_ID) =
                    1.5 * phi(i, j, k, E_ID) * K_B * phi(i, j, k, ETEMP_ID);
            }
        }
    }
}
#endif

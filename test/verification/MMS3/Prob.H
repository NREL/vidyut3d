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
    pp.query("axial_dir", h_prob_parm.axdir);

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
                Real alpha = ECHARGE / EPS0;
                IntVect iv(AMREX_D_DECL(i, j, k));
                Real lenx =
                    prob_hi[prob_parm->axdir] - prob_lo[prob_parm->axdir];

                for (int c = 0; c < NVAR; c++)
                {
                    phi(i, j, k, c) = 0.0;
                }
                phi(i, j, k, CMASK_ID) = 1.0;

                Real x = prob_lo[prob_parm->axdir] +
                         (0.5 + iv[prob_parm->axdir]) * dx[prob_parm->axdir];

                // little twiddle so the electron density is not 0
                phi(i, j, k, E_ID) = prob_parm->n0 * (1 + sin(PI * x / lenx));
                phi(i, j, k, HE_ID) = 1.0;
                phi(i, j, k, HEp_ID) = prob_parm->n0 * (1 + sin(PI * x / lenx));
                phi(i, j, k, EEN_ID) = prob_parm->n0 * (1 + sin(PI * x / lenx));
                phi(i, j, k, ETEMP_ID) = 0.666;
            }
        }
    }
}
#endif

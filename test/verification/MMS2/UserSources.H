#ifndef _REACTIONS_H_
#define _REACTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace user_sources {
AMREX_GPU_DEVICE AMREX_INLINE void add_user_react_sources(
    int i,
    int j,
    int k,
    Array4<Real> const& phi,
    Array4<Real> const& rxnsrc,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    amrex::Real Tg,
    amrex::Real Pg)
{
    int bigdir = 0;
    amrex::Real src = 0.0;
    amrex::Real xlen = prob_hi[0] - prob_lo[0];
#if AMREX_SPACEDIM > 1
    amrex::Real ylen = prob_hi[1] - prob_lo[1];
#if AMREX_SPACEDIM == 3
    amrex::Real zlen = prob_hi[2] - prob_lo[2];
#endif
#endif

    amrex::Real maxlen = xlen;
#if AMREX_SPACEDIM > 1
    maxlen = std::max(maxlen, ylen);
#if AMREX_SPACEDIM == 3
    maxlen = std::max(maxlen, zlen);
#endif
#endif

    if (maxlen == xlen) bigdir = 0;
#if AMREX_SPACEDIM > 1
    if (maxlen == ylen) bigdir = 1;
#if AMREX_SPACEDIM == 3
    if (maxlen == zlen) bigdir = 2;
#endif
#endif

    IntVect iv{AMREX_D_DECL(i, j, k)};
    amrex::Real x = prob_lo[bigdir] + (iv[bigdir] + 0.5) * dx[bigdir];
    amrex::Real n0 = prob_parm.n0;
    amrex::Real alpha = ECHARGE / EPS0;
    amrex::Real Ee = (x * x * x) / alpha + n0;
    amrex::Real ne = (x * x * x) / alpha + n0;
    amrex::Real dnedx = 3.0 * x * x / alpha;
    amrex::Real d2nedx2 = 6.0 * x / alpha;
    amrex::Real ni = (x * x * x) / (2.0 * alpha) + n0;
    amrex::Real dnidx = 3.0 * x * x / (2.0 * alpha);
    amrex::Real d2nidx2 = 6.0 * x / (2.0 * alpha);

    amrex::Real E = (-1 / 40.0) * (5 * x * x * x * x - 1.0);
    amrex::Real dEdx = -x * x * x / 2.0;

    amrex::Real Te = (2.0 / 3.0) * (Ee / ne) / K_B;
    amrex::Real Ta = 2.0 / K_B;
    amrex::Real mue = -1.0;
    amrex::Real De = 1.0;
    amrex::Real mui = 0.5;
    amrex::Real Di = 0.5;

    amrex::Real gama_e = mue * ne * E - De * dnedx;
    amrex::Real dgamaedx = mue * dnedx * E + mue * ne * dEdx - De * d2nedx2;
    amrex::Real dgamaEedx = (5.0 / 3.0) * dgamaedx;

    amrex::Real gama_i = mui * ni * E - Di * dnidx;
    amrex::Real dgamaidx = mui * dnidx * E + mui * ni * dEdx - Di * d2nidx2;
    amrex::Real izrate = 5.0 * std::exp(-Ta / Te) * ne * phi(i, j, k, HE_ID);

    rxnsrc(i, j, k, E_ID) += (dgamaedx - izrate);
    rxnsrc(i, j, k, HEp_ID) += (dgamaidx - izrate);

    amrex::Real jheat = -ECHARGE * gama_e * E;
    amrex::Real inel_loss = 4.0 * izrate;
    amrex::Real nu = 10000.0;
    amrex::Real el_loss =
        3.0 / 2.0 * K_B * ne * (Te - Tg) * nu * (2.0 * ME / (4.0 * M_AMU));
    rxnsrc(i, j, k, EEN_ID) += (dgamaEedx - jheat + inel_loss + el_loss);
}

AMREX_GPU_DEVICE AMREX_INLINE void add_user_potential_sources(
    int i,
    int j,
    int k,
    Array4<Real> const& phi,
    Array4<Real> const& source,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm)
{
    // add any user specific sources
    // source(i,j,k)+=phi(i,j,k,AR_ID);
}
} // namespace user_sources
#endif

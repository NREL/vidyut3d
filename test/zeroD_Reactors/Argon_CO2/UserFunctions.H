#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace user_transport {
AMREX_GPU_DEVICE AMREX_INLINE void get_photoion_acoeff(
    int i,
    int j,
    int k,
    int sph_id,
    Array4<Real> const& sb_arr,
    Array4<Real> const& acoeff,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    Real Tg,
    Real Pg)
{}

AMREX_GPU_DEVICE AMREX_INLINE void get_photoion_rhs(
    int i,
    int j,
    int k,
    int sph_id,
    Array4<Real> const& sb_arr,
    Array4<Real> const& rhs_arr,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    Real Tg,
    Real Pg)
{}

AMREX_GPU_DEVICE AMREX_INLINE void bc_ib(
    IntVect face,
    int dir,
    int sgn,
    int solved_comp,
    int rhs_comp,
    Array4<Real> const& sb_arr,
    Array4<Real> const& acoeff,
    Array4<Real> const& rhs,
    GpuArray<int, AMREX_SPACEDIM> domlo,
    GpuArray<int, AMREX_SPACEDIM> domhi,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    Real Tg,
    Real Pg)
{}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int is_dielectric(
    int i,
    int j,
    int k,
    int dir,
    int sign,
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    const GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    const GpuArray<Real, AMREX_SPACEDIM> dx,
    Real time,
    ProbParm const& prob_parm)
{
    return (0);
}

AMREX_GPU_DEVICE AMREX_INLINE void potential_bc(
    int i,
    int j,
    int k,
    int dir,
    int sgn,
    Array4<Real> const& phi,
    Array4<Real> const& bc_arr,
    Array4<Real> const& robin_a,
    Array4<Real> const& robin_b,
    Array4<Real> const& robin_f,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    amrex::Real Tg,
    amrex::Real Pg,
    amrex::Real app_voltage)
{
    // Empty
}

AMREX_GPU_DEVICE AMREX_INLINE void photoionization_bc(
    int i,
    int j,
    int k,
    int dir,
    int sgn,
    Array4<Real> const& phi,
    Array4<Real> const& bc_arr,
    Array4<Real> const& robin_a,
    Array4<Real> const& robin_b,
    Array4<Real> const& robin_f,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    amrex::Real Tg,
    amrex::Real Pg)
{
    // Empty
}

AMREX_GPU_DEVICE AMREX_INLINE void species_bc(
    int i,
    int j,
    int k,
    int dir,
    int sgn,
    int spec_id,
    int bcspec_id,
    Array4<Real> const& phi,
    Array4<Real> const& bc_arr,
    Array4<Real> const& robin_a,
    Array4<Real> const& robin_b,
    Array4<Real> const& robin_f,
    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    amrex::Real Tg,
    amrex::Real Pg)

{
    // Empty
}

AMREX_GPU_DEVICE AMREX_INLINE amrex::Real compute_vel(
    IntVect iv,
    int dir,
    int specid,
    Array4<Real> const& phi,
    GpuArray<Real, AMREX_SPACEDIM> dx,
    const Real time,
    ProbParm const& prob_parm,
    amrex::Real Tg,
    amrex::Real Pg)
{
    return (0);
}
} // namespace user_transport
#endif

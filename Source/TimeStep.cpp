#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <UserFunctions.H>

// a wrapper for EstTimeStep
void Vidyut::ComputeDt(
    amrex::Real cur_time,
    amrex::Real dt_delay,
    amrex::Real dt_edrift,
    amrex::Real dt_ediff,
    amrex::Real dt_diel_relax)
{
    BL_PROFILE("vidyut::ComputeDt()");
    if (adaptive_dt && cur_time > dt_delay)
    {
        amrex::Real old_dt = dt[0];
        amrex::Real trans_min_dt = std::numeric_limits<Real>::max();
        amrex::Real adp_dt = std::numeric_limits<Real>::max();
        if (do_transport)
            trans_min_dt =
                std::min(dt_edrift * advective_cfl, dt_ediff * diffusive_cfl);
        adp_dt = std::min(trans_min_dt, dt_diel_relax * dielectric_cfl);
        if (adp_dt > dt_max) adp_dt = dt_max;
        if (adp_dt < dt_min) adp_dt = dt_min;
        dt[0] = std::min(old_dt * dt_stretch, adp_dt);
    } else
    {
        dt[0] = fixed_dt; // from input file
    }

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        dt[lev] = dt[lev - 1];
    }
}

// compute dt from CFL considerations
void Vidyut::find_time_scales(
    int lev,
    amrex::Real& dt_edrift,
    amrex::Real& dt_ediff,
    amrex::Real& dt_diel_relax)
{
    BL_PROFILE("Vidyut::find_time_scales()");
    dt_edrift = std::numeric_limits<Real>::max();
    dt_diel_relax = std::numeric_limits<Real>::max();
    dt_ediff = std::numeric_limits<Real>::max();

    const auto dx = geom[lev].CellSizeArray();
    MultiFab& S_new = phi_new[lev];

    amrex::Real captured_gastemp = gas_temperature;
    amrex::Real captured_gaspres = gas_pressure;
    int eidx = E_IDX;

    MultiFab edriftvel(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab ediff(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab mue_ne(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        auto state_arrays = S_new.arrays();
        auto edriftvel_arrays = edriftvel.arrays();
        auto ediff_arrays = ediff.arrays();
        auto mue_ne_arrays = mue_ne.arrays();
        auto prob_lo = geom[lev].ProbLoArray();

        amrex::ParallelFor(
            S_new, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                auto state_array = state_arrays[nbx];
                auto edriftvel_array = edriftvel_arrays[nbx];
                auto ediff_array = ediff_arrays[nbx];
                auto mue_ne_array = mue_ne_arrays[nbx];

                amrex::Real etemp = state_array(i, j, k, ETEMP_ID);
                amrex::Real Esum = 0.0;
                for (int dim = 0; dim < AMREX_SPACEDIM; dim++)
                    Esum += amrex::Math::powi<2>(
                        state_array(i, j, k, EFX_ID + dim));
                amrex::Real efield_mag = std::sqrt(Esum);

                amrex::Real ndens = 0.0;
                for (int sp = 0; sp < NUM_SPECIES; sp++)
                    ndens += state_array(i, j, k, sp);

                amrex::Real dcoeff =
                    specDiff(eidx, etemp, ndens, efield_mag, captured_gastemp);

                amrex::Real mu =
                    specMob(eidx, etemp, ndens, efield_mag, captured_gastemp);

                edriftvel_array(i, j, k) = amrex::Math::abs(mu) * efield_mag;
                ediff_array(i, j, k) = dcoeff;
                mue_ne_array(i, j, k) =
                    amrex::Math::abs(mu) * state_array(i, j, k, eidx);
            });
    }

    amrex::Real max_edriftvel = edriftvel.norm0(0, 0, true);
    amrex::Real max_ediff = ediff.norm0(0, 0, true);
    amrex::Real max_muene = mue_ne.norm0(0, 0, true);
    if (max_edriftvel > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            dt_edrift = std::min(dt_edrift, dx[i] / max_edriftvel);
        }
    }
    if (max_ediff > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            dt_ediff = std::min(
                dt_ediff, 0.5 * dx[i] * dx[i] / max_ediff / AMREX_SPACEDIM);
        }
    }
    if (max_muene > 0)
    {
        dt_diel_relax = EPS0 / ECHARGE / max_muene;
    }

    ParallelDescriptor::ReduceRealMin(dt_edrift);
    ParallelDescriptor::ReduceRealMin(dt_ediff);
    ParallelDescriptor::ReduceRealMin(dt_diel_relax);
}

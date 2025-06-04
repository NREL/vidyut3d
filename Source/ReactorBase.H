#ifndef REACTORBASE_H
#define REACTORBASE_H

#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>
#include <sundials/sundials_context.h>
#include "Factory.H"
// // #include "PelePhysics.H"
#include "Chemistry.H"
#include "ReactorUtils.H"

namespace pele::physics::reactions {

class ReactorBase : public Factory<ReactorBase>
{
public:
  static std::string base_identifier() { return "ReactorBase"; }

  virtual int init(int reactor_type, int ncells) = 0;

  virtual void close() = 0;

  virtual void print_final_stats(void* sundials_mem) = 0;

  virtual int react(
    const amrex::Box& box,
    amrex::Array4<amrex::Real> const& nspec_in,
    amrex::Array4<amrex::Real> const& Ue_in,
    amrex::Array4<amrex::Real> const& Te_in,
    amrex::Real Tgas_in,
    amrex::Array4<amrex::Real> const& EN_in,
    amrex::Real& dt_react,
    amrex::Real& time
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) = 0;

  virtual int surf_react(
    const amrex::Box& box,
    amrex::Array4<amrex::Real> const& nspec_in,
    amrex::Real Tgas_in,
    amrex::Real& dt_react,
    amrex::Real& time
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) = 0;

  virtual void flatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& nspec,
    amrex::Array4<const amrex::Real> const& Ue,
    amrex::Array4<const amrex::Real> const& EN,
    amrex::Real* EN_vect,
    amrex::Real* y_vect) = 0;

  virtual void unflatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<amrex::Real> const& nspec,
    amrex::Array4<amrex::Real> const& Ue,
    amrex::Array4<amrex::Real> const& Te,
    amrex::Real* y_vect,
    amrex::Real dt) = 0;

  virtual void surf_flatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& nspec,
    amrex::Real* y_vect) = 0;

  virtual void surf_unflatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<amrex::Real> const& nspec,
    amrex::Real* y_vect,
    amrex::Real dt) = 0;

  void set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals);

  ~ReactorBase() override = default;

protected:
  int verbose{0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 1> m_typ_vals = {0.0};
};
} // namespace pele::physics::reactions
#endif

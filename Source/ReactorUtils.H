#ifndef REACTORUTILS_H
#define REACTORUTILS_H

#include <utility>
#include <AMReX_Gpu.H>
#include "Chemistry.H"
// #include "PelePhysics.H"
#include "ReactorTypes.H"
#include "AMReX_Sundials.H"

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <cvode/cvode.h>

#include <nvector/nvector_serial.h>

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_SYCL
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele::physics::reactions::utils {
int check_flag(void* flagvalue, const char* funcname, int opt);

struct YCOrder
{
  // default: [ [Y_0,...Y_n]_{cell 0}, ..., [Y_0,...Y_n]_{cell m}]
  using ordering_type = YCOrder;
};

struct CYOrder
{
  // [ [Y_{cell 0},...Y_m]_{species 0}, ..., [Y_{cell 0},...Y_m]_{species n}]
  using ordering_type = CYOrder;
};

template <typename T>
struct is_ordering_type
{
  static constexpr bool value = std::is_base_of_v<typename T::ordering_type, T>;
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * (NUM_SPECIES + 1) + specie_idx;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * NUM_SPECIES + specie_idx;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return specie_idx * ncells + icell;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return vec_index<CYOrder>(specie_idx, icell, ncells);
}

template <typename OrderType>
void
set_sundials_solver_tols(
  // cppcheck-suppress constParameter
  ::sundials::Context& sunctx,
  void* sundials_mem,
  const int ncells,
  const amrex::Real relTol,
  const amrex::Real absTol,
  const amrex::GpuArray<amrex::Real, NUM_SPECIES + 1>& typ_vals,
  const std::string& solvername,
  const int verbose)
{
  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  const int neq_tot = (NUM_SPECIES + 1) * ncells;

#if defined(AMREX_USE_CUDA)
  N_Vector atol = N_VNewWithMemHelp_Cuda(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Cuda(atol);
#elif defined(AMREX_USE_HIP)
  N_Vector atol = N_VNewWithMemHelp_Hip(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Hip(atol);
#elif defined(AMREX_USE_SYCL)
  N_Vector atol = N_VNewWithMemHelp_Sycl(
    neq_tot, /*use_managed_mem=*/false,
    *amrex::sundials::The_SUNMemory_Helper(),
    &amrex::Gpu::Device::streamQueue(), sunctx);
  amrex::Real* ratol = N_VGetDeviceArrayPointer_Sycl(atol);
#else
  N_Vector atol = N_VNew_Serial(neq_tot, sunctx);
  amrex::Real* ratol = N_VGetArrayPointer(atol);
#endif

  if (typ_vals[0] > 0.0) {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances with TypVals rtol = " << relTol
                     << " atolfact = " << absTol << " in PelePhysics \n";
    }

    amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) {
      for (int k = 0; k < NUM_SPECIES + 1; k++) {
        ratol[vec_index<OrderType>(k, icell, ncells)] = typ_vals[k] * absTol;
      }
    });
  } else {
    // cppcheck-suppress knownConditionTrueFalse
    if ((verbose > 0) && (omp_thread == 0)) {
      amrex::Print() << " Setting " << solvername
                     << " tolerances rtol = " << relTol << " atol = " << absTol
                     << " in PelePhysics \n";
    }
    amrex::ParallelFor(
      neq_tot, [=] AMREX_GPU_DEVICE(int i) { ratol[i] = absTol; });
  }

  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  int flag;
  if (solvername == "cvode") {
    flag = CVodeSVtolerances(sundials_mem, relTol, atol);
  } else if ((solvername == "arkstep") || (solvername == "erkstep")) {
    flag = ARKodeSVtolerances(sundials_mem, relTol, atol);
  } else {
    amrex::Abort("setSundialsSolverTols not implemented for this solver type");
  }
  if (utils::check_flag(&flag, "SVtolerances", 1) != 0) {
    amrex::Abort("Problem in setSundialsSolverTols");
  }

  N_VDestroy(atol);
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_flatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  amrex::Array4<const amrex::Real> const& nspec,
  amrex::Array4<const amrex::Real> const& Ue,
  amrex::Array4<const amrex::Real> const& EN,
  amrex::Real* EN_vect,
  amrex::Real* y_vect)
{
  // amrex::Real mass_frac[NUM_SPECIES];
  // amrex::Real rho = 0.0;
  // amrex::Real rho_old = 0.0;
  // if (clean_init_massfrac) {
  //   for (int n = 0; n < NUM_SPECIES; n++) {
  //     rho_old += rhoY(i, j, k, n);
  //   }
  // }
  // for (int n = 0; n < NUM_SPECIES; n++) {
  //   const int idx2 = spec_index<OrderType>(n, icell, ncells);
  //   const amrex::Real rhoYn =
  //     (clean_init_massfrac)
  //       ? amrex::min<amrex::Real>(
  //           amrex::max<amrex::Real>(0.0, rhoY(i, j, k, n)), rho_old)
  //       : rhoY(i, j, k, n);
  //   src_vect[idx2] = frcExt(i, j, k, n);
  //   rho += rhoYn;
  //   mass_frac[n] = rhoYn;
  // }
  // const amrex::Real rho_inv = 1.0 / rho;
  // for (amrex::Real& Yn : mass_frac) {
  //   Yn *= rho_inv;
  // }

  for (int n = 0; n < NUM_SPECIES; n++) {
    const int idx = vec_index<OrderType>(n, icell, ncells);
    y_vect[idx] = nspec(i, j, k, n);
  }

  // Convert to electron energy
  y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] = Ue(i,j,k,0);
  EN_vect[icell] = EN(i, j, k, 0);

  // amrex::Real Enrg_loc = vect_energy[icell] * rho_inv;
  // auto eos = pele::physics::PhysicsType::eos();
  // if (reactor_type == ReactorTypes::e_reactor_type) {
  //   eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  // } else if (reactor_type == ReactorTypes::h_reactor_type) {
  //   eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  // } else {
  //   amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  // }
  // y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] = temp;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_unflatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  amrex::Array4<amrex::Real> const& nspec,
  amrex::Array4<amrex::Real> const& Ue,
  amrex::Array4<amrex::Real> const& Te,
  const amrex::Real* y_vect,
  amrex::Real dt)
{
  // Get the old Ys & rho as well as the clipped version
  // that went into CVODE
  // amrex::Real rhoY_oldclipped[NUM_SPECIES];
  // amrex::Real rho_old = 0.0;
  // amrex::Real rho_oldclipped = 0.0;
  // if (clean_init_massfrac) {
  //   for (int n = 0; n < NUM_SPECIES; n++) {
  //     rho_old += rhoY(i, j, k, n);
  //   }
  //   for (int n = 0; n < NUM_SPECIES; n++) {
  //     rhoY_oldclipped[n] = amrex::min<amrex::Real>(
  //       amrex::max<amrex::Real>(0.0, rhoY(i, j, k, n)), rho_old);
  //     rho_oldclipped += rhoY_oldclipped[n];
  //   }
  // }

  // amrex::Real mass_frac[NUM_SPECIES];
  // amrex::Real rho = 0.0;
  // for (int n = 0; n < NUM_SPECIES; n++) {
  //   // With clean_init_massfrac, new rhoY equal old one + CVODE increment
  //   const amrex::Real rhoYn =
  //     (clean_init_massfrac)
  //       ? rhoY(i, j, k, n) + (y_vect[vec_index<OrderType>(n, icell, ncells)] -
  //                             rhoY_oldclipped[n]) *
  //                              rho_old / rho_oldclipped
  //       : y_vect[vec_index<OrderType>(n, icell, ncells)];
  //   rhoY(i, j, k, n) = rhoYn;
  //   rho += rhoYn;
  //   mass_frac[n] = rhoYn;
  // }
  // const amrex::Real rho_inv = 1.0 / rho;
  // for (amrex::Real& Yn : mass_frac) {
  //   Yn *= rho_inv;
  // }

  for (int n = 0; n < NUM_SPECIES; n++) nspec(i,j,k,n) = y_vect[vec_index<OrderType>(n, icell, ncells)];
  Ue(i, j, k, 0) = y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];
  Te(i, j, k, 0) = 1.5 * Ue(i,j,k,0) * K_B * nspec(i,j,k,E_ID);

  // amrex::Real temp = y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];
  // rhoE(i, j, k, 0) = vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  // amrex::Real Enrg_loc = rhoE(i, j, k, 0) * rho_inv;
  // auto eos = pele::physics::PhysicsType::eos();
  // if (reactor_type == ReactorTypes::e_reactor_type) {
  //   eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  // } else if (reactor_type == ReactorTypes::h_reactor_type) {
  //   eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  // } else {
  //   amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  // }
  // temperature(i, j, k, 0) = temp;
  // FC_in(i, j, k, 0) = static_cast<amrex::Real>(FCunt);
}

template <typename OrderType>
class FlattenOps
{
public:
  static_assert(
    is_ordering_type<OrderType>::value,
    "FlattenOps must have ordering type as its template argument");

  void flatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& nspec,
    amrex::Array4<const amrex::Real> const& Ue,
    amrex::Array4<const amrex::Real> const& EN,
    amrex::Real* EN_vect,
    amrex::Real* y_vect)
  {
    BL_PROFILE("Pele::flatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_flatten<OrderType>(
        icell, i, j, k, ncells, nspec, Ue, EN, EN_vect, y_vect);
    });
  }

  void unflatten(
    const amrex::Box& box,
    const int ncells,
    amrex::Array4<amrex::Real> const& nspec,
    amrex::Array4<amrex::Real> const& Ue,
    amrex::Array4<amrex::Real> const& Te,
    amrex::Real* y_vect,
    amrex::Real dt)
  {
    BL_PROFILE("Pele::unflatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_unflatten<OrderType>(
        icell, i, j, k, ncells, nspec, Ue, Te,  y_vect, dt);
    });
  }
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
fKernelSpec(
  const int icell,
  const int ncells,
  const amrex::Real dt_save,
  const amrex::Real* EN_in,
  const amrex::Real Tgas_in,
  const amrex::Real* yvec_d,
  amrex::Real* ydot_d // NOLINT(readability-non-const-parameter)
)
{
  amrex::Real NA = 6.02214076e23;
  amrex::Real specconc[NUM_SPECIES];
  amrex::Real wdot[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; n++) {
    specconc[n] = yvec_d[vec_index<OrderType>(n, icell, ncells)] / NA;
  }

  amrex::Real Ue_pt =
    yvec_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];

  amrex::Real ne_pt =
    yvec_d[vec_index<OrderType>(E_ID, icell, ncells)];

  amrex::Real EN_pt = EN_in[icell];

  // Calculate Te from Ue
  amrex::Real Te_in = Ue_pt / (1.5 * K_B * ne_pt);

  amrex::Real ener_exch = 0.0;
  CKWC(Tgas_in, specconc, wdot, Te_in, EN_pt, &ener_exch);

  for (int n = 0; n < NUM_SPECIES; n++) ydot_d[vec_index<OrderType>(n, icell, ncells)] = wdot[n];

  ydot_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] = ener_exch;
}

#ifdef AMREX_USE_GPU
N_Vector
setNVectorGPU(int nvsize, int atomic_reduction, amrex::gpuStream_t stream);
#endif
} // namespace pele::physics::reactions::utils
#endif

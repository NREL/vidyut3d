#ifndef REACTORUTILSCVODE_H
#define REACTORUTILSCVODE_H
#include "ReactorTypes.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "ReactorUtils.H"

#ifdef AMREX_USE_CUDA
#include <cusolverSp.h>
#endif

#ifdef PELE_USE_KLU
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#endif

#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA)
#define CVODE_NB_THREADS 32
#elif defined(AMREX_USE_HIP)
#define CVODE_NB_THREADS 64
#else
#define CVODE_NB_THREADS 32
#endif
#endif

namespace pele::physics::reactions {

struct CVODEUserData
{
  amrex::Real dt_save;     // Internal cvode dt holder
  int ncells;              // Number of cells in solve
  int verbose;             // Verbose
  int maxOrder;            // CVODE max order
  int reactor_type;        // Either HP (LM) or UV (C)
  int solve_type;          // Type of linear solve for Newton direction
  int NNZ; // Number of non-zero entry for sparse representations
  int* mask =
    nullptr; // Masking tagging cells where integration should not be performed
  int* FCunt = nullptr; // Number of RHS evaluations (not used on GPU)
  amrex::Real gamma;    // System Chem. jacobian coefficient
  int nbBlocks;         // GPU kernel launch parameter
  int nbThreads;        // GPU kernel launch parameter
  amrex::Real Tgas_vect; 
  amrex::Real* EN_vect = nullptr;
};

namespace cvode {

#ifdef AMREX_USE_GPU
enum linSolveType {
  fixedPoint,
  sparseDirect,
  customDirect,
  magmaDirect,
  GMRES,
  precGMRES
};
enum precondType { sparseSimpleAJac };
#else
enum linSolveType {
  fixedPoint,
  denseFDDirect,
  denseDirect,
  sparseDirect,
  customDirect,
  GMRES,
  precGMRES,
  hackDumpSparsePattern
};
enum precondType { denseSimpleAJac, sparseSimpleAJac, customSimpleAJac };
#endif

// Error function for CVODE
void cvodeErrHandler(
  int error_code,
  const char* /*module*/,
  const char* /*function*/,
  char* msg,
  void* /*eh_data*/);

} // namespace cvode
} // namespace pele::physics::reactions
#endif

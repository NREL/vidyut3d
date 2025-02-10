#include "ReactorCvode.H"

#include <iostream>

namespace pele::physics::reactions {

int
ReactorCvode::init(int reactor_type, int /*ncells*/)
{
  BL_PROFILE("Pele::ReactorCvode::init()");

  amrex::Print() << "Initializing CVODE:\n";

  // Only parsing/checks are performed here, no actual initialization of
  // the SUNDIALs CVode object.
  m_reactor_type = reactor_type;
  ReactorTypes::check_reactor_type(m_reactor_type);
  amrex::ParmParse pp("vidyut_ode");
  pp.query("verbose", verbose);
  pp.query("rtol", relTol);
  pp.query("atol", absTol);
  pp.query("atomic_reductions", atomic_reductions);
  pp.query("max_nls_iters", max_nls_iters);
  pp.query("max_fp_accel", max_fp_accel);
  pp.query("print_profiling", m_print_profiling);

  // Query CVODE options
  amrex::ParmParse ppcv("cvode");
  ppcv.query("max_order", m_cvode_maxorder);
  ppcv.query("max_substeps", m_cvode_maxstep);
  std::string linear_solve_type;
  ppcv.query("solve_type", linear_solve_type);
  std::string precondJFNK_type;
  ppcv.query("precond_type", precondJFNK_type);

  // Checks
  checkCvodeOptions(
    linear_solve_type, precondJFNK_type, m_solve_type, m_analytical_jacobian,
    m_precond_type);

  if (verbose > 0) {
    if (atomic_reductions != 0) {
      amrex::Print() << "  Using atomic reductions\n";
    } else {
      amrex::Print() << "  Using LDS reductions\n";
    }
  }

  return (0);
}

#ifdef AMREX_USE_GPU
int
ReactorCvode::initCvode(
  N_Vector& a_y,
  SUNMatrix& a_A,
  CVODEUserData* a_udata,
  SUNNonlinearSolver& a_NLS,
  SUNLinearSolver& a_LS,
  void* a_cvode_mem,
  amrex::gpuStream_t stream,
  const amrex::Real& a_time,
  const int ncells)
{
  int flag = CVodeSetUserData(a_cvode_mem, static_cast<void*>(a_udata));

  // Call CVodeInit to initialize the integrator memory and specify the user's
  // right hand side function, the initial time, and initial dependent variable
  // vector a_y.
  flag = CVodeInit(a_cvode_mem, cF_RHS, a_time, a_y);
  if (utils::check_flag(&flag, "CVodeInit", 1)) {
    return (1);
  }

  // NOTE: For simplicity, GMRES is only LS option in Vidyut
  // Solver data
  if (a_udata->solve_type == cvode::GMRES) {
    a_LS = SUNLinSol_SPGMR(
      a_y, SUN_PREC_NONE, 0, *amrex::sundials::The_Sundials_Context());
    if (utils::check_flag(static_cast<void*>(a_LS), "SUNLinSol_SPGMR", 0)) {
      return (1);
    }
    flag = CVodeSetLinearSolver(a_cvode_mem, a_LS, nullptr);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1)) {
      return (1);
    }
    flag = CVodeSetJacTimes(a_cvode_mem, nullptr, nullptr);
    if (utils::check_flag(&flag, "CVodeSetJacTimes", 1)) {
      return (1);
    }
  }

  // CVODE runtime options
  flag = CVodeSetMaxNonlinIters(a_cvode_mem, max_nls_iters);
  if (utils::check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) {
    return (1);
  }
  flag = CVodeSetMaxNumSteps(a_cvode_mem, m_cvode_maxstep);
  if (utils::check_flag(&flag, "CVodeSetMaxNumSteps", 1)) {
    return (1);
  }
  flag = CVodeSetMaxOrd(a_cvode_mem, m_cvode_maxorder);
  if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1)) {
    return (1);
  }
  if (a_LS != nullptr) {
    flag = CVodeSetJacEvalFrequency(a_cvode_mem, 100); // Max Jac age
    if (utils::check_flag(&flag, "CVodeSetJacEvalFrequency", 1) != 0) {
      return (1);
    }
  }

  return (0);
}

#else

int
ReactorCvode::initCvode(
  N_Vector& a_y,
  SUNMatrix& a_A,
  CVODEUserData* a_udata,
  SUNNonlinearSolver& a_NLS,
  SUNLinearSolver& a_LS,
  void* a_cvode_mem,
  const amrex::Real& a_time,
  const int ncells,
  const amrex::Real Tgas_in,
  const amrex::Real Te_in,
  const amrex::Real EN_in)
{
  // Solution vector
  int neq_tot = (NUM_SPECIES + 1) * ncells;
  a_y = N_VNew_Serial(neq_tot, *amrex::sundials::The_Sundials_Context());
  if (utils::check_flag(static_cast<void*>(a_y), "N_VNew_Serial", 0) != 0) {
    return (1);
  }

  // Populate the userData
  allocUserData(a_udata, ncells, Tgas_in, Te_in, EN_in);
  if (utils::check_flag(static_cast<void*>(a_udata), "allocUserData", 2) != 0) {
    return (1);
  }

  // Set the pointer to user-defined data
  int flag = CVodeSetUserData(a_cvode_mem, a_udata);
  if (utils::check_flag(&flag, "CVodeSetUserData", 1) != 0) {
    return (1);
  }

  // Call CVodeInit to initialize the integrator memory and specify the user's
  // right hand side function, the initial time, and initial dependent variable
  // vector a_y.
  flag = CVodeInit(a_cvode_mem, cF_RHS, a_time, a_y);
  if (utils::check_flag(&flag, "CVodeInit", 1) != 0) {
    return (1);
  }

  // Linear solver data
  if (a_udata->solve_type == cvode::GMRES) {
    // Create the GMRES linear solver object
    a_LS = SUNLinSol_SPGMR(
      a_y, SUN_PREC_NONE, 0, *amrex::sundials::The_Sundials_Context());
    if (
      utils::check_flag(static_cast<void*>(a_LS), "SUNLinSol_SPGMR", 0) != 0) {
      return (1);
    }

    // Set CVode linear solver to LS
    flag = CVodeSetLinearSolver(a_cvode_mem, a_LS, nullptr);
    if (utils::check_flag(&flag, "CVodeSetLinearSolver", 1) != 0) {
      return (1);
    }
  } else {
    amrex::Abort("Wrong choice of linear solver");
  }

  // CVODE runtime options
  flag = CVodeSetMaxNonlinIters(a_cvode_mem, max_nls_iters); // Max newton iter.
  if (utils::check_flag(&flag, "CVodeSetMaxNonlinIters", 1) != 0) {
    return (1);
  }
  flag = CVodeSetMaxErrTestFails(a_cvode_mem, 100); // Max Err.test failure
  if (utils::check_flag(&flag, "CVodeSetMaxErrTestFails", 1) != 0) {
    return (1);
  }
  flag = CVodeSetMaxNumSteps(a_cvode_mem, 10000); // Max substeps
  if (utils::check_flag(&flag, "CVodeSetMaxNumSteps", 1) != 0) {
    return (1);
  }
  flag = CVodeSetMaxOrd(a_cvode_mem, m_cvode_maxorder); // Max order
  if (utils::check_flag(&flag, "CVodeSetMaxOrd", 1) != 0) {
    return (1);
  }
  if (a_LS != nullptr) {
    flag = CVodeSetJacEvalFrequency(a_cvode_mem, 100); // Max Jac age
    if (utils::check_flag(&flag, "CVodeSetJacEvalFrequency", 1) != 0) {
      return (1);
    }
  }
  return (0);
}

#endif // End check GPU for initCvode method

void
ReactorCvode::checkCvodeOptions(
  const std::string& a_solve_type_str,
  const std::string& a_precond_type_str,
  int& a_solve_type,
  int& a_ajac,
  int& a_precond_type) const
{
  if (verbose > 0) {
    amrex::Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
#ifdef PELE_CVODE_FORCE_YCORDER
    amrex::Print() << "Using YCOrder\n";
#else
    amrex::Print() << "Using CYOrder\n";
#endif
  }

  //-------------------------------------------------------------
  // Shared CPU/GPU options
  //-------------------------------------------------------------
  if (a_solve_type_str == "GMRES") {
    a_solve_type = cvode::GMRES;
    if (verbose > 0) {
      amrex::Print() << " Using a JFNK GMRES linear solve\n";
    }
  } else {
#ifdef AMREX_USE_GPU
    amrex::Abort(
      "Wrong solve_type. Options are: 'GMRES'");
#else
    amrex::Abort(
      "Wrong solve_type. Options are: 'GMRES'");
#endif
  }
}

void
ReactorCvode::allocUserData(
  CVODEUserData* udata,
  int a_ncells,
#ifdef AMREX_USE_GPU
  ,
  SUNMatrix& a_A,
  amrex::gpuStream_t stream
#endif
) const
{
  // Pass options to udata
  udata->solve_type = m_solve_type;

  const int HP =
    static_cast<int>(m_reactor_type == ReactorTypes::h_reactor_type);
  int nspec_tot = (NUM_SPECIES)*a_ncells;
  udata->reactor_type = m_reactor_type;
  udata->ncells = a_ncells;
  udata->verbose = verbose;
#ifdef AMREX_USE_GPU
  udata->nbThreads = CVODE_NB_THREADS;
  udata->nbBlocks = std::max(1, a_ncells / udata->nbThreads);
  udata->stream = stream;
#endif

  // Alloc internal udata solution/forcing containers
  udata->Ue_init = static_cast<amrex::Real*>(
    amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real)));
  udata->Tgas_init = static_cast<amrex::Real*>(
    amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real)));
  udata->Te_init = static_cast<amrex::Real*>(
    amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real)));
  udata->EN_init = static_cast<amrex::Real*>(
    amrex::The_Arena()->alloc(a_ncells * sizeof(amrex::Real)));
  udata->mask =
    static_cast<int*>(amrex::The_Arena()->alloc(a_ncells * sizeof(int)));

#ifndef AMREX_USE_GPU
  udata->FCunt =
    static_cast<int*>(amrex::The_Arena()->alloc(a_ncells * sizeof(int)));
#endif
}

int
ReactorCvode::react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& nspec_in,
  amrex::Array4<amrex::Real> const& Tgas_in,
  amrex::Array4<amrex::Real> const& Te_in,
  amrex::Array4<amrex::Real> const& EN_in,
  amrex::Array4<amrex::Real> const& Ue_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time,
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
)
{
  BL_PROFILE("Pele::ReactorCvode::react()");

  // CPU and GPU version are very different such that most of the function
  // is split between a GPU region and a CPU region

  amrex::Real time_start = time;
  amrex::Real time_final = time + dt_react;
  amrex::Real CvodeActual_time_final = 0.0;

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  SUNProfiler sun_profiler = nullptr;
  SUNContext_GetProfiler(
    *amrex::sundials::The_Sundials_Context(), &sun_profiler);
#endif

  // Set of SUNDIALs objects needed for Cvode
  SUNMatrix A = nullptr;             // Jacobian matrix
  auto* udata = new CVODEUserData{}; // Userdata container
  SUNNonlinearSolver NLS = nullptr;  // Non-linear solver
  SUNLinearSolver LS = nullptr;      // Linear solver

  // Call CVodeCreate to create the solver memory and specify the Backward
  // Differentiation Formula and the use of a Newton iteration
  void* cvode_mem =
    CVodeCreate(CV_BDF, *amrex::sundials::The_Sundials_Context());
  ; // Internal Cvode memory

  //----------------------------------------------------------
  // GPU Region
  //----------------------------------------------------------

#ifdef AMREX_USE_GPU
  const int ncells = box.numPts();
  const int neq_tot = (NUM_SPECIES + 1) * ncells;

  // Solution vector and execution policy
  auto y = utils::setNVectorGPU(neq_tot, atomic_reductions, stream);

  // Solution data array
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y);

  // Populate the userData
  amrex::Gpu::streamSynchronize();
  allocUserData(udata, ncells, A, stream);

  // Fill data
  flatten(
    box, ncells, nspec_in, Tgas_in, Te_in, EN_in, Ue_in, yvec_d, udata->Ue_init);

#ifdef AMREX_USE_OMP
  amrex::Gpu::Device::streamSynchronize();
#endif

  initCvode(y, A, udata, NLS, LS, cvode_mem, stream, time_start, ncells, Tgas_in, Te_in, EN_in);

  // Setup tolerances with typical values
  utils::set_sundials_solver_tols<Ordering>(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, udata->ncells, relTol,
    absTol, m_typ_vals, "cvode", verbose);

  // Actual CVODE solve
  BL_PROFILE_VAR("Pele::ReactorCvode::react():CVode", AroundCVODE);
  int flag =
    CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
  if (utils::check_flag(&flag, "CVode", 1)) {
    return (1);
  }
  BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
  dt_react =
    time_start - CvodeActual_time_final; // Actual dt_react performed by Cvode
  time += dt_react;                      // Increment time in reactor mode
#endif

#ifdef AMREX_USE_OMP
  amrex::Gpu::Device::streamSynchronize();
#endif

  // Get workload estimate
  long int nfe;
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

  amrex::Gpu::DeviceVector<long int> v_nfe(ncells, nfe);
  long int* d_nfe = v_nfe.data();
  unflatten(
    box, ncells, nspec_in, Te_in, Ue_in, FC_in, yvec_d,
    udata->Ue_init, d_nfe, dt_react);

  if (udata->verbose > 1) {
    print_final_stats(cvode_mem, LS != nullptr);
  }

#else
  //----------------------------------------------------------
  // CPU Region
  //----------------------------------------------------------

  N_Vector y = nullptr; // Solution vector

  // Perform integration one cell at a time
  const int icell = 0;
  const int ncells = 1;

  int omp_thread = 0;
#ifdef AMREX_USE_OMP
  omp_thread = omp_get_thread_num();
#endif

  initCvode(y, A, udata, NLS, LS, cvode_mem, time_start, ncells);

  // Update TypicalValues
  // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
  utils::set_sundials_solver_tols<Ordering>(
    *amrex::sundials::The_Sundials_Context(), cvode_mem, udata->ncells, relTol,
    absTol, m_typ_vals, "cvode", verbose);

  const auto captured_reactor_type = m_reactor_type;
  const auto captured_clean_init_massfrac = m_clean_init_massfrac;
  ParallelFor(
    box, [=, &CvodeActual_time_final] AMREX_GPU_DEVICE(
           int i, int j, int k) noexcept {
      if (mask(i, j, k) != -1) {

        amrex::Real* yvec_d = N_VGetArrayPointer(y);
        utils::box_flatten<Ordering>(
          icell, i, j, k, ncells, nspec_in, Tgas_in, Te_in, EN_in, Ue_in, yvec_d, udata->Ue_init);

        // ReInit CVODE is faster
        CVodeReInit(cvode_mem, time_start, y);

        BL_PROFILE_VAR("Pele::ReactorCvode::react():CVode", AroundCVODE);
        CVode(cvode_mem, time_final, y, &CvodeActual_time_final, CV_NORMAL);
        BL_PROFILE_VAR_STOP(AroundCVODE);

        // cppcheck-suppress knownConditionTrueFalse
        if ((udata->verbose > 1) && (omp_thread == 0)) {
          amrex::Print() << "Additional verbose info --\n";
          print_final_stats(cvode_mem, LS != nullptr);
          amrex::Print() << "\n -------------------------------------\n";
        }

        amrex::Real actual_dt = CvodeActual_time_final - time_start;

        // Get estimate of how hard the integration process was
        long int nfe = 0;
        long int nfeLS = 0;
        CVodeGetNumRhsEvals(cvode_mem, &nfe);
        if (LS != nullptr) {
          CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
        }
        const long int nfe_tot = nfe + nfeLS;

        utils::box_unflatten<Ordering>(
          icell, i, j, k, ncells, nspec_in, Te_in, Ue_in,
          FC_in, yvec_d, udata->Ue_init, nfe_tot, dt_react);

        // cppcheck-suppress knownConditionTrueFalse
        if ((udata->verbose > 3) && (omp_thread == 0)) {
          amrex::Print() << "END : time curr is " << CvodeActual_time_final
                         << " and actual dt_react is " << actual_dt << "\n";
        }
      } else {
        FC_in(i, j, k, 0) = 0.0;
      }
    });

#ifdef MOD_REACTOR
  dt_react =
    time_start -
    time_final; // In this case, assumes all individual CVODE calls nailed it.
  time += dt_react; // Increment time in reactor mode
#endif

  long int nfe =
    20; // Dummy, the return value is no longer used for this function.

#endif // End GPU check

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  if (m_print_profiling) {
    SUNProfiler_Print(sun_profiler, stdout);
  }
#endif

  // Clean up
  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  if (LS != nullptr) {
    SUNLinSolFree(LS);
  }
  if (NLS != nullptr) {
    SUNNonlinSolFree(NLS);
  }
  if (A != nullptr) {
    SUNMatDestroy(A);
  }
  freeUserData(udata);

  return static_cast<int>(nfe);
}

int
ReactorCvode::cF_RHS(
  sunrealtype t, N_Vector y_in, N_Vector ydot_in, void* user_data)
{
  BL_PROFILE("Pele::ReactorCvode::cF_RHS()");
#ifdef AMREX_USE_GPU
  amrex::Real* yvec_d = N_VGetDeviceArrayPointer(y_in);
  amrex::Real* ydot_d = N_VGetDeviceArrayPointer(ydot_in);
#else
  amrex::Real* yvec_d = N_VGetArrayPointer(y_in);
  amrex::Real* ydot_d = N_VGetArrayPointer(ydot_in);
#endif

  auto* udata = static_cast<CVODEUserData*>(user_data);
  udata->dt_save = t;

  const auto ncells = udata->ncells;
  const auto dt_save = udata->dt_save;
  const auto reactor_type = udata->reactor_type;
  auto* Ue_init = udata->Ue_init;
  amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(int icell) noexcept {
    utils::fKernelSpec<Ordering>(
      icell, ncells, dt_save, yvec_d, ydot_d, Ue_init, Tgas_in, Te_in, EN_in);
  });
  amrex::Gpu::Device::streamSynchronize();
  return 0;
}

void
ReactorCvode::freeUserData(CVODEUserData* data_wk)
{
  amrex::The_Arena()->free(data_wk->Ue_init);
  amrex::The_Arena()->free(data_wk->mask);

#ifdef AMREX_USE_GPU

  if (data_wk->solve_type == cvode::sparseDirect) {
#ifdef AMREX_USE_CUDA
    amrex::The_Pinned_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Pinned_Arena()->free(data_wk->csr_col_index_h);
    cusolverStatus_t cusolver_status =
      cusolverSpDestroy(data_wk->cusolverHandle);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusparseStatus_t cusparse_status = cusparseDestroy(data_wk->cuSPHandle);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
#endif
  } else if (data_wk->solve_type == cvode::customDirect) {
#ifdef AMREX_USE_CUDA
    amrex::The_Pinned_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Pinned_Arena()->free(data_wk->csr_col_index_h);
    cusparseStatus_t cusparse_status = cusparseDestroy(data_wk->cuSPHandle);
    AMREX_ASSERT(cusparse_status == CUSPARSE_STATUS_SUCCESS);
#endif
  }
  // Preconditioner analytical Jacobian data
  if (data_wk->precond_type == cvode::sparseSimpleAJac) {
#ifdef AMREX_USE_CUDA
    amrex::The_Pinned_Arena()->free(data_wk->csr_row_count_h);
    amrex::The_Pinned_Arena()->free(data_wk->csr_col_index_h);
    amrex::The_Arena()->free(data_wk->csr_row_count_d);
    amrex::The_Arena()->free(data_wk->csr_col_index_d);
    amrex::The_Arena()->free(data_wk->csr_val_d);
    amrex::The_Arena()->free(data_wk->csr_jac_d);
    cusolverStatus_t cusolver_status =
      cusolverSpDestroy(data_wk->cusolverHandle);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusolver_status = cusolverSpDestroyCsrqrInfo(data_wk->info);
    AMREX_ASSERT(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cudaFree(data_wk->buffer_qr);
#endif
  }
  delete data_wk;

#else
  amrex::The_Arena()->free(data_wk->FCunt);

  // Direct solver Jac. data
  if (data_wk->solve_type == cvode::sparseDirect) {
#ifdef PELE_USE_KLU
    delete[] data_wk->colPtrs;
    delete[] data_wk->rowVals;
    delete[] data_wk->Jdata;
    SUNMatDestroy(A);
    SUNMatDestroy((data_wk->PS)[0]);
    delete[] (data_wk->PS);
#endif
  } else if (data_wk->solve_type == cvode::customDirect) {
    SUNMatDestroy(data_wk->PSc);
  }

  // Preconditionner Jac. data
  if (data_wk->precond_type == cvode::denseSimpleAJac) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      SUNDlsMat_destroyMat((data_wk->P)[i][i]);
      SUNDlsMat_destroyMat((data_wk->Jbd)[i][i]);
      SUNDlsMat_destroyArray((data_wk->pivot)[i][i]);
    }
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[] (data_wk->P)[i];
      delete[] (data_wk->Jbd)[i];
      delete[] (data_wk->pivot)[i];
    }
    delete[] (data_wk->P);
    delete[] (data_wk->Jbd);
    delete[] (data_wk->pivot);
  } else if (data_wk->precond_type == cvode::sparseSimpleAJac) {
#ifdef PELE_USE_KLU
    delete[] data_wk->colPtrs;
    delete[] data_wk->rowVals;
    delete[] data_wk->Jdata;
    for (int i = 0; i < data_wk->ncells; ++i) {
      klu_free_symbolic(&(data_wk->Symbolic[i]), &(data_wk->Common[i]));
      klu_free_numeric(&(data_wk->Numeric[i]), &(data_wk->Common[i]));
      delete[] data_wk->JSPSmat[i];
      SUNMatDestroy((data_wk->PS)[i]);
    }
    delete[] data_wk->JSPSmat;
    delete[] data_wk->Common;
    delete[] data_wk->Symbolic;
    delete[] data_wk->Numeric;
    delete[] data_wk->PS;
#endif
  } else if (data_wk->precond_type == cvode::customSimpleAJac) {
    for (int i = 0; i < data_wk->ncells; ++i) {
      delete[] data_wk->JSPSmat[i];
      SUNMatDestroy((data_wk->PS)[i]);
    }
    delete[] data_wk->colVals;
    delete[] data_wk->rowPtrs;
    delete[] data_wk->PS;
    delete[] data_wk->JSPSmat;
  }

  delete data_wk;
#endif
}

void
ReactorCvode::close()
{
}

void
ReactorCvode::print_final_stats(void* cvodemem, bool print_ls_stats) // NOLINT
{
  long int nst, nfe, nsetups, nni, ncfn, netf, nje;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  // CVODE stats
  flag = CVodeGetNumSteps(cvodemem, &nst);
  utils::check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumErrTestFails(cvodemem, &netf);
  utils::check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumRhsEvals(cvodemem, &nfe);
  utils::check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  // Nonlinear solver stats
  flag = CVodeGetNumNonlinSolvIters(cvodemem, &nni);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodemem, &ncfn);
  utils::check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  if (print_ls_stats) {
    // Linear solver stats
    flag = CVodeGetNumLinSolvSetups(cvodemem, &nsetups);
    utils::check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumJacEvals(cvodemem, &nje);
    utils::check_flag(&flag, "CVodeGetNumJacEvals", 1);
    flag = CVodeGetNumLinIters(cvodemem, &nli);
    utils::check_flag(&flag, "CVodeGetNumLinIters", 1);
    flag = CVodeGetNumLinConvFails(cvodemem, &ncfl);
    utils::check_flag(&flag, "CVodeGetNumLinConvFails", 1);
    flag = CVodeGetNumLinRhsEvals(cvodemem, &nfeLS);
    utils::check_flag(&flag, "CVodeGetNumLinRhsEvals", 1);
    // Preconditioner stats
    flag = CVodeGetNumPrecEvals(cvodemem, &npe);
    utils::check_flag(&flag, "CVodeGetNumPrecEvals", 1);
    flag = CVodeGetNumPrecSolves(cvodemem, &nps);
    utils::check_flag(&flag, "CVodeGetNumPrecSolves", 1);
  }

#ifdef AMREX_USE_OMP
  amrex::Print() << "\nFinal Statistics: "
                 << "(thread:" << omp_get_thread_num() << ", ";
  amrex::Print() << "cvode_mem:" << cvodemem << ")\n";
#else
  amrex::Print() << "\nFinal Statistics:\n";
#endif
  // CVODE stats
  amrex::Print() << "  nSteps       = " << nst << "\n";
  amrex::Print() << "  nErrtf       = " << netf << "\n";
  amrex::Print() << "  nRHSeval     = " << nfe << "\n";
  // NLS stats
  amrex::Print() << "  nnLinIt      = " << nni << "\n";
  amrex::Print() << "  nConvfail    = " << ncfn << "\n";
  // LS stats
  if (print_ls_stats) {
    amrex::Print() << "  nLinsetups   = " << nsetups << "\n";
    amrex::Print() << "  nJeval       = " << nje << "\n";
    amrex::Print() << "  nLinIt       = " << nli << "\n";
    amrex::Print() << "  nLinConvfail = " << ncfl << "\n";
    amrex::Print() << "  nLinRHSeval  = " << nfeLS << "\n";
    // Prec
    amrex::Print() << "  nPreceval    = " << npe << "\n";
    amrex::Print() << "  nPrecsolve   = " << nps << "\n";
  }
}

} // namespace pele::physics::reactions

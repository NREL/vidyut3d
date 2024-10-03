#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <UserFunctions.H>
#include <UserSources.H>
#include <compute_explicit_flux.H>
#include <AMReX_MLABecLaplacian.H>

void Vidyut::compute_dsdt(int lev, int specid, 
                            Array<MultiFab,AMREX_SPACEDIM>& flux, 
                            MultiFab& rxn_src,
                            MultiFab& dsdt,
                            Real time, Real dt)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int captured_specid = specid;
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real> rxn_arr = rxn_src.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);

        GpuArray<Array4<Real>, AMREX_SPACEDIM> flux_arr{
            AMREX_D_DECL(flux[0].array(mfi), 
                         flux[1].array(mfi), flux[2].array(mfi))};

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            dsdt_arr(i, j, k) = (flux_arr[0](i, j, k) - flux_arr[0](i + 1, j, k)) / dx[0] 
            + rxn_arr(i,j,k,captured_specid);
#if AMREX_SPACEDIM > 1
            dsdt_arr(i,j,k) += (flux_arr[1](i, j, k) - flux_arr[1](i, j + 1, k)) / dx[1];
#if AMREX_SPACEDIM == 3
            dsdt_arr(i,j,k) += (flux_arr[2](i, j, k) - flux_arr[2](i, j, k + 1)) / dx[2]; 
#endif
#endif
        });
    }
}

void Vidyut::update_explsrc_at_all_levels(int specid, Vector<MultiFab>& Sborder,
                                            Vector<Array<MultiFab,AMREX_SPACEDIM>>& flux,
                                            Vector<MultiFab>& rxn_src, 
                                            Vector<MultiFab>& expl_src, 
                                            Vector<int>& bc_lo, Vector<int>& bc_hi,
                                            amrex::Real cur_time)
{
    for(int lev=0; lev <= finest_level; lev++)
    {
        expl_src[lev].setVal(0.0);
        flux[lev][0].setVal(0.0);
#if AMREX_SPACEDIM > 1
        flux[lev][1].setVal(0.0);
#if AMREX_SPACEDIM == 3
        flux[lev][2].setVal(0.0);
#endif
#endif
    }

    if(do_transport)
    {
        for(int lev=0; lev <= finest_level; lev++)
        {
            compute_scalar_transport_flux(lev, Sborder[lev], 
                                          flux[lev], bc_lo, bc_hi, 
                                          cur_time, specid);
        }
    }

    // =======================================================
    // Average down the fluxes before using them to update phi
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
        average_down_faces(amrex::GetArrOfConstPtrs(flux[lev  ]),
                           amrex::GetArrOfPtrs(flux[lev-1]),
                           refRatio(lev-1), Geom(lev-1));
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        //FIXME: need to avoid this fillpatch
        compute_dsdt(lev, specid, 
                     flux[lev], rxn_src[lev], expl_src[lev], 
                     cur_time, dt[lev]);
    }

    // Additional source terms for axisymmetric geometry
    if(geom[0].IsRZ()){
        for (int lev = 0; lev <= finest_level; lev++)
        {
            compute_axisym_correction(lev, Sborder[lev], expl_src[lev], cur_time, specid);   
        }
    }
}

void Vidyut::update_rxnsrc_at_all_levels(Vector<MultiFab>& Sborder,
                                         Vector<MultiFab>& rxn_src, 
                                         amrex::Real cur_time)
{
    amrex::Real time = cur_time;
    ProbParm const* localprobparm = d_prob_parm;

    // Zero out reactive source MFs
    for(int lev=0; lev <= finest_level; lev++)
    {
        rxn_src[lev].setVal(0.0);
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        amrex::Real captured_gastemp=gas_temperature;
        amrex::Real captured_gaspres=gas_pressure;
        const auto dx = geom[lev].CellSizeArray();
        auto prob_lo = geom[lev].ProbLoArray();
        auto prob_hi = geom[lev].ProbHiArray();

        for (MFIter mfi(rxn_src[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> sborder_arr = Sborder[lev].array(mfi);
            Array4<Real> rxn_arr = rxn_src[lev].array(mfi);

            // update residual
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Create array with species concentrations (1/m3 -> mol/cm3)
                amrex::Real spec_C[NUM_SPECIES];
                amrex::Real spec_wdot[NUM_SPECIES];
                amrex::Real Te = sborder_arr(i,j,k,ETEMP_ID);
                amrex::Real EN = sborder_arr(i,j,k,REF_ID);
                amrex::Real ener_exch = 0.0;
                for(int sp=0; sp<NUM_SPECIES; sp++) spec_C[sp] = sborder_arr(i,j,k,sp) * 1.0e-6 / N_A;

                // Get molar production rates
                CKWC(captured_gastemp, spec_C, spec_wdot, Te, EN, &ener_exch);

                // Convert from mol/cm3-s to 1/m3-s and add to scalar react source MF
                for(int sp = 0; sp<NUM_SPECIES; sp++) rxn_arr(i,j,k,sp) = spec_wdot[sp] * N_A * 1.0e6;
                rxn_arr(i,j,k,NUM_SPECIES) = ener_exch;

                // Add on user-defined reactive sources
                user_sources::add_user_react_sources
                (i, j, k, sborder_arr, rxn_arr,
                 prob_lo, prob_hi, dx, time, *localprobparm,
                 captured_gastemp,
                 captured_gaspres);
            });
        }
    }
}

void Vidyut::compute_scalar_transport_flux(int lev, MultiFab& Sborder, 
                                           Array<MultiFab,AMREX_SPACEDIM>& flux, 
                                           Vector<int>& bc_lo, Vector<int>& bc_hi,
                                           Real current_time,int specid)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int captured_specid = specid;
    //class member variable
    int captured_hyporder = hyp_order; 
    int captured_wenoscheme = weno_scheme;
    int userdefvel = user_defined_vel;

    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;
    amrex::Real lev_dt=dt[lev];

    // Get the boundary ids
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();

    GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
    GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};
    
    GpuArray<int,AMREX_SPACEDIM> bclo={AMREX_D_DECL(bc_lo[0], bc_lo[1], bc_lo[2])};
    GpuArray<int,AMREX_SPACEDIM> bchi={AMREX_D_DECL(bc_hi[0], bc_hi[1], bc_hi[2])};


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Box bx_x = convert(bx, {AMREX_D_DECL(1, 0, 0)});
#if AMREX_SPACEDIM > 1
            Box bx_y = convert(bx, {AMREX_D_DECL(0, 1, 0)});
#if AMREX_SPACEDIM == 3
            Box bx_z = convert(bx, {0, 0, 1});
#endif
#endif
            
            Real time = current_time; // for GPU capture

            Array4<Real> sborder_arr = Sborder.array(mfi);

            GpuArray<Array4<Real>, AMREX_SPACEDIM> 
            flux_arr{AMREX_D_DECL(flux[0].array(mfi), 
                                  flux[1].array(mfi), flux[2].array(mfi))};

            //amrex::Print()<<"bx:"<<bx<<"\n";
            //amrex::Print()<<"bx_x:"<<bx_x<<"\n";
            amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 0, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[0], 
                             captured_gastemp,captured_gaspres,
                             time, dx, lev_dt, *localprobparm, captured_hyporder,
                             userdefvel,captured_wenoscheme); 
            });

#if AMREX_SPACEDIM > 1
            amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 1, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[1], 
                             captured_gastemp,captured_gaspres,
                             time, dx, lev_dt, *localprobparm, captured_hyporder,
                             userdefvel,captured_wenoscheme); 
            });

#if AMREX_SPACEDIM == 3
            amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 2, captured_specid, sborder_arr, 
                             bclo, bchi, domlo, domhi, flux_arr[2], 
                             captured_gastemp, captured_gaspres,
                             time, dx, lev_dt, *localprobparm, captured_hyporder,
                             userdefvel,weno_scheme);
            });
#endif
#endif
        }
    }
}

void Vidyut::compute_axisym_correction(int lev, MultiFab& Sborder,MultiFab& dsdt,
                                       Real time,int specid)
{
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;
    const auto dx = geom[lev].CellSizeArray();

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real> s_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);

        // Evaluate cell-centered axisymmetric source terms (Gamma_k / r)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          // calculate r
          amrex::Real rval = amrex::Math::abs((i+0.5) * dx[0]);

          // Calculate the advective source term component
          amrex::Real etemp = s_arr(i,j,k,ETEMP_ID);
          amrex::Real ndens = 0.0;
          amrex::Real Esum = 0.0;
          for(int sp=0; sp<NUM_SPECIES; sp++) ndens += s_arr(i,j,k,sp);
          for (int dim = 0; dim < AMREX_SPACEDIM; dim++) Esum += std::pow(s_arr(i,j,k,EFX_ID+dim),2.0);
          amrex::Real efield_mag=std::sqrt(Esum);
          amrex::Real mu = specMob(specid, etemp, ndens, efield_mag,captured_gastemp);  
          dsdt_arr(i,j,k) -= mu * s_arr(i,j,k,specid) * s_arr(i,j,k,EFX_ID) / rval;
        });
    }
}

void Vidyut::implicit_solve_scalar(Real current_time, Real dt, int spec_id, 
                                   Vector<MultiFab>& Sborder, 
                                   Vector<MultiFab>& Sborder_old, 
                                   Vector<MultiFab>& dsdt_expl, 
                                   Vector<int>& bc_lo, Vector<int>& bc_hi,
                                   Vector<Array<MultiFab,AMREX_SPACEDIM>>& grad_fc)
{
    BL_PROFILE("Vidyut::implicit_solve_species(" + std::to_string( spec_id ) + ")");


    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int linsolve_verbose = 1;
    int captured_spec_id=spec_id;
    int electron_flag=(spec_id==E_IDX)?1:0;
    int electron_energy_flag=(spec_id==EEN_ID)?1:0;

    //==================================================
    // amrex solves
    // read small a as alpha, b as beta

    //(A a - B del.(b del)) phi = f
    //
    // A and B are scalar constants
    // a and b are scalar fields
    // f is rhs
    // in this case: A=0,a=0,B=1,b=conductivity
    // note also the negative sign
    //====================================================
    ProbParm const* localprobparm = d_prob_parm;

    const Real tol_rel = linsolve_reltol;
    const Real tol_abs = linsolve_abstol;

    // set A and B, A=1/dt, B=1
    Real ascalar = 1.0;
    Real bscalar = 1.0;
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;
    int userdefspec = user_defined_species;
    int eidx = E_IDX;

#ifdef AMREX_USE_HYPRE
    if(use_hypre)
    {
        amrex::Print()<<"using hypre\n";
    }
#endif

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_lo 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_hi 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (bc_lo[idim] == PERBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (bc_lo[idim] == DIRCBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_lo[idim] == HNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }
        if (bc_lo[idim] == IHNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_lo[idim] == ROBINBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }

        //higher side bcs
        if (bc_hi[idim] == PERBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Periodic;
        }
        if (bc_hi[idim] == DIRCBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_hi[idim] == HNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
        if (bc_hi[idim] == IHNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_hi[idim] == ROBINBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
    }

    Vector<MultiFab> specdata(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<MultiFab> bcoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        specdata[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        bcoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    linsolve_ptr.reset(new MLABecLaplacian(Geom(0,finest_level), 
                                           boxArray(0,finest_level), 
                                           DistributionMap(0,finest_level), info));
    MLMG mlmg(*linsolve_ptr);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(linsolve_verbose);

#ifdef AMREX_USE_HYPRE
        if (use_hypre)
        {
            mlmg.setHypreOptionsNamespace("vidyut.hypre");
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        }
#endif
    linsolve_ptr->setDomainBC(bc_linsolve_lo, bc_linsolve_hi);
    linsolve_ptr->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        // Copy args (FabArray<FAB>& dst, FabArray<FAB> const& src, 
        // int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
        specdata[ilev].setVal(0.0);
        amrex::Copy(specdata[ilev], Sborder_old[ilev], captured_spec_id, 
                    0, 1, num_grow);

        acoeff[ilev].setVal(1.0/dt);
        bcoeff[ilev].setVal(1.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        rhs[ilev].setVal(0.0);

        //adding U^n/dt and explicit sources
        MultiFab::LinComb(rhs[ilev], 1.0/dt, specdata[ilev], 0, 1.0, 
                          dsdt_expl[ilev], 0, 0, 1, 0);
        
        amrex::Copy(specdata[ilev], Sborder[ilev], captured_spec_id, 
                    0, 1, num_grow);

        solution[ilev].setVal(0.0);
        amrex::MultiFab::Copy(solution[ilev], specdata[ilev], 0, 0, 1, 0);

        // fill cell centered diffusion coefficients and rhs
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Real time = current_time; // for GPU capture

            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Array4<Real> acoeff_arr = acoeff[ilev].array(mfi);
            Array4<Real> bcoeff_arr = bcoeff[ilev].array(mfi);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                //FIXME:may be use updated efields here
                    amrex::Real Esum = 0.0;
                    for(int dim=0; dim<AMREX_SPACEDIM; dim++) Esum += std::pow(sb_arr(i,j,k,EFX_ID+dim),2.0);
                    amrex::Real efield_mag=std::sqrt(Esum);
                
                    amrex::Real ndens = 0.0;
                    for(int sp=0; sp<NUM_SPECIES; sp++) ndens += sb_arr(i,j,k,sp);
                    bcoeff_arr(i,j,k)=specDiff(captured_spec_id, 
                                               sb_arr(i,j,k,ETEMP_ID), ndens,
                                               efield_mag, 
                                               captured_gastemp);

            });
        }



        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoeff[ilev].boxArray(), 
                                                IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, bcoeff[ilev].DistributionMap(), 1, 0);
        }
        // true argument for harmonic averaging
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoeff), 
                                          bcoeff[ilev], geom[ilev], true);


        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> bc_arr = specdata[ilev].array(mfi);
            Array4<Real> sb_arr = Sborder[ilev].array(mfi);
            Real time = current_time; // for GPU capture

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (!geom[ilev].isPeriodic(idim))
                {
                    //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                    //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                    //so the ghost cell index at left side is i-1 while it is i on the right
                    if (bx.smallEnd(idim) == domain.smallEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            if(userdefspec == 1){
                                user_transport::species_bc(i, j, k, idim, -1, 
                                                             captured_spec_id, sb_arr, bc_arr, robin_a_arr,
                                                             robin_b_arr, robin_f_arr, 
                                                             prob_lo, prob_hi, dx, time, *localprobparm,
                                                             captured_gastemp,captured_gaspres);
                            } else {
                                plasmachem_transport::species_bc(i, j, k, idim, -1, 
                                                             captured_spec_id, sb_arr, bc_arr, robin_a_arr,
                                                             robin_b_arr, robin_f_arr, 
                                                             prob_lo, prob_hi, dx, time, *localprobparm,
                                                             captured_gastemp,captured_gaspres);
                            }
                        });
                    }
                    if (bx.bigEnd(idim) == domain.bigEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            if(userdefspec == 1){
                                user_transport::species_bc(i, j, k, idim, +1, 
                                                             captured_spec_id, sb_arr, bc_arr, robin_a_arr, 
                                                             robin_b_arr, robin_f_arr,
                                                             prob_lo, prob_hi, dx, time, *localprobparm,
                                                             captured_gastemp,captured_gaspres);
                            } else {
                                plasmachem_transport::species_bc(i, j, k, idim, +1, 
                                                             captured_spec_id, sb_arr, bc_arr, robin_a_arr, 
                                                             robin_b_arr, robin_f_arr,
                                                             prob_lo, prob_hi, dx, time, *localprobparm,
                                                             captured_gastemp,captured_gaspres);
                            }
                        });
                    }
                }
            }
        }

        linsolve_ptr->setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_ptr->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells
        if(mixedbc)
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]), &(robin_a[ilev]), 
                                     &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]));
        }
    }

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    //bound species density
    if(bound_specden && !electron_energy_flag)
    { 
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            amrex::Real minelecden=min_electron_density; 
            amrex::Real minspecden=min_species_density; 
            for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> soln_arr = solution[ilev].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    if(electron_flag){
                      if(soln_arr(i,j,k) < minelecden)
                      {
                        soln_arr(i,j,k)=minelecden;
                      } 
                    } else {
                      if(soln_arr(i,j,k) < minspecden)
                      {
                        soln_arr(i,j,k)=minspecden;
                      } 
                    }
                });
            }
        }
    }
    if(electron_flag)
    {
        mlmg.getGradSolution(GetVecOfArrOfPtrs(grad_fc));
    }

    // copy solution back to phi_new
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, spec_id, 1, 0);
    }
    
    Print()<<"Solved species:"<<allvarnames[spec_id]<<"\n";

    if(electron_energy_flag)
    {
        /*for(int ilev=0; ilev <= finest_level; ilev++)
          {
          phi_new[ilev].setVal(1.0,ETEMP_ID,1);
          amrex::MultiFab::Multiply(phi_new[ilev],solution[ilev],EEN_ID, ETEMP_ID, 1, 0);
          amrex::MultiFab::Divide(phi_new[ilev],phi_new[ilev],EDN_ID, ETEMP_ID, 1, 0);
          phi_new[ilev].mult(twothird/K_B, ETEMP_ID, 1);
          }*/

        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            amrex::Real minetemp=min_electron_temp; 
            for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> phi_arr = phi_new[ilev].array(mfi);
                Array4<Real> sb_arr = Sborder[ilev].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    phi_arr(i,j,k,ETEMP_ID)=twothird/K_B*phi_arr(i,j,k,EEN_ID)/sb_arr(i,j,k,eidx);
                    if(phi_arr(i,j,k,ETEMP_ID) < minetemp)
                    {
                        phi_arr(i,j,k,ETEMP_ID)=minetemp;
                        phi_arr(i,j,k,EEN_ID)=1.5*K_B*phi_arr(i,j,k,eidx)*minetemp;
                    }
                });
            }
        }
    }

    //clean-up
    specdata.clear();
    acoeff.clear();
    bcoeff.clear();
    solution.clear();
    rhs.clear();

    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}


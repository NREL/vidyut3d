#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <Vidyut.H>
#include <Tagging.H>
#include <Chemistry.H>
#include <PlasmaChem.H>
#include <ProbParm.H>
#include <stdio.h>
#include <VarDefines.H>

using namespace amrex;

ProbParm* Vidyut::h_prob_parm = nullptr;
ProbParm* Vidyut::d_prob_parm = nullptr;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
Vidyut::Vidyut()
{
    ReadParameters();
    h_prob_parm = new ProbParm{};
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));
    amrex_probinit(*h_prob_parm, *d_prob_parm);

    plasma_param_names.resize(NUM_PLASMAVARS);
    plasma_param_names[0]="Electron_energy";
    plasma_param_names[1]="Electron_Temp";
    plasma_param_names[2]="Potential";
    plasma_param_names[3]="Efieldx";
    plasma_param_names[4]="Efieldy";
    plasma_param_names[5]="Efieldz";
    plasma_param_names[6]="Electron_Jheat";
    plasma_param_names[7]="Electron_inelasticHeat";
    plasma_param_names[8]="Electron_elasticHeat";
    plasma_param_names[9]="ReducedEF";
    plasma_param_names[10]="SurfaceCharge";
    plasma_param_names[11]="PhotoIon_Src";
    

    allvarnames.resize(NVAR);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        allvarnames[i] = plasmachem::specnames[i];
    }
    for(int i=0;i<NUM_PLASMAVARS;i++)
    {
        allvarnames[i+NUM_SPECIES] = plasma_param_names[i];
    }

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev)
    {
        nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    for(int lev=0;lev<nlevs_max;lev++)
    {
        dt[lev]=fixed_dt;
    }

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    
    ParmParse pp("vidyut");
    pp.queryarr("pot_bc_lo", pot_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("pot_bc_hi", pot_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("eden_bc_lo", eden_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("eden_bc_hi", eden_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("eenrg_bc_lo", eenrg_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("eenrg_bc_hi", eenrg_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("ion_bc_lo", ion_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("ion_bc_hi", ion_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("neutral_bc_lo", neutral_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("neutral_bc_hi", neutral_bc_hi, 0, AMREX_SPACEDIM);

    pp.queryarr("photoion_bc_lo", photoion_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("photoion_bc_hi", photoion_bc_hi, 0, AMREX_SPACEDIM);    

    //foextrap all states as bcs imposed
    //through linear solver
    bcspec.resize(NVAR);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {

        int bctype=(geom[0].isPeriodic(idim))?BCType::int_dir:BCType::foextrap;

        for (int sp=0; sp < NVAR; sp++) 
        {
            bcspec[sp].setLo(idim, bctype);
            bcspec[sp].setHi(idim, bctype);
        }
    }

    // Find the electron index, throw error if not found
    int E_idx = (plasmachem::find_id("E") != -1) ? plasmachem::find_id("E"):
                (plasmachem::find_id("E-") != -1)? plasmachem::find_id("E-"):
                (plasmachem::find_id("e") != -1) ? plasmachem::find_id("e"):
                plasmachem::find_id("e-");

    if(E_idx != -1){
      E_IDX = E_idx;
    } else{
      amrex::Abort("Electron not found in chemistry mechanism!\n");
    }

    if(multicompsolves)
    {
        if(   ion_bc_lo[0] == ROBINBC || ion_bc_hi[0] == ROBINBC 
#if AMREX_SPACEDIM > 1
          ||  ion_bc_lo[1] == ROBINBC || ion_bc_lo[1] == ROBINBC
#if AMREX_SPACEDIM == 3
           || ion_bc_lo[2] == ROBINBC || ion_bc_lo[2] == ROBINBC
#endif
#endif
           )
        {
            amrex::Print()<<"cannot do multicomponent solves with Robin BC for ions**\n";
            amrex::Print()<<"comp_ion_chunks set to 1***\n";
            comp_ion_chunks=1;
        }
        
        if(   neutral_bc_lo[0] == ROBINBC || neutral_bc_hi[0] == ROBINBC 
#if AMREX_SPACEDIM > 1
          ||  neutral_bc_lo[1] == ROBINBC || neutral_bc_lo[1] == ROBINBC
#if AMREX_SPACEDIM == 3
           || neutral_bc_lo[2] == ROBINBC || neutral_bc_lo[2] == ROBINBC
#endif
#endif
           )
        {
            amrex::Print()<<"cannot do multicomponent solves with Robin BC for neutrals**\n";
            amrex::Print()<<"comp_neutral_chunks set to 1***\n";
            comp_neutral_chunks=1;
        }
    }

    //Check inputs for axisymmetric geometry
    //only needed if one boundary is at r=0 and 
    //the user will set the condition accordingly
    /*if(geom[0].IsRZ()){
        if(AMREX_SPACEDIM != 2) amrex::Abort("AMREX_SPACEDIM should be 2 for axisymmetric coordinates");
        // Axisymmetric implementation assumes x-low boundary is the axis of symmatry
        if(pot_bc_lo[0] != HNEUBC || eden_bc_lo[0] != HNEUBC || ion_bc_lo[0] != HNEUBC || neutral_bc_lo[0] != HNEUBC 
           || eenrg_bc_lo[0] != HNEUBC || photoion_bc_lo[0] != HNEUBC)
        {
            if(pot_bc_lo[0] != AXISBC || eden_bc_lo[0] != AXISBC || ion_bc_lo[0] != AXISBC 
               || neutral_bc_lo[0] != AXISBC || eenrg_bc_lo[0] != AXISBC || photoion_bc_lo[0] != AXISBC)
            {
                amrex::Abort("All x_lo boundaries must be Homogenous Neumann (equal to 2) or axis (equal to 5)");
            }
        }
    }*/

}

Vidyut::~Vidyut()
{
    delete h_prob_parm;
    The_Arena()->free(d_prob_parm);
}
// initializes multilevel data
void Vidyut::InitData()
{
    BL_PROFILE("Vidyut::InitData()");
    ProbParm* localprobparm = d_prob_parm;

    if (restart_chkfile == "")
    {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0 || chk_time > 0.0)
        {
            WriteCheckpointFile(0);
        }

    } 
    else
    {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0 || plot_time > 0.0)
    {
        WritePlotFile(amrex::Math::floor
                      (amrex::Real(istep[0])/amrex::Real(plot_int)));
    }
    if(monitor_file_int > 0)
    {
        WriteMonitorFile(0.0);
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void Vidyut::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow)
{
    BL_PROFILE("Vidyut::ErrorEst()");
    static bool first = true;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        ParmParse pp("vidyut");
        if (pp.contains("tagged_vars"))
        {
            int nvars = pp.countval("tagged_vars");
            refine_phi.resize(nvars);
            refine_phigrad.resize(nvars);
            refine_phi_comps.resize(nvars);
            std::string varname;
            for (int i = 0; i < nvars; i++)
            {
                pp.get("tagged_vars", varname, i);
                pp.get((varname + "_refine").c_str(), refine_phi[i]);
                pp.get((varname + "_refinegrad").c_str(), refine_phigrad[i]);

                int spec_id = plasmachem::find_id(varname);
                if (spec_id == -1)
                {
                    int varname_id=-1;
                    auto it=std::find(plasma_param_names.begin(),plasma_param_names.end(),varname);
                    if(it != plasma_param_names.end())
                    {
                        varname_id=it-plasma_param_names.begin();
                    }

                    if (varname_id == -1)
                    {
                        Print() << "Variable name:" << varname << " not found for tagging\n";
                        amrex::Abort("Invalid tagging variable");
                    }
                    else
                    {
                        refine_phi_comps[i] = varname_id+NUM_SPECIES;
                    }
                }
                else
                {
                    refine_phi_comps[i] = spec_id;
                }
            }
        }
    }

    if (refine_phi.size() == 0) return;

    //    const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;

    const MultiFab& state = phi_new[lev];
    MultiFab Sborder(grids[lev], dmap[lev], state.nComp(), 1);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        auto sb_arrays = Sborder.arrays();
        auto tags_arrays = tags.arrays();

        amrex::Real* refine_phi_dat = refine_phi.data();
        amrex::Real* refine_phigrad_dat = refine_phigrad.data();
        int* refine_phi_comps_dat = refine_phi_comps.data();
        int ntagged_comps = refine_phi_comps.size();

        amrex::ParallelFor(Sborder, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            auto statefab = sb_arrays[nbx];
            auto tagfab = tags_arrays[nbx];
            state_based_refinement(i, j, k, tagfab, statefab, refine_phi_dat, refine_phi_comps_dat, ntagged_comps, tagval);
        });

        amrex::ParallelFor(Sborder, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            auto statefab = sb_arrays[nbx];
            auto tagfab = tags_arrays[nbx];
            stategrad_based_refinement(i, j, k, tagfab, statefab, refine_phigrad_dat, refine_phi_comps_dat, ntagged_comps, tagval);
        });
    }
}

// read in some parameters from inputs file
void Vidyut::ReadParameters()
{
    BL_PROFILE("Vidyut::ReadParameters()");
    {
        ParmParse pp; // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.
        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        plot_int_old=plot_int;
        pp.query("plot_int_old", plot_int_old);
        pp.query("plot_time", plot_time);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        chk_int_old=chk_int;
        pp.query("chk_int_old", chk_int_old);
        pp.query("chk_time", chk_time);
        pp.query("restart", restart_chkfile);
    }

    {
        ParmParse pp("vidyut");

        pp.query("dt", fixed_dt);
        pp.query("adaptive_dt", adaptive_dt);
        if(adaptive_dt){
            pp.query("advective_cfl", advective_cfl);
            pp.query("diffusive_cfl", diffusive_cfl);
            pp.query("dielectric_cfl", dielectric_cfl);
            pp.query("dt_min", dt_min);
            pp.query("dt_max", dt_max);
            pp.query("adaptive_dt_delay", adaptive_dt_delay);
            pp.query("dt_stretch", dt_stretch);
        }

        pp.query("linsolve_reltol",linsolve_reltol);
        pp.query("linsolve_abstol",linsolve_abstol);
        pp.query("linsolve_bot_reltol",linsolve_bot_reltol);
        pp.query("linsolve_bot_abstol",linsolve_bot_abstol);

        pp.query("linsolve_num_pre_smooth",linsolve_num_pre_smooth);
        pp.query("linsolve_num_post_smooth",linsolve_num_post_smooth);
        pp.query("linsolve_num_final_smooth",linsolve_num_final_smooth);
        pp.query("linsolve_num_bottom_smooth",linsolve_num_bottom_smooth);

        pp.query("linsolve_maxiter",linsolve_maxiter);
        pp.query("linsolve_max_coarsening_level",linsolve_max_coarsening_level);
        pp.query("bound_specden", bound_specden);
        pp.query("min_species_density",min_species_density);
        pp.query("min_electron_density",min_electron_density);
        pp.query("min_electron_temp",min_electron_temp);
        pp.query("elecenergy_solve",elecenergy_solve);
        pp.query("hyp_order",hyp_order);
        pp.query("do_reactions",do_reactions);
        pp.query("do_transport",do_transport);
        pp.query("do_spacechrg",do_spacechrg);
        pp.query("user_defined_potential", user_defined_potential);
        pp.query("user_defined_species", user_defined_species);
        pp.query("user_defined_vel", user_defined_vel);
        pp.query("do_bg_reactions",do_bg_reactions);
        pp.query("do_photoionization",do_photoionization);
        pp.query("photoion_ID",photoion_ID);
        pp.query("multicompsolves",multicompsolves);
        pp.query("comp_ion_chunks",comp_ion_chunks);
        pp.query("comp_neutral_chunks",comp_neutral_chunks);

        pp.query("gas_temperature",gas_temperature);
        pp.query("gas_pressure",gas_pressure);
        bg_specid_list.resize(0);
        pp.queryarr("bg_species_ids",bg_specid_list);
        
        pp.query("weno_scheme",weno_scheme);
        pp.query("track_surf_charge",track_surf_charge);

        if(hyp_order==1) //first order upwind
        {
            ngrow_for_fillpatch=1;
        }
        else if(hyp_order==2) //second-order flux limited
        {
            ngrow_for_fillpatch=2;
        }
        else if(hyp_order==5)  //weno 5
        {
            ngrow_for_fillpatch=3;
            //amrex::Abort("hyp_order 5 not implemented yet");
        }
        else
        {
            amrex::Abort("Specified hyp_order not implemented yet");
        }

        // Voltage options
        pp.query("voltage_profile", voltage_profile);
        pp.queryarr("voltage_amp_lo", voltage_amp_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("voltage_amp_hi", voltage_amp_hi, 0, AMREX_SPACEDIM);
        if(voltage_profile == 1){
            pp.get("voltage_freq", voltage_freq);
        } else if (voltage_profile == 2) {
            pp.get("voltage_dur", voltage_dur);
            pp.get("voltage_center", voltage_center);
        } 

        pp.query("monitor_file_int", monitor_file_int);
        pp.query("num_timestep_correctors",num_timestep_correctors);
        pp.query("efield_limiter",efield_limiter);

        pp.query("cs_technique",cs_technique);
        if(cs_technique)
        {
            pp.get("cs_ncharges",cs_ncharges);
            cs_locx.resize(cs_ncharges);
            cs_locy.resize(cs_ncharges);
            cs_locz.resize(cs_ncharges);
            cs_pin_locx.resize(cs_ncharges);
            cs_pin_locy.resize(cs_ncharges);
            cs_pin_locz.resize(cs_ncharges);
            cs_voltages.resize(cs_ncharges);

            pp.query("cs_2d",cs_2d);
            pp.getarr("cs_locx",cs_locx);
            pp.getarr("cs_locy",cs_locy);
            pp.getarr("cs_locz",cs_locz);
            pp.getarr("cs_pin_locx",cs_pin_locx);
            pp.getarr("cs_pin_locy",cs_pin_locy);
            pp.getarr("cs_pin_locz",cs_pin_locz);
            pp.getarr("cs_voltages",cs_voltages);
        }
        pp.query("floor_jh",floor_jh);

#ifdef AMREX_USE_HYPRE
        pp.query("use_hypre",use_hypre);
#endif

    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void Vidyut::GetData(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    BL_PROFILE("Vidyut::GetData()");
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    } else
    {
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

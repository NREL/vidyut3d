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
#include <PlasmaChem.H>
#include <AMReX_MLABecLaplacian.H>

// advance solution to final time
void Vidyut::Evolve()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    Real plottime = 0.0;
    Real chktime = 0.0;

    //there is a slight issue when restart file is not a multiple
    //a plot file may get the same number with an "old" file generated
    int plotfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(plot_int));
    int chkfilenum=amrex::Math::floor(amrex::Real(istep[0])/amrex::Real(chk_int));
    if(plot_time > 0.0) plotfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(plot_time));
    if(chk_time > 0.0) chkfilenum=amrex::Math::floor(amrex::Real(cur_time)/amrex::Real(chk_time));
    amrex::Real dt_edrift,dt_ediff,dt_diel_relax;
    amrex::Real dt_edrift_lev,dt_ediff_lev,dt_diel_relax_lev;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..." << std::endl;

        dt_edrift = std::numeric_limits<Real>::max();
        dt_diel_relax = std::numeric_limits<Real>::max();
        dt_ediff = std::numeric_limits<Real>::max();

        for(int lev=0;lev<=finest_level;lev++)
        {
            find_time_scales(lev,dt_edrift_lev,dt_ediff_lev,dt_diel_relax_lev);
            amrex::Print()<<"electron drift, diffusion and dielectric relaxation time scales at level (sec):"<<lev<<"\t"<<
            dt_edrift_lev<<"\t"<<dt_ediff_lev<<"\t"<<dt_diel_relax_lev<<"\n";

            if(dt_edrift_lev < dt_edrift)
            {
                dt_edrift = dt_edrift_lev;
            }
            if(dt_ediff_lev < dt_ediff)
            {
                dt_ediff = dt_ediff_lev;
            }
            if(dt_diel_relax_lev < dt_diel_relax)
            {
                dt_diel_relax = dt_diel_relax_lev;
            }
        }

        amrex::Print()<<"global minimum electron drift, diffusion and dielectric relaxation time scales (sec):"<<
        dt_edrift<<"\t"<<dt_ediff<<"\t"<<dt_diel_relax<<"\n";

        ComputeDt(cur_time, adaptive_dt_delay, dt_edrift, dt_ediff, dt_diel_relax);

        if (max_level > 0 && regrid_int > 0)  // We may need to regrid
        {
            if (istep[0] % regrid_int == 0)
            {
                regrid(0, cur_time);
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << "ADVANCE with time = " << t_new[lev]
            << " dt = " << dt[0] << std::endl;
        }
        amrex::Real dt_common=dt[0]; //no subcycling

        //ngrow fillpatch set in Vidyut.cpp
        //depending on hyperbolic order
        int num_grow=ngrow_for_fillpatch; 

        Vector< Array<MultiFab,AMREX_SPACEDIM> > flux(finest_level+1);

        //face centered solution gradients
        Vector< Array<MultiFab,AMREX_SPACEDIM> > gradne_fc(finest_level+1);
        Vector< Array<MultiFab,AMREX_SPACEDIM> > grad_fc(finest_level+1);

        // Solution and sources MFs
        Vector<MultiFab> expl_src(finest_level+1);
        Vector<MultiFab> rxn_src(finest_level+1);
        Vector<MultiFab> Sborder(finest_level+1);
        Vector<MultiFab> Sborder_old(finest_level+1);
        Vector<MultiFab> phi_tmp(finest_level+1);

        // edge centered efield
        Vector< Array<MultiFab,AMREX_SPACEDIM> > efield_ec(finest_level+1);

        //copy new to old and update time
        for(int lev=0;lev<=finest_level;lev++)
        {
            phi_tmp[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), phi_new[lev].nGrow());
            phi_tmp[lev].setVal(0.0);
            amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            amrex::MultiFab::Copy(phi_tmp[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            t_old[lev] = t_new[lev];
            t_new[lev] += dt_common;
        }

        //allocate flux, expl_src, Sborder
        for(int lev=0;lev<=finest_level;lev++)
        {
            Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder[lev].setVal(0.0);
            
            Sborder_old[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder_old[lev].setVal(0.0);
           
            FillPatch(lev, cur_time, Sborder_old[lev], 0, Sborder_old[lev].nComp());

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                BoxArray ba = grids[lev];
                ba.surroundingNodes(idim);

                flux[lev][idim].define(ba, dmap[lev], 1, 0);
                flux[lev][idim].setVal(0.0);

                efield_ec[lev][idim].define(ba, dmap[lev], 1, 0);
                efield_ec[lev][idim].setVal(0.0);

                gradne_fc[lev][idim].define(ba, dmap[lev], 1, 0);
                gradne_fc[lev][idim].setVal(0.0);

                grad_fc[lev][idim].define(ba, dmap[lev], 1, 0);
                grad_fc[lev][idim].setVal(0.0);
            }
            
            expl_src[lev].define(grids[lev], dmap[lev], 1, 0);
            expl_src[lev].setVal(0.0);

            rxn_src[lev].define(grids[lev], dmap[lev], NUM_SPECIES+1, 0);
            rxn_src[lev].setVal(0.0);
        }
               

        for(int niter=0;niter<num_timestep_correctors;niter++)
        {
            //reset all
            for(int lev=0;lev<=finest_level;lev++)
            {
                Sborder[lev].setVal(0.0);
            
                //grab phi_new all the time
                //at first iter phi new and old are same
                FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    flux[lev][idim].setVal(0.0);
                    gradne_fc[lev][idim].setVal(0.0);
                    grad_fc[lev][idim].setVal(0.0);
                    efield_ec[lev][idim].setVal(0.0);
                }
                expl_src[lev].setVal(0.0);
                rxn_src[lev].setVal(0.0);
            }

            solve_potential(cur_time, Sborder, pot_bc_lo, pot_bc_hi, efield_ec);

            if(cs_technique)
            {
                update_cs_technique_potential(); 
            }

            //fillpatching here to get the latest potentials in 
            //sborder so that it can be used in efield calc
            for(int lev=0;lev<=finest_level;lev++)
            {
                Sborder[lev].setVal(0.0);
                FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
            }

            //update cell-centered electric fields using alternative method if cs_technique is used
            if(cs_technique){
                update_cc_efields(Sborder);
                //fillpatching here to get the latest efields 
                //in sborder so that it can be used in drift vel calcs
                //may be there is a clever way to improve performance 
                for(int lev=0;lev<=finest_level;lev++)
                {
                    Sborder[lev].setVal(0.0);
                    FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
                }
            }

            // Calculate the reactive source terms for all species/levels
            if(do_reactions)
            {
                update_rxnsrc_at_all_levels(Sborder, rxn_src, cur_time);
            }

            //electron density solve
            update_explsrc_at_all_levels(E_IDX, Sborder, flux, rxn_src, expl_src, eden_bc_lo, eden_bc_hi, cur_time);
            implicit_solve_scalar(cur_time,dt_common,E_IDX,Sborder,Sborder_old,expl_src,eden_bc_lo,eden_bc_hi, gradne_fc);

            //electron energy solve
            if(elecenergy_solve)
            {
                update_explsrc_at_all_levels(EEN_ID, Sborder, flux, rxn_src, expl_src, 
                                             eenrg_bc_lo,eenrg_bc_hi,
                                             cur_time);
                
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    compute_elecenergy_source(lev, Sborder[lev],
                                              rxn_src[lev], 
                                              efield_ec[lev],
                                              gradne_fc[lev],
                                              expl_src[lev], cur_time, dt_common,floor_jh);
                }

                implicit_solve_scalar(cur_time,dt_common,EEN_ID, Sborder,Sborder_old, 
                                      expl_src,eenrg_bc_lo,eenrg_bc_hi, grad_fc);
            }

            //all species except electrons solve
            for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
            {
                bool solveflag=true;
                auto it=std::find(bg_specid_list.begin(),bg_specid_list.end(),ind);
                if(it != bg_specid_list.end())
                {
                    solveflag=false;
                }
                if(ind==E_IDX)
                {
                   solveflag=false;
                }

                if(solveflag)
                {
                    //ions
                    if(plasmachem::get_charge(ind)!=0)
                    {
                        update_explsrc_at_all_levels(ind, Sborder, flux, rxn_src,
                                                     expl_src, ion_bc_lo, ion_bc_hi, cur_time);
                        
                        implicit_solve_scalar(cur_time, dt_common, ind, Sborder, Sborder_old,
                                              expl_src,ion_bc_lo,ion_bc_hi, grad_fc);
                    }
                    //neutrals
                    else
                    {
                        update_explsrc_at_all_levels(ind, Sborder, flux, rxn_src, expl_src, 
                                                     neutral_bc_lo, neutral_bc_hi, cur_time);

                        implicit_solve_scalar(cur_time, dt_common, ind, Sborder, Sborder_old,
                                              expl_src,neutral_bc_lo,neutral_bc_hi, grad_fc);
                    }
                } 
                else if (do_bg_reactions)
                {
                    for (int ilev = 0; ilev <= finest_level; ilev++)
                    {
                        amrex::Real minspecden=min_species_density; 
                        int boundspecden = bound_specden;
                        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                        {
                            const Box& bx = mfi.tilebox();
                            Array4<Real> phi_arr = phi_new[ilev].array(mfi);
                            Array4<Real> rxn_arr = rxn_src[ilev].array(mfi);
                            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                phi_arr(i,j,k,ind) += rxn_arr(i,j,k,ind)*dt_common;
                                if(phi_arr(i,j,k,ind) < minspecden && boundspecden)
                                {
                                    phi_arr(i,j,k,ind) = minspecden;
                                }
                            });
                        }
                    }
                }
            }

            if(niter<num_timestep_correctors-1)
            {
                //copy new to old and update time
                for(int lev=0;lev<=finest_level;lev++)
                {
                    amrex::Print()<<"averaging state at iter:"<<niter<<"\n";
                    MultiFab::LinComb(phi_tmp[lev], 0.5, phi_old[lev], 0, 0.5, 
                                      phi_new[lev], 0, 0, phi_new[lev].nComp(), 0);

                    amrex::MultiFab::Copy(phi_new[lev], phi_tmp[lev], 
                                          0, 0, phi_new[lev].nComp(), 0);
                }
            }
            amrex::Print()<<"\n================== timestep iter:"<<niter<<" ================\n";
        }


        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
            ++istep[lev];

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }

        cur_time += dt_common;
        plottime += dt_common;
        chktime += dt_common;

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt_common << std::endl;

        if (plot_time > 0){
            if(plottime > plot_time){
                last_plot_file_step = step + 1;
                plotfilenum++;
                WritePlotFile(plotfilenum);
                plottime = 0.0;
            }
        }
        else if (plot_int > 0 && (step + 1) % plot_int == 0)
        {
            last_plot_file_step = step + 1;
            plotfilenum++;
            WritePlotFile(plotfilenum);
        }

        if(chk_time > 0){
            if(chktime > chk_time){
                chkfilenum++;
                WriteCheckpointFile(chkfilenum);
                chktime = 0.0;
            }
        }
        else if (chk_int > 0 && (step + 1) % chk_int == 0)
        {
            chkfilenum++;
            WriteCheckpointFile(chkfilenum);
        }

        if (monitor_file_int > 0 && (step + 1) % monitor_file_int == 0)
        {
            WriteMonitorFile(cur_time);
        }

        if (cur_time >= stop_time - 1.e-6 * dt_common) break;


        //local cleanup
        flux.clear();
        gradne_fc.clear();
        grad_fc.clear();
        expl_src.clear();
        Sborder.clear();
        Sborder_old.clear();
        phi_tmp.clear();
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step)
    {
        plotfilenum++;
        WritePlotFile(plotfilenum);
    }
}

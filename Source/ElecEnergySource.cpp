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
#include <compute_explicit_flux.H>
#include <AMReX_MLABecLaplacian.H>

void Vidyut::compute_elecenergy_source(int lev, 
                            MultiFab& Sborder, 
                            MultiFab& rxnsrc, 
                            Array<MultiFab,AMREX_SPACEDIM>& efield, 
                            Array<MultiFab,AMREX_SPACEDIM>& gradne, 
                            MultiFab& dsdt,
                            Real time, Real dt, int floor_jh)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int ncomp = Sborder.nComp();
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;
    int consteletrans = const_ele_trans;
    int captured_cs_technique=cs_technique;
    amrex::Real elemob = ele_mob;
    amrex::Real elediff = ele_diff;
    int eidx = E_IDX;
    
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();
        
    GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
    GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);

        Array4<Real> sborder_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);
        Array4<Real> phi_arr = phi_new[lev].array(mfi);
        Array4<Real> rxn_arr = rxnsrc.array(mfi);
        
        GpuArray<Array4<Real>, AMREX_SPACEDIM> 
        gradne_arr{AMREX_D_DECL(gradne[0].array(mfi), 
                                gradne[1].array(mfi), 
                                gradne[2].array(mfi))};

        GpuArray<Array4<Real>, AMREX_SPACEDIM>
        ef_arr{AMREX_D_DECL(efield[0].array(mfi),
                            efield[1].array(mfi),
                            efield[2].array(mfi))};

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    
            //Joule heating
            amrex::Real mu,dcoeff,etemp,ne;
            amrex::Real efield_x,efield_y,efield_z,efield_face,gradne_face;
            amrex::Real charge=plasmachem::get_charge(eidx)*ECHARGE;
            amrex::Real current_density;
            amrex::Real elec_jheat=0.0;

            //FIXME: This can be done more efficiently sweeping over
            //faces
            for(int idim=0;idim<AMREX_SPACEDIM;idim++)
            {
                //left,right,top,bottom,front,back
                for(int f=0;f<2;f++)
                {
                    IntVect face{AMREX_D_DECL(i,j,k)};
                    face[idim]+=f;

                    if(face[idim]==domlo[idim])
                    {
                       //do the interior face
                       face[idim]+=1;
                    }
                    
                    if(face[idim]==(domhi[idim]+1))
                    {
                       //do the interior face
                        face[idim]-=1;
                    }
                    
                    IntVect lcell{AMREX_D_DECL(face[0],face[1],face[2])};
                    IntVect rcell{AMREX_D_DECL(face[0],face[1],face[2])};
                    lcell[idim]-=1;
                    
                    etemp=0.5*(sborder_arr(lcell,ETEMP_ID) 
                               + sborder_arr(rcell,ETEMP_ID));

                    //FIXME:use face centered updated efields here?
                    Real Esum = 0.0;
                    Real efieldvec_face[AMREX_SPACEDIM];
                    efieldvec_face[0]=0.5*(sborder_arr(lcell,EFX_ID) 
                                  + sborder_arr(rcell,EFX_ID));
                    Esum += std::pow(efieldvec_face[0], 2.0);
#if AMREX_SPACEDIM > 1
                    efieldvec_face[1]=0.5*(sborder_arr(lcell,EFY_ID) 
                                  + sborder_arr(rcell,EFY_ID));
                    Esum += std::pow(efieldvec_face[1], 2.0);
#if AMREX_SPACEDIM == 3
                    efieldvec_face[2]=0.5*(sborder_arr(lcell,EFZ_ID) 
                                  + sborder_arr(rcell,EFZ_ID));
                    Esum += std::pow(efieldvec_face[2], 2.0);
#endif
#endif
            
                    amrex::Real efield_mag=std::sqrt(Esum);

                    ne = 0.5*(sborder_arr(lcell,eidx) 
                              + sborder_arr(rcell,eidx));

                    //efield_face=ef_arr[idim](face);
                    efield_face = (captured_cs_technique) ? efieldvec_face[idim] : ef_arr[idim](face);
                    gradne_face=gradne_arr[idim](face);

                    amrex::Real ndens = 0.0;
                    for(int sp=0; sp<NUM_SPECIES; sp++) ndens += 0.5 * (sborder_arr(lcell,sp) + sborder_arr(rcell,sp));

                    mu = (consteletrans) ? elemob/ndens:specMob(eidx, etemp, ndens,
                                               efield_mag,captured_gastemp);

                    dcoeff = (consteletrans) ? elediff/ndens:specDiff(eidx, etemp, ndens,
                                               efield_mag,captured_gastemp);

                    current_density = charge*(mu*ne*efield_face-dcoeff*gradne_face);
                    elec_jheat += current_density*efield_face;
                }
            }

            //why only 0.5? 
            //because left and right contribute to JxEx, 
            //top and bottom contribute to JyEy,
            //front and back contribute to JzEz,
            //so we are averaging each component, so it is 1/2 and not 1/6
            //read Deconinck, T, S. Mahadevan, and L. L. Raja. 
            //"Discretization of the Joule heating term for plasma discharge 
            //fluid models in unstructured meshes." 
            //Journal of computational physics 228.12 (2009): 4435-4443.

            elec_jheat*=0.5;

            //a switch to make sure joule heating is 
            //heating the electrons and not cooling them
            if(floor_jh)
            {
               if(elec_jheat < 0.0)
               {
                  elec_jheat=0.0;
               }
            }
            
            //inelastic term already added through reaction source
            dsdt_arr(i, j, k) += (elec_jheat);
            // TODO: Adjust reactive source calculations to split into elastic/inelastic
            phi_arr(i,j,k,EJH_ID)=elec_jheat;
            phi_arr(i,j,k,EIH_ID)=rxn_arr(i,j,k,EEN_ID); //EEN_ID is same as NUM_SPECIES+1
            // phi_arr(i,j,k,EEH_ID)=elec_elastic_coll_term;
        });
    }
}

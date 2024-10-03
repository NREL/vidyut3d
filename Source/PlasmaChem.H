using namespace amrex;

namespace plasmachem
{
    extern amrex::Vector<std::string> specnames;
    extern AMREX_GPU_DEVICE_MANAGED  int spec_chrg[NUM_SPECIES];

    void init();
    void close();
    int find_id(std::string specname);
    
    AMREX_GPU_HOST_DEVICE int get_charge(int specid);
}

namespace plasmachem_transport
{

    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                      int dir, int sgn,
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
                      amrex::Real Tg,amrex::Real Pg,
                      amrex::Real app_voltage)
    {
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;



        if(sgn == -1)
        { // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = app_voltage;
            bc_arr(ghost_cell) = app_voltage;
        }
        else if(sgn == 1)
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = app_voltage;
            bc_arr(ghost_cell) = app_voltage;
        }
        else
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
        }
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn,
                    int spec_id, Array4<Real> const &phi,
                    Array4<Real> const& bc_arr,
                    Array4<Real> const& robin_a,
                    Array4<Real> const& robin_b,
                    Array4<Real> const& robin_f,
                    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                    GpuArray<Real, AMREX_SPACEDIM> dx,
                    const Real time,
                    ProbParm const& prob_parm,
                    amrex::Real Tg,amrex::Real Pg)
    {
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=Real(sgn);

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        //default
        robin_a(ghost_cell) = 0.0;
        robin_b(ghost_cell) = 1.0;
        robin_f(ghost_cell) = 0.0;

        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cell_int,sp);

        // Electron density and enerygy fluxes assumed to be driven by electron thermal velocity
        if(spec_id==E_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));

            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            amrex::Real flx = ne*cbar/4.0*outward_normal[dir];

            //why -flx/D?
            //because diffusive flux -D dc/dx = flx
            //so dc/dx = -flx/D
            bc_arr(ghost_cell)=(-flx/D);
        }
        else if(spec_id==EEN_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));

            amrex::Real Te=phi(cell_int,ETEMP_ID);

            amrex::Real D=specDiff(spec_id, Te, ndens, efieldmag, Tg);

            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);

            amrex::Real flx = (2.0*K_B*Te)*ne*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
        }

        // Ion fluxes determined by drift velocity
        else if(plasmachem::get_charge(spec_id) != 0)
        {
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));

            amrex::Real mu = specMob(spec_id,phi(cell_int,ETEMP_ID), ndens,
                                     efieldmag,
                                     Tg);

            amrex::Real D = specDiff(spec_id,phi(cell_int,ETEMP_ID), ndens,
                                     efieldmag,
                                     Tg);

            amrex::Real efield=phi(cell_int,EFX_ID+dir);
            amrex::Real ion_drift_flux = mu*phi(cell_int,spec_id)*efield;
            amrex::Real ion_drift_flux_dotn = ion_drift_flux*outward_normal[dir];

            //if driftvel is pointing towards the wall
            //then flux=v n - D dndx = v n
            //else  flux = v n - D dndx = 0
            bc_arr(ghost_cell)=(ion_drift_flux_dotn >= 0.0)?(-ion_drift_flux/D):0.0;
        }
        // Neutral species fluxes determined by thermal velocity
        else
        {
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));

            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);

            amrex::Real cbar=sqrt(8.0*K_B*Tg/PI/(mw(spec_id)*M_AMU));
            amrex::Real flx = phi(i,j,k,spec_id)*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
        }
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_vel(IntVect cellid, int dir,
                            int specid,
                            Array4<Real> const& phi,
                            GpuArray<Real, AMREX_SPACEDIM> dx,
                            const Real time,
                            ProbParm const& prob_parm,
                            amrex::Real Tg,amrex::Real Pg)
    {
        amrex::Real efieldmag=std::sqrt( std::pow(phi(cellid,EFX_ID),2.0)+
                                        std::pow(phi(cellid,EFY_ID),2.0)+
                                        std::pow(phi(cellid,EFZ_ID),2.0) );

        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cellid,sp);

        amrex::Real mu = specMob(specid,phi(cellid,ETEMP_ID),ndens,
                                 efieldmag,Tg);

        return(mu*phi(cellid,EFX_ID+dir)); 
    }
}

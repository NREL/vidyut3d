#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace user_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void get_photoion_acoeff(int i,int j,int k,int sph_id,
               Array4<Real> const& sb_arr,
               Array4<Real> const& acoeff,
               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
               GpuArray<Real, AMREX_SPACEDIM> dx,
               const Real time,
               ProbParm const& prob_parm,
               Real Tg, Real Pg)
    {
    
    }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    void get_photoion_rhs(int i,int j,int k,int sph_id,
               Array4<Real> const& sb_arr,
               Array4<Real> const& rhs_arr,
               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
               GpuArray<Real, AMREX_SPACEDIM> dx,
               const Real time,
               ProbParm const& prob_parm,
               Real Tg, Real Pg)
    {
    
    }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    void bc_ib(IntVect face,
               int dir, int sgn,
               int solved_comp,
               int rhs_comp,
               Array4<Real> const& sb_arr,
               Array4<Real> const& acoeff,
               Array4<Real> const& rhs,
               GpuArray<int, AMREX_SPACEDIM> domlo,
               GpuArray<int, AMREX_SPACEDIM> domhi,
               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
               GpuArray<Real, AMREX_SPACEDIM> dx,
               const Real time,
               ProbParm const& prob_parm,
               Real Tg, Real Pg)
    {

    }
    
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    int is_dielectric(int i, int j, int k, int dir, int sgn, 
                      const GpuArray<Real,AMREX_SPACEDIM> prob_lo, 
                      const GpuArray<Real,AMREX_SPACEDIM> prob_hi, 
                      const GpuArray<Real,AMREX_SPACEDIM> dx,
                      Real time,
                      ProbParm const& prob_parm)
    {
        int val=0;
        if(sgn == 1 && dir==0)
        {
            val=1;
        }
        return(val);
    } 

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

        amrex::Real y=prob_lo[1]+(cell_int[1]+0.5)*dx[1];

        if(sgn == 1 && dir==0) 
        { 
            // hi sides
            if(y>=prob_parm.y1 && y<=(prob_parm.y1+prob_parm.l1))
            {
               robin_a(ghost_cell) = 1.0;
               robin_b(ghost_cell) = 0.0;
               amrex::Real V=prob_parm.rfv1*std::sin(2.0*PI*prob_parm.freq1*time);
               amrex::Real phic=phi(cell_int,POT_ID);
               amrex::Real sig=phi(cell_int,SRFCH_ID);
               amrex::Real delta=dx[0];
               amrex::Real t=prob_parm.dielthick;
               amrex::Real K=prob_parm.Kdiel;
               amrex::Real phib=(phic + sig*delta/EPS0 
                                 + V*K*delta/t)/(K*delta/t+1.0);
               robin_f(ghost_cell)=phib;
            }
            else if(y>=prob_parm.y2 && y<=(prob_parm.y2+prob_parm.l2))
            {
               robin_a(ghost_cell) = 1.0;
               robin_b(ghost_cell) = 0.0;
               amrex::Real V=prob_parm.rfv2*std::sin(2.0*PI*prob_parm.freq2*time);
               amrex::Real phic=phi(cell_int,POT_ID);
               amrex::Real sig=phi(cell_int,SRFCH_ID);
               amrex::Real delta=dx[0];
               amrex::Real t=prob_parm.dielthick;
               amrex::Real K=prob_parm.Kdiel;
               amrex::Real phib=(phic + sig*delta/EPS0 
                                 + V*K*delta/t)/(K*delta/t+1.0);
               robin_f(ghost_cell)=phib;
            
            }
            else
            {
               robin_a(ghost_cell) = 0.0;
               robin_b(ghost_cell) = 1.0;
               //n is outward normal
               //E.n = -sig/eps0
               //-dphidn = -sig/eps0
               //dphidn = sig/eps0
               robin_f(ghost_cell) = phi(cell_int,SRFCH_ID)/EPS0;
               //robin_f(ghost_cell) = 0.0;
            }
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
    void photoionization_bc(int i, int j, int k,
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
                      amrex::Real Tg,amrex::Real Pg)
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

        if(sgn == -1 && dir==1) 
        { // lo sides
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
        }
        else if(sgn == 1 && dir==1)
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
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
                    int spec_id, int bcspec_id,
                    Array4<Real> const &phi, 
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
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;
        
        amrex::Real y=prob_lo[1]+(cell_int[1]+0.5)*dx[1];
        
        //default
        robin_a(ghost_cell,bcspec_id) = 0.0;
        robin_b(ghost_cell,bcspec_id) = 1.0;
        robin_f(ghost_cell,bcspec_id) = 0.0;
        bc_arr(ghost_cell,bcspec_id) = 0.0;

        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cell_int,sp);
        
        if(spec_id==E_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(amrex::Math::powi<2>(phi(cell_int,EFX_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFY_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFZ_ID)));

            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            amrex::Real flx = ne*cbar/4.0*outward_normal[dir];

            //why -flx/D?
            //because diffusive flux -D dc/dx = flx
            //so dc/dx = -flx/D
            bc_arr(ghost_cell,bcspec_id)=(-flx/D);
        }

        if(spec_id==HEp_ID)
        {
            amrex::Real efieldmag=std::sqrt(amrex::Math::powi<2>(phi(cell_int,EFX_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFY_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFZ_ID)));

            amrex::Real mu = specMob(spec_id,phi(cell_int,ETEMP_ID),
                                     ndens, efieldmag,
                                     Tg);

            amrex::Real D = specDiff(spec_id,phi(cell_int,ETEMP_ID),
                                     ndens, efieldmag,
                                     Tg);

            amrex::Real efield=phi(cell_int,EFX_ID+dir);
            amrex::Real ion_drift_flux = mu*phi(cell_int,spec_id)*efield;
            amrex::Real ion_drift_flux_dotn = ion_drift_flux*outward_normal[dir];

            //if driftvel is pointing towards the wall
            //then flux=v n - D dndx = v n
            //else  flux = v n - D dndx = 0 
            bc_arr(ghost_cell,bcspec_id)=(ion_drift_flux_dotn >= 0.0)?(-ion_drift_flux/D):0.0;
        }

        if(spec_id==HEm_ID || spec_id==HE2m_ID)
        {
            amrex::Real efieldmag=std::sqrt(amrex::Math::powi<2>(phi(cell_int,EFX_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFY_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFZ_ID)));


            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);

            amrex::Real cbar=sqrt(8.0*K_B*Tg/PI*imw(spec_id));
            amrex::Real flx = phi(i,j,k,spec_id)*cbar/4.0*outward_normal[dir];
            bc_arr(ghost_cell,bcspec_id)=(-flx/D);
        }

        if(spec_id==EEN_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(amrex::Math::powi<2>(phi(cell_int,EFX_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFY_ID))+
                                            amrex::Math::powi<2>(phi(cell_int,EFZ_ID)));

            amrex::Real Te=phi(cell_int,ETEMP_ID);

            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);

            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);

            amrex::Real flx = (2.0*K_B*Te)*ne*cbar/4.0*outward_normal[dir];
            bc_arr(ghost_cell,bcspec_id)=(-flx/D);
        }

    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_vel(IntVect iv, int dir,
                            int specid,
                            Array4<Real> const& phi,
                            GpuArray<Real, AMREX_SPACEDIM> dx,
                            const Real time,
                            ProbParm const& prob_parm,
                            amrex::Real Tg,amrex::Real Pg)
    {
        amrex::Real efieldmag=std::sqrt( amrex::Math::powi<2>(phi(iv,EFX_ID))+
                                        amrex::Math::powi<2>(phi(iv,EFY_ID))+
                                        amrex::Math::powi<2>(phi(iv,EFZ_ID)) );
        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(iv,sp);

        amrex::Real mu = specMob(specid,phi(iv,ETEMP_ID),
                                 ndens,efieldmag,Tg);

        return(mu*phi(iv,EFX_ID+dir));
    }
}
#endif

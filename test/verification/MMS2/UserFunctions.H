#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>
#include <VarDefines.H>

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
        IntVect lcell=face;
        IntVect rcell=face;
        lcell[dir]-=1;
        amrex::Real xloc;
        xloc=prob_lo[prob_parm.axdir]+face[prob_parm.axdir]*dx[prob_parm.axdir];
        
        //interior cell
        IntVect intcell=(sgn==1)?lcell:rcell;
        
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        amrex::Real dx2=dx[dir]*dx[dir];


        if(solved_comp==POT_ID)
        {
            //We are solving del^2 (phi)=source
            //in vidyut
            //so del(phi)|i+1/2 - del(phi)|i-1/2 = source
            
            //for both left or right being the interior cell 
            //Dirichlet condition works out to the same sign
            amrex::Real phi_dirc=0.0;
            acoeff(intcell)-= 1.0/(0.5*dx2);
            rhs(intcell,rhs_comp)   -= phi_dirc/(0.5*dx2);
        }
        else if(solved_comp==E_ID || solved_comp==HEp_ID || 
                solved_comp==EEN_ID)
        {
            amrex::Real dircval_left=prob_parm.n0;
            amrex::Real alpha=ECHARGE/EPS0;
            amrex::Real dircval_right=1.0/alpha+prob_parm.n0;
            int left_side=(xloc<0.5)?1:0;
            if(solved_comp==HEp_ID)
            {
                dircval_right=1.0/(2.0*alpha)+prob_parm.n0;
            }
            
            amrex::Real ndens = 0.0;
            for(int sp=0; sp<NUM_SPECIES; sp++) ndens += sb_arr(intcell,sp);
            
            amrex::Real emag=std::sqrt(amrex::Math::powi<2>(sb_arr(intcell,EFX_ID))+
                    amrex::Math::powi<2>(sb_arr(intcell,EFY_ID))+
                    amrex::Math::powi<2>(sb_arr(intcell,EFZ_ID)));
        
            amrex::Real Dcoeff=specDiff(solved_comp,sb_arr(intcell,ETEMP_ID),
                    ndens, emag, Tg);

            amrex::Real mu = specMob(solved_comp,sb_arr(intcell,ETEMP_ID), ndens,
                                           emag,
                                           Tg);
            amrex::Real vel=mu*sb_arr(intcell,EFX_ID+dir);
            
            acoeff(intcell) += Dcoeff/(0.5*dx2);

            amrex::Real dircval=(left_side)?dircval_left:dircval_right;

            rhs(intcell,rhs_comp)   += Dcoeff*dircval/(0.5*dx2);
            
            if(sgn==1) //left cell is internal
            {
                amrex::Real upwindval=(vel>0.0)?sb_arr(intcell,solved_comp):dircval;
                rhs(intcell,rhs_comp) -= upwindval*vel/dx[dir];
            }
            else //right cell is internal
            {
                amrex::Real upwindval=(vel>0.0)?dircval:sb_arr(intcell,solved_comp);
                rhs(intcell,rhs_comp) += upwindval*vel/dx[dir];
            }
        }
        else if(solved_comp==HE_ID)
        {
            //Homogenous Neumann
        }
        else
        {
        
        }
    }
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    int is_dielectric(int i, int j, int k, int dir, int sign, 
                      const GpuArray<Real,AMREX_SPACEDIM> prob_lo, 
                      const GpuArray<Real,AMREX_SPACEDIM> prob_hi, 
                      const GpuArray<Real,AMREX_SPACEDIM> dx,
                      Real time,
                      ProbParm const& prob_parm)
    {
        return(0);
    } 

    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                      int dir, int sgn,
                      Array4<Real> const& phi,
                      Array4<Real> const& bcarr,
                      Array4<Real> const& robin_a,
                      Array4<Real> const& robin_b,
                      Array4<Real> const& robin_f,
                      GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                      GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                      GpuArray<Real, AMREX_SPACEDIM> dx,
                      const Real time,
                      ProbParm const& prob_parm,
                      Real Tg, Real Pg, Real app_voltage)
    {
        //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        if(sgn == -1) 
        { // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
        }
        else
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
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
                    Array4<Real> const& bcarr,
                    Array4<Real> const& robin_a,
                    Array4<Real> const& robin_b,
                    Array4<Real> const& robin_f,
                    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                    GpuArray<Real, AMREX_SPACEDIM> dx,
                    const Real time,
                    ProbParm const& prob_parm,
                    Real Tg, Real Pg)

    {
        //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        gcell_adjust= (sgn==-1)?-1:0;
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        if(sgn == -1) 
        { 
            // lo sides
            // default homogenous Neumann
            robin_a(ghost_cell,bcspec_id) = 0.0;
            robin_b(ghost_cell,bcspec_id) = 1.0;
            robin_f(ghost_cell,bcspec_id) = 0.0;

            if(spec_id==E_ID || spec_id==HEp_ID || spec_id==EEN_ID)
            {
                bcarr(ghost_cell,bcspec_id)=0.0+prob_parm.n0;
            }

        }
        else
        {
            robin_a(ghost_cell,bcspec_id) = 0.0;
            robin_b(ghost_cell,bcspec_id) = 1.0;
            robin_f(ghost_cell,bcspec_id) = 0.0;

            amrex::Real alpha=ECHARGE/EPS0;

            if(spec_id==E_ID)
            {
                bcarr(ghost_cell,bcspec_id)=1.0/alpha+prob_parm.n0;
            }

            if(spec_id==HEp_ID)
            {
                bcarr(ghost_cell,bcspec_id)=1.0/(2.0*alpha)+prob_parm.n0;
            }

            if(spec_id==EEN_ID)
            {
                bcarr(ghost_cell,bcspec_id)=1.0/alpha+prob_parm.n0;
            }
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
        amrex::Real vel=0.0;
        if(dir==prob_parm.axdir)
        {
            amrex::Real efieldmag=std::sqrt( amrex::Math::powi<2>(phi(cellid,EFX_ID))+
                                            amrex::Math::powi<2>(phi(cellid,EFY_ID))+
                                            amrex::Math::powi<2>(phi(cellid,EFZ_ID)));

            amrex::Real ndens=1.0; //dummy
            amrex::Real mu = specMob(specid,phi(cellid,ETEMP_ID), ndens,
                                     efieldmag,
                                     Tg);
            vel=mu*phi(cellid,EFX_ID+dir);
        }

        return(vel);
    }
}
#endif

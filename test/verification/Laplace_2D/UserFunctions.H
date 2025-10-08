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

        IntVect lcell=face;
        IntVect rcell=face;
        lcell[dir]-=1;
        amrex::Real xloc[2];
        xloc[0]=prob_lo[0]+(lcell[0]+0.5)*dx[0];
        xloc[1]=prob_lo[1]+(lcell[1]+0.5)*dx[1];
        xloc[dir]+=0.5*dx[dir];
        amrex::Real rad = sqrt(xloc[0]*xloc[0] + xloc[1]*xloc[1]);

        // interior cell
        IntVect intcell=(sgn==1)?lcell:rcell;
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        amrex::Real dx2=dx[dir]*dx[dir];

        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += sb_arr(intcell,sp);

        if(solved_comp==POT_ID)
        {
          // We are solving -B del.(beta * grad(phi))=source
          // in vidyut where B=1 and beta=-1
          // This is unlike the species equations where
          // B=1 and beta > 0
          // FVM for cell i at the boundary
          //(-B beta grad(phi).\hat{n}|internal
          //-B beta grad(phi).\hat{n}|boundary)/h = source
          // now grad(phi).\hat{n} is always (phi_d-phi_i)/(h/2)
          // therefore, diagonal contribution is always +B beta/(0.5*h^2)
          // rhs contribution is always +B beta phi_d/(0.5*h^2)
          // remember B*beta=-1 for Poisson
          // and B*beta=beta>0 for other equations

          // for both left or right being the interior cell
          // Dirichlet condition works out to the same sign
          amrex::Real phi_dirc = 0.0;
          if (0.5*(prob_parm.inner_radius+prob_parm.outer_radius) < rad){
            phi_dirc = prob_parm.outer_voltage;
          }
          else{
            phi_dirc = prob_parm.inner_voltage;
          }
          acoeff(intcell)-= 1.0/(0.5*dx2);
          rhs(intcell,rhs_comp)   -= phi_dirc/(0.5*dx2);
        }
        else
        {
            //zero flux
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

    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn,
                    int spec_id,
                    int bcspec_id,
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
                                        amrex::Math::powi<2>(phi(iv,EFZ_ID)));

        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(iv,sp);

        amrex::Real mu = specMob(specid,phi(iv,ETEMP_ID),
                                 ndens,efieldmag,Tg);

        return(mu*phi(iv,EFX_ID+dir));
    }
}
#endif

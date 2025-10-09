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
          if (0.5*(prob_parm.inner_radius+prob_parm.outer_radius) < rad){
            if (prob_parm.outer_dirichlet) {
              amrex::Real phi_dirc = prob_parm.outer_voltage;
              acoeff(intcell) -= 1.0 / (0.5 * dx2);
              rhs(intcell, rhs_comp) -= phi_dirc / (0.5 * dx2);
            } else {
              amrex::Real dphidr =
                  (prob_parm.outer_voltage - prob_parm.inner_voltage) /
                  std::log(prob_parm.outer_radius / prob_parm.inner_radius) /
                  prob_parm.outer_radius;
              amrex::Real dphidn =
                  dphidr * std::fabs(xloc[dir]) / prob_parm.outer_radius;
              rhs(intcell, rhs_comp) -= dphidn / dx[dir];
            }
          } else {
            if (prob_parm.inner_dirichlet) {
              amrex::Real phi_dirc = prob_parm.inner_voltage;
              acoeff(intcell) -= 1.0 / (0.5 * dx2);
              rhs(intcell, rhs_comp) -= phi_dirc / (0.5 * dx2);
            } else {
              amrex::Real dphidr =
                  (prob_parm.outer_voltage - prob_parm.inner_voltage) /
                  std::log(prob_parm.outer_radius / prob_parm.inner_radius) /
                  prob_parm.inner_radius;
              amrex::Real dphidn =
                  -dphidr * std::fabs(xloc[dir]) / prob_parm.inner_radius;
              rhs(intcell, rhs_comp) -= dphidn / dx[dir];
            }
          }
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

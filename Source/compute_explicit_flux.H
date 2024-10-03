#ifndef _compute_explicit_flux_H_
#define _compute_explicit_flux_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <UnivConstants.H>
#include <VarDefines.H>
#include <PlasmaChem.H>

using namespace amrex;

#define heps 1.0e-30

//===========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real 
gk(amrex::Real wk,amrex::Real gama_k)
{
    amrex::Real gkval=wk*(gama_k+gama_k*gama_k
                          -3.0*gama_k*wk+wk*wk)
    /(gama_k*gama_k+wk*(1-2.0*gama_k));
    return(gkval);
}
//===========================================================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real weno(amrex::Real phi_im2,amrex::Real phi_im1,amrex::Real phi_i,
                 amrex::Real phi_ip1,amrex::Real phi_ip2,int scheme=1)
{

    //see description of all flavors of WENO here
    //Rathan, Samala, and G. Naga Raju. "A modified fifth-order WENO scheme 
    //for hyperbolic conservation laws." Computers & Mathematics with 
    //Applications 75.5 (2018): 1531-1549.

    amrex::Real eps=1e-20;

    amrex::Real phi_iphalf_1 = (1.0/3.0)*phi_im2 - (7.0/6.0)*phi_im1 + (11.0/6.0)*phi_i;
    amrex::Real phi_iphalf_2 = -(1.0/6.0)*phi_im1 + (5.0/6.0)*phi_i + (1.0/3.0)*phi_ip1;
    amrex::Real phi_iphalf_3 = (1.0/3.0)*phi_i + (5.0/6.0)*phi_ip1 - (1.0/6.0)*phi_ip2;

    amrex::Real beta_1 = (13.0/12.0)*std::pow((phi_im2-2.0*phi_im1+phi_i),2.0) + 
    (1.0/4.0)*std::pow((phi_im2-4.0*phi_im1+3.0*phi_i),2.0);

    amrex::Real beta_2 = (13.0/12.0)*std::pow((phi_im1-2.0*phi_i+phi_ip1),2.0) 
    + (1.0/4.0)*std::pow((phi_im1-phi_ip1),2.0);

    amrex::Real beta_3 = (13.0/12.0)*std::pow((phi_i-2.0*phi_ip1+phi_ip2),2.0) 
    + (1.0/4.0)*std::pow((3.0*phi_i-4.0*phi_ip1+phi_ip2),2.0);

    amrex::Real gama_1 = (1.0/10.0);
    amrex::Real gama_2 = (6.0/10.0);
    amrex::Real gama_3 = (3.0/10.0);

    amrex::Real w1,w2,w3;
    amrex::Real alpha_1,alpha_2,alpha_3;

    //this would return to fifth order scheme
    /*alpha_1 = gama_1;
      alpha_2 = gama_2;
      alpha_3 = gama_3;*/

    //WENO JS   
    if(scheme==1)
    { 
        alpha_1 = gama_1/std::pow((beta_1+eps),2.0);
        alpha_2 = gama_2/std::pow((beta_2+eps),2.0);
        alpha_3 = gama_3/std::pow((beta_3+eps),2.0);
    }
    //WENO M
    else if(scheme==2)
    {
        //do JS stuff
        alpha_1 = gama_1/std::pow((beta_1+eps),2.0);
        alpha_2 = gama_2/std::pow((beta_2+eps),2.0);
        alpha_3 = gama_3/std::pow((beta_3+eps),2.0);

        amrex::Real sum_alpha=(alpha_1+alpha_2+alpha_3);

        w1=alpha_1/sum_alpha;
        w2=alpha_2/sum_alpha;
        w3=alpha_3/sum_alpha;

        alpha_1 = gk(w1,gama_1);
        alpha_2 = gk(w2,gama_2);
        alpha_3 = gk(w3,gama_3);

    }
    //WENO Z
    else if(scheme==3)
    {
        amrex::Real tau = std::abs(beta_1 - beta_3);
        alpha_1=gama_1*(1.0+std::pow(tau/(beta_1+eps),2.0));
        alpha_2=gama_2*(1.0+std::pow(tau/(beta_2+eps),2.0));
        alpha_3=gama_3*(1.0+std::pow(tau/(beta_3+eps),2.0));

    }
    else
    {
        amrex::Abort("Unknown WENO scheme\n");
    }

    amrex::Real sum_alpha=(alpha_1+alpha_2+alpha_3);

    w1=alpha_1/sum_alpha;
    w2=alpha_2/sum_alpha;
    w3=alpha_3/sum_alpha;

    amrex::Real phi_iphalf=w1*phi_iphalf_1+w2*phi_iphalf_2+w3*phi_iphalf_3;

    return(phi_iphalf);
}
//====================================================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void weno_reconstruct(amrex::Real uim2,
                      amrex::Real uim1,
                      amrex::Real ui,
                      amrex::Real uip1,
                      amrex::Real uip2,
                      amrex::Real uip3,
                      amrex::Real &umhalf,
                      amrex::Real &uphalf,
                      int weno_scheme)
{
    uphalf=weno(uim2,uim1,ui,uip1,uip2,weno_scheme);
    umhalf = weno(uip3,uip2,uip1,ui,uim1,weno_scheme);
}
//====================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
amrex::Real get_slope_lim(amrex::Real c,amrex::Real diff_left,
                          amrex::Real diff_mid,
                          amrex::Real diff_right)
{
    //MINBEE vesion of limiter for
    //WAF scheme, see page 502, 
    //Toro's Riemann solver textbook
    amrex::Real r,lim;
    r=0.0;
    r=(c>=0.0)?diff_left/(diff_mid+heps):diff_right/(diff_mid+heps);

    lim=0.0;
    if(r<=0.0)
    {
        lim=1.0;
    }
    else if(r>0.0 && r<=1.0)
    {
        //NOTE: something funky can happen
        //if c or advective CFL is greater than 1
        //But we should be safe with 
        //explicit advection :)
        lim=1.0-(1.0-amrex::Math::abs(c))*r;
    }
    else
    {
        lim=amrex::Math::abs(c);
    }
    return(lim);
}
//====================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_secondorder_WAF_flux(amrex::Real vel_L,amrex::Real vel_R,
                                     amrex::Real phi_Lm1,amrex::Real phi_L,
                                     amrex::Real phi_R,amrex::Real phi_Rp1,
                                     amrex::Real dx, amrex::Real dt)
{
    //See TVD Version of WAF schemes
    //Page 501, Toro's Riemann solver texbook

    amrex::Real diff_left=phi_L-phi_Lm1;
    amrex::Real diff_mid=phi_R-phi_L;
    amrex::Real diff_right=phi_Rp1-phi_R;
    amrex::Real vel_mid=0.5*(vel_L+vel_R);
    amrex::Real vel_sgn = std::copysign(1.0, vel_mid);
    amrex::Real c=vel_mid*dt/dx;
    amrex::Real fL=phi_L*vel_L;
    amrex::Real fR=phi_R*vel_R;

    amrex::Real lim=get_slope_lim(c,diff_left,diff_mid,diff_right);

    amrex::Real flux=0.5*(fL+fR)-0.5*vel_sgn*lim*(fR-fL);

    return(flux);
}
//====================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_firstorder_upwind_flux(amrex::Real vel_L,amrex::Real vel_R,
                                       amrex::Real phi_L,amrex::Real phi_R)
{
    amrex::Real fL=phi_L*vel_L;
    amrex::Real fR=phi_R*vel_R;
    amrex::Real vel_mid=0.5*(vel_L+vel_R);
    amrex::Real vel_sgn = std::copysign(1.0, vel_mid);
    amrex::Real flux = fL *  0.5 * (1.0 + vel_sgn) 
    + fR *  0.5 * (1.0 - vel_sgn);
    return(flux);
}
//====================================================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void compute_flux(
    int i,
    int j,
    int k,
    int sweepdir,
    int specid,
    Array4<Real> const& phi, // state variable
    const GpuArray<int, AMREX_SPACEDIM> bclo,
    const GpuArray<int, AMREX_SPACEDIM> bchi,
    const GpuArray<int, AMREX_SPACEDIM> domlo,
    const GpuArray<int, AMREX_SPACEDIM> domhi,
    Array4<Real> const& fx,
    amrex::Real Tg,
    amrex::Real Pg,
    amrex::Real time,
    const GpuArray<Real, AMREX_SPACEDIM>& dx,
    amrex::Real dt,
    ProbParm const& prob_parm,int hyp_order,
    int user_defined_vel,int weno_scheme)
{

    //|--left1--|--left--|--right--|--right1--|
    //|--(i-2)--|--(i-1)--|--i--|--(i+1)--|
    //remember: face i separates (i-1) and i cells

    IntVect face{AMREX_D_DECL(i,j,k)};
    IntVect left{AMREX_D_DECL(i, j, k)};
    IntVect right{AMREX_D_DECL(i, j, k)};
    left[sweepdir] -= 1;

    Real phi_L = phi(left, specid);
    Real phi_R = phi(right, specid);
    Real vel_L;
    Real vel_R;

    if(user_defined_vel == 1){
        vel_L = user_transport::compute_vel(left,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
        vel_R = user_transport::compute_vel(right,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
    } else {
        vel_L = plasmachem_transport::compute_vel(left,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
        vel_R = plasmachem_transport::compute_vel(right,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
    }

    if(hyp_order==1)
    {
        fx(face)=get_firstorder_upwind_flux(vel_L,vel_R,phi_L,phi_R);
    }
    else if(hyp_order==2)
    {
        IntVect left1{AMREX_D_DECL(i, j, k)};
        IntVect right1{AMREX_D_DECL(i, j, k)};
        left1[sweepdir] -= 2;
        right1[sweepdir] += 1;
        amrex::Real phi_Lm1 = phi(left1, specid);
        amrex::Real phi_Rp1 = phi(right1, specid);
        amrex::Real dxmin=amrex::min(AMREX_D_DECL(dx[0],dx[1],dx[2]));

        fx(face)=get_secondorder_WAF_flux(vel_L,vel_R,phi_Lm1,phi_L,
                                          phi_R,phi_Rp1,dxmin,dt);
    }
    else if(hyp_order==5)
    {
        IntVect left1{AMREX_D_DECL(i, j, k)};
        IntVect right1{AMREX_D_DECL(i, j, k)};
        IntVect left2{AMREX_D_DECL(i, j, k)};
        IntVect right2{AMREX_D_DECL(i, j, k)};
        left1[sweepdir] -= 2;
        right1[sweepdir] += 1;
        left2[sweepdir] -= 3;
        right2[sweepdir] += 2;

        amrex::Real vel_Lm1, vel_Lm2, vel_Rp1, vel_Rp2;

        if(user_defined_vel == 1){
            vel_Lm1 = user_transport::compute_vel(left1,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Lm2 = user_transport::compute_vel(left2,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Rp1 = user_transport::compute_vel(right1,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Rp2 = user_transport::compute_vel(right2,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
        } else {
            vel_Lm1 = plasmachem_transport::compute_vel(left1,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Lm2 = plasmachem_transport::compute_vel(left2,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Rp1 = plasmachem_transport::compute_vel(right1,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
            vel_Rp2 = plasmachem_transport::compute_vel(right2,sweepdir,specid,phi,dx,time,prob_parm,Tg,Pg);
        }

        amrex::Real ulm2=phi(left2,specid);
        amrex::Real ulm1=phi(left1,specid);
        amrex::Real ul=phi(left,specid);
        amrex::Real ur=phi(right,specid);
        amrex::Real urp1=phi(right1,specid);
        amrex::Real urp2=phi(right2,specid);

        amrex::Real flm2=ulm2*vel_Lm2;
        amrex::Real flm1=ulm1*vel_Lm1;
        amrex::Real fl=ul*vel_L;
        amrex::Real fr=ur*vel_R;
        amrex::Real frp1=urp1*vel_Rp1;
        amrex::Real frp2=urp2*vel_Rp2;

        amrex::Real fmhalf=0.0;
        amrex::Real fphalf=0.0;
        amrex::Real umhalf=0.0;
        amrex::Real uphalf=0.0;
        amrex::Real lambda_max=amrex::Math::abs(0.5*(vel_L+vel_R));

        weno_reconstruct(flm2,flm1,fl,fr,frp1,frp2,fmhalf,fphalf,weno_scheme);
        weno_reconstruct(ulm2,ulm1,ul,ur,urp1,urp2,umhalf,uphalf,weno_scheme);

        fx(face)=0.5*(fphalf+lambda_max*uphalf) +
        0.5*(fmhalf-lambda_max*umhalf);
    }
    else
    {
        amrex::Abort("hyp_order not implemented yet\n");
    }

    if( (face[sweepdir]==domlo[sweepdir]) && (bclo[sweepdir]==IHNEUBC) )
    {
        fx(face)=0.0;
    }
    if( (face[sweepdir]==(domhi[sweepdir]+1)) && (bchi[sweepdir]==IHNEUBC) )
    {
        fx(face)=0.0;
    }
}
//====================================================================
#endif

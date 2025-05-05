#ifndef _helperfuncs_H_
#define _helperfuncs_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;


AMREX_GPU_HOST_DEVICE 
AMREX_FORCE_INLINE 
amrex::Real get_efield_alongdir(int i,int j,int k,int dir,
          const GpuArray<int,AMREX_SPACEDIM> domlo,
          const GpuArray<int,AMREX_SPACEDIM> domhi, 
          const GpuArray<Real,AMREX_SPACEDIM> dx, 
          Array4<Real>const &s_arr)
{
    IntVect cellid{AMREX_D_DECL(i, j, k)};
    Real efield_dir=0.0;

    if(cellid[dir]>domlo[dir] && cellid[dir]<domhi[dir])
    {
        IntVect ip1{AMREX_D_DECL(i,j,k)};
        IntVect im1{AMREX_D_DECL(i,j,k)};

        ip1[dir]+=1;
        im1[dir]-=1;

        efield_dir=-0.5*(s_arr(ip1,POT_ID)-s_arr(im1,POT_ID))/dx[dir];
    }
    else if(cellid[dir]==domlo[dir])
    {
        IntVect ip1{AMREX_D_DECL(i,j,k)};
        IntVect ip2{AMREX_D_DECL(i,j,k)};
        ip1[dir]+=1;
        ip2[dir]+=2;
    
        //fdash = (-f_{i+2}+4f_{i+1}-3f_i)/(2 dx)
        if(domhi[dir] - domlo[dir] > 2){
            efield_dir=-0.5*(  -s_arr(ip2,POT_ID) 
                      + 4.0*s_arr(ip1,POT_ID)
                      - 3.0*s_arr(cellid,POT_ID))/dx[dir];
        } else {
            efield_dir=-(s_arr(ip1,POT_ID) - s_arr(cellid,POT_ID))/dx[dir];
        }
    }
    else if(cellid[dir]==domhi[dir])
    {
        IntVect im1{AMREX_D_DECL(i,j,k)};
        IntVect im2{AMREX_D_DECL(i,j,k)};
        im1[dir]-=1;
        im2[dir]-=2;
    
        //fdash = (f_{i-2}-4f_{i-1}+3f_i)/(2 dx)
        if(domhi[dir] - domlo[dir] > 2){
            efield_dir=-0.5*(     s_arr(im2,POT_ID)
                        - 4.0*s_arr(im1,POT_ID)
                        + 3.0*s_arr(cellid,POT_ID))/dx[dir];
        } else {
            efield_dir=-(s_arr(cellid,POT_ID) - s_arr(im1,POT_ID))/dx[dir];
        }
    }
    else
    {
        //there shouldnt be another option!
    }
    return(efield_dir);

}

AMREX_GPU_HOST_DEVICE 
AMREX_FORCE_INLINE
amrex::Real get_applied_potential(Real current_time, int domain_end, int vprof, Real v1, Real v2, Real vfreq, Real vdur, Real vcen, int np)
{
    Real voltage;
    Real voltage_amp = (domain_end == -1) ? v1:v2;
    if(vprof == 1) {  // Sinusoidal pulse shape
        voltage = sin(2.0*PI*vfreq*current_time)*voltage_amp;
    } else if (vprof == 2) {    // Triangular pulses
        voltage = 0.0;
        amrex::Real pulse_time_tmp;
        for(int i=0; i<np; i++){
            pulse_time_tmp = vcen + (i)*(1.0/vfreq);
            if(current_time <= pulse_time_tmp) {
                voltage += (pulse_time_tmp - current_time < vdur/2.0) ? (1.0 - (pulse_time_tmp - current_time)/(vdur/2.0))*voltage_amp:0.0;
            } else {
                voltage += (current_time - pulse_time_tmp < vdur/2.0) ? (1.0 - (current_time - pulse_time_tmp)/(vdur/2.0))*voltage_amp:0.0;
            }
        }
    } else if (vprof == 3) {    // Gaussian profile pulses
        voltage = 0.0;
        amrex::Real pulse_time_tmp;
        amrex::Real pulse_sigma = vdur / (2.0 * sqrt(2.0*log(2.0)));
        for(int i=0; i<np; i++){
            pulse_time_tmp = vcen + (i)*(1.0/vfreq);
            voltage += voltage_amp * exp(-0.5 * pow( (current_time - pulse_time_tmp) / pulse_sigma, 2) );
        }
    } else {
        voltage = voltage_amp;
    }

    return voltage;
}

AMREX_GPU_HOST_DEVICE 
AMREX_FORCE_INLINE
amrex::Real get_experimental_potential(Real cur_time, Real dt, Real prev_voltage, std::string exp_volt_name, Real diel_thickness, Real eps_r, Real L_gap, Real vdur, Real pfreq, int np, Real curr_dens)
{
  
  amrex::Real voltage = 0.0;

  // First, figure out if we are in the middle of a pulse
  bool in_pulse = false;
  amrex::Real pulse_gap = 1.0 / pfreq;
  amrex::Real time_tmp = fmod(cur_time, pulse_gap);
  amrex::Real final_time = (np-1)*pulse_gap + vdur;
  if(time_tmp < vdur && cur_time <= final_time) in_pulse = true;

  // Exit if we are between pulses
  if(!in_pulse || time_tmp == 0.0) return voltage;

  // Find where to index into the voltage array 
  int ind_lo = 0;
  int ind_hi = 0;
  if(exp_volt_name == "yang_2022"){
    while(voltages::yang_2022[ind_hi][0] <= time_tmp){
        ind_hi++;
    } 
    ind_lo = ind_hi - 1;

    // Calculate first-order dU_app/dt
    amrex::Real dUappdt = (voltages::yang_2022[ind_hi][1] - voltages::yang_2022[ind_lo][1]) / (voltages::yang_2022[ind_hi][0] - voltages::yang_2022[ind_lo][0]);

    // Calculate dUgap/dt
    amrex::Real dUgapdt = 1.0/(1.0 + 2.0*diel_thickness / (eps_r*L_gap)) * (dUappdt - 2.0*diel_thickness/(eps_r*EPS0)*curr_dens);
    voltage = prev_voltage + dt*(dUgapdt);
    // voltage = prev_voltage + dt*(dUappdt);
  }

  return voltage;
}


AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
amrex::Real minmod_limiter(amrex::Real left_difference,amrex::Real right_difference)
{
    amrex::Real phi=1.0;
    if(amrex::Math::abs(right_difference) > 0.0)
    {
        amrex::Real r = left_difference/right_difference;
        phi = std::max(0.0,std::min(1.0,r));
    }
    else
    {
        phi = 1.0;
    }
    return(phi);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real efieldlimiter(int i,int j,int k, int sweepdir,
                          Array4<Real> const& state_arr)
{
    IntVect iv{AMREX_D_DECL(i,j,k)};
    IntVect ivl{AMREX_D_DECL(i,j,k)};
    IntVect ivr{AMREX_D_DECL(i,j,k)};

    ivl[sweepdir]-=1;
    ivr[sweepdir]+=1;

    amrex::Real ldiff=state_arr(iv,POT_ID)-state_arr(ivl,POT_ID);
    amrex::Real rdiff=state_arr(ivr,POT_ID)-state_arr(iv,POT_ID);

    amrex::Real lim=minmod_limiter(ldiff,rdiff);

    return(lim);
}

#endif

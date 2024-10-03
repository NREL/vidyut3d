#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <UnivConstants.H>
#include<VarDefines.H>

/* Elements
0  E
1  He
*/

// Species
#define E_ID 0
#define HE_ID 1
#define HEp_ID 2
#define HEm_ID 3
#define HE2p_ID 4
#define HE2m_ID 5

#define NUM_ELEMENTS 2
#define NUM_SPECIES 6
#define NUM_IONS 3
#define NUM_REACTIONS 9

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real *awt);
//  MISC
void CKAWT(amrex::Real *awt);
void CKNCF(int *ncf);
void CKSYME_STR(amrex::Vector<std::string> &ename);
void CKSYMS_STR(amrex::Vector<std::string> &kname);
void GET_RMAP(int *_rmap);
void CKINU(const int i, int &nspec, int *ki, int *nu);
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real *x,
            amrex::Real *q_f, amrex::Real *q_r);
void progressRateFR(amrex::Real *q_f, amrex::Real *q_r, amrex::Real *sc,
                    amrex::Real T);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int &mm, int &kk, int &ii,
                                                     int &nfit) {
  mm = 2;
  kk = 6;
  ii = 9;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[6] = {
    1822.8884868472639482, // E
    0.2498374807187924,    // HE
    0.2498717270938452,    // HEp
    0.2498374807187924,    // HEm
    0.1249273013664111,    // HE2p
    0.1249187403593962,    // HE2m
};
#endif
const amrex::Real h_global_imw[6] = {
    1822.8884868472639482, // E
    0.2498374807187924,    // HE
    0.2498717270938452,    // HEp
    0.2498374807187924,    // HEm
    0.1249273013664111,    // HE2p
    0.1249187403593962,    // HE2m
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[6] = {
    0.000549, // E
    4.002602, // HE
    4.002053, // HEp
    4.002602, // HEm
    8.004655, // HE2p
    8.005204, // HE2m
};
#endif
const amrex::Real h_global_mw[6] = {
    0.000549, // E
    4.002602, // HE
    4.002053, // HEp
    4.002602, // HEm
    8.004655, // HE2p
    8.005204, // HE2m
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_imw(amrex::Real *imw_new) {
  imw_new[0] = 1822.8884868472639482; // E
  imw_new[1] = 0.2498374807187924;    // HE
  imw_new[2] = 0.2498717270938452;    // HEp
  imw_new[3] = 0.2498374807187924;    // HEm
  imw_new[4] = 0.1249273013664111;    // HE2p
  imw_new[5] = 0.1249187403593962;    // HE2m
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real imw(const int n) {
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_mw(amrex::Real *mw_new) {
  mw_new[0] = 0.000549; // E
  mw_new[1] = 4.002602; // HE
  mw_new[2] = 4.002053; // HEp
  mw_new[3] = 4.002602; // HEm
  mw_new[4] = 8.004655; // HE2p
  mw_new[5] = 8.005204; // HE2m
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real mw(const int n) {
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real &ru, amrex::Real &ruc, amrex::Real &pa) {
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real *species,
                                                    const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = -7.459783900000000e+02 * invT + 1.423710750000000e+01 -
                 2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 3: HEm
    species[3] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 4: HE2p
    species[4] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 5: HE2m
    species[5] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
  } else {
    // species 0: E
    species[0] = -7.459784500000000e+02 * invT + 1.423710750000000e+01 -
                 2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 3: HEm
    species[3] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 4: HE2p
    species[4] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 5: HE2m
    species[5] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
  }
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT(amrex::Real wt[]) {
  get_mw(wt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specMob(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
    amrex::Real mob=0.0;
    //mobility is scaled by charge

    //reduced mobility obtained from Yuan and Raja,
    //IEEE. Trans. Plasma Sci.,31,4,2003
    //they say these are 393K, may be temperature
    //scaling is necessary
    //amrex::Real Pres_ratio = (P_NTP/Pg);
    //amrex::Real elecmob=-0.1132*Pres_ratio;

    amrex::Real Pg = ndens * K_B * T;
    amrex::Real EbyN=emag/ndens;
    amrex::Real meanE = 1.5*Te/eV;

    //mobility from 
    //Turner, Miles M., et al. "Simulation benchmarks for low-pressure
    //plasmas: Capacitive discharges." 
    //Physics of Plasmas 20.1 (2013): 013507.
    amrex::Real Hepmob=2.69*std::pow((1.0+1.2e-3*std::pow(EbyN,2.0)+4.2e-8*std::pow(EbyN,4.0)),-0.125);
    amrex::Real He2pmob=Hepmob;

    //computed by Taaresh Taneja (U Minnesota) using Turner's cross sections
    amrex::Real elecmob = (-1.0)*std::exp(55.0 + 0.3942*std::log(meanE) + 2.134/meanE
               -0.6433/std::pow(meanE,2.0) + (0.7112e-1)/std::pow(meanE,3.0)) / ndens;

    if(specid==E_ID)
    {
        mob=elecmob;
    }
    if(specid==HEp_ID)
    {
        //mob=0.001482*Pres_ratio;
        mob=Hepmob;
    }
    if(specid==HE2p_ID)
    {
        //mob=0.002403*Pres_ratio;
        mob=He2pmob;
    }
    if(specid==NUM_SPECIES)
    {
        mob=fivebythree*elecmob;
    }

    return(mob);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specDiff(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
    amrex::Real Pg = ndens * K_B * T;
    amrex::Real Pres_ratio = (P_NTP/Pg);
    amrex::Real elecdcoeff=0.1737*(Te/17406.0)*Pres_ratio;
    amrex::Real dcoeff=0.0;
    if(specid==E_ID)
    {
        dcoeff=elecdcoeff;
    }
    if(specid==HEp_ID)
    {
        dcoeff=5.026e-5*Pres_ratio;
        //dcoeff=0.0;
    }
    if(specid==HE2p_ID)
    {
        dcoeff=8.148e-5*Pres_ratio;
    }
    if(specid==HEm_ID)
    {
       dcoeff=4.116e-4*Pres_ratio;
    }
    if(specid==HE2m_ID)
    {
        dcoeff=2.029e-4*Pres_ratio;
    }
    if(specid==NUM_SPECIES)
    {
        dcoeff=fivebythree*elecdcoeff;
    }
        return(dcoeff);
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(amrex::Real *qf, amrex::Real *qr, const amrex::Real *sc,
          const amrex::Real * /*sc_qss*/, const amrex::Real *tc,
          const amrex::Real invT) {

  // reaction 2: HE + E => E + HEp + E
  qf[0] = sc[0] * sc[1];
  qr[0] = 0.0;

  // reaction 3: HEm + E => E + HEp + E
  qf[1] = sc[0] * sc[3];
  qr[1] = 0.0;

  // reaction 4: HE2m + E => E + HE2p + E
  qf[2] = sc[0] * sc[5];
  qr[2] = 0.0;

  // reaction 7: HE + HEm + HE => HE2m + HE
  qf[3] = (sc[1] * sc[1]) * sc[3];
  qr[3] = 0.0;

  // reaction 8: HE + HEp + HE => HE2p + HE
  qf[4] = (sc[1] * sc[1]) * sc[2];
  qr[4] = 0.0;

  // reaction 0: E + HE => E + HEm
  qf[5] = sc[0] * sc[1];
  qr[5] = 0.0;

  // reaction 1: E + HEm => E + HE
  qf[6] = sc[0] * sc[3];
  qr[6] = 0.0;

  // reaction 5: E + HE2p => HE + HEm
  qf[7] = sc[0] * sc[4];
  qr[7] = 0.0;

  // reaction 6: 2 HEm => E + HE + HEp
  qf[8] = (sc[3] * sc[3]);
  qr[8] = 0.0;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 6; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[6];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;

  // reaction 0:  E + HE => E + HEm
  k_f = 139000000;
  qf[5] *= k_f;
  qr[5] *= k_f * exp(-(g_RT[0] - g_RT[0] + g_RT[1] - g_RT[3]));
  // reaction 1:  E + HEm => E + HE
  k_f = 6619000;
  qf[6] *= k_f;
  qr[6] *= k_f * exp(-(g_RT[0] - g_RT[0] - g_RT[1] + g_RT[3]));
  // reaction 2:  HE + E => E + HEp + E
  k_f = 1556000;
  qf[0] *= k_f;
  qr[0] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[2])) * (refCinv);
  // reaction 3:  HEm + E => E + HEp + E
  k_f = 466100000;
  qf[1] *= k_f;
  qr[1] *= k_f * exp(-(-g_RT[0] - g_RT[2] + g_RT[3])) * (refCinv);
  // reaction 4:  HE2m + E => E + HE2p + E
  k_f = 763700;
  qf[2] *= k_f;
  qr[2] *= k_f * exp(-(-g_RT[0] - g_RT[4] + g_RT[5])) * (refCinv);
  // reaction 5:  E + HE2p => HE + HEm
  k_f = 324400000000;
  qf[7] *= k_f;
  qr[7] *= k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[3] + g_RT[4]));
  // reaction 6:  2 HEm => E + HE + HEp
  k_f = 162600000;
  qf[8] *= k_f;
  qr[8] *= k_f * exp(-(-g_RT[0] - g_RT[1] - g_RT[2] + 2.000000 * g_RT[3])) *
           (refCinv);
  // reaction 7:  HE + HEm + HE => HE2m + HE
  k_f = 471.6;
  qf[3] *= k_f;
  qr[3] *= k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[5])) * (refC);
  // reaction 8:  HE + HEp + HE => HE2p + HE
  k_f = 36280;
  qf[4] *= k_f;
  qr[4] *= k_f * exp(-(g_RT[1] + g_RT[2] - g_RT[4])) * (refC);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_ener_exch(amrex::Real qf, amrex::Real qr, const amrex::Real * sc, amrex::Real kf, int rxntype, amrex::Real eexci, int
 elidx, amrex::Real * enerExch, amrex::Real Ue, amrex::Real T, const amrex::Real Te)
{
  // (0)='no electron'
  // (1)='excitation/de-excitation'
  // (2)='ionization'
  // (3)='recombination'
  // (4)='attachment'
  // (5)='elastic exchange'
  // (6)='chemi-ioniz'
  // (7)='detachment'

  amrex::Real eV2J = 1.602176620800000e-19;   // eV -> J
  amrex::Real NA = 6.02214085774e23;          // 1/mol
  amrex::Real kB = 1.380649e-23;              // J/K
  amrex::Real ne = sc[E_ID] * NA;             // 1/m3
  amrex::Real mee = Ue / ne;                  // J

  // Molecular weight used to calculate elastic collision energy exchange
  amrex::Real wt[NUM_SPECIES];
  get_mw(wt);

  // Superelastic factor (keeps Te from dropping below gas temperarture)
  amrex::Real sefact = (tanh((Te - T) / 0.5) + 1.0) / 2.0 ;

  // *enerExch = 0.0;
  if(rxntype == 1 || rxntype == 2){
    *enerExch -= (eexci*eV2J)*NA*(qf - qr) * sefact;
  }
  if(rxntype == 3){
    *enerExch -= mee*NA*(qf - qr);
  }
  // if(rxntype[r] == 4){
  //   *enerExch -= mee*NA*(qf[r] - qr[r]);
  // }
  if(rxntype == 5){
    amrex::Real qq = sc[E_ID] * sc[elidx] * kf;
    *enerExch -= 3.0 * kB * (wt[E_ID] / wt[elidx]) * (Te - T) * qq * NA;
  }
  // if(rxntype[r] == 7){
  //   *enerExch += mee*NA*(qf[r] - qr[r]);
  // }

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, const amrex::Real* sc, const amrex::Real T, const amrex::Real Te, amrex::Real * enerExch)
{
  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 6; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 6; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[6];
  gibbs(g_RT, tc);

  // Precalculating values for electron energy exchange evaluation
  amrex::Real ne = sc[E_ID] * 6.02214085774e23;
  amrex::Real Ue = 1.5 * Te * ne * 1.380649e-23;

  {
    // reaction 2:  HE + E => E + HEp + E
    const amrex::Real k_f = 1556000*std::pow(Te,0.68)*std::exp(-2.854092e5/Te);
    const amrex::Real qf = k_f * (sc[0] * sc[1]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    int rxntype = 2;
    amrex::Real eexci = 24.6;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  {
    // reaction 3:  HEm + E => E + HEp + E
    const amrex::Real k_f = 466100000*std::pow(Te,0.6)*std::exp(-5.546e4/Te);
    const amrex::Real qf = k_f * (sc[0] * sc[3]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    int rxntype = 2;
    amrex::Real eexci = 4.78;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  {
    // reaction 4:  HE2m + E => E + HE2p + E
    const amrex::Real k_f = 763700*std::pow(Te,0.71)*std::exp(-3.945e4/Te);
    const amrex::Real qf = k_f * (sc[0] * sc[5]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    int rxntype = 2;
    amrex::Real eexci = 3.4;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  {
    // reaction 7:  HE + HEm + HE => HE2m + HE
    const amrex::Real k_f = 471.6;
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[3]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
  }

  {
    // reaction 8:  HE + HEp + HE => HE2p + HE
    const amrex::Real k_f = 36280;
    const amrex::Real qf = k_f * ((sc[1] * sc[1]) * sc[2]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
  }

  {
    // reaction 0:  E + HE => E + HEm
    const amrex::Real k_f = 139000000*std::pow(Te,0.31)*std::exp(-2.297e5/Te);
    const amrex::Real qf = k_f * (sc[0] * sc[1]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    int rxntype = 1;
    amrex::Real eexci = 19.8;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  {
    // reaction 1:  E + HEm => E + HE
    const amrex::Real k_f = 6619000*std::pow(Te,0.31);
    const amrex::Real qf = k_f * (sc[0] * sc[3]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    int rxntype = 1;
    amrex::Real eexci = -19.8;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  {
    // reaction 5:  E + HE2p => HE + HEm
    const amrex::Real k_f = 324400000000*std::pow(Te,-0.5);
    const amrex::Real qf = k_f * (sc[0] * sc[4]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
  }

  {
    // reaction 6:  2 HEm => E + HE + HEp
    const amrex::Real k_f = 162600000;
    const amrex::Real qf = k_f * ((sc[3] * sc[3]));
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[3] -= 2.000000 * qdot;
    int rxntype = 1;
    amrex::Real eexci = -15.0;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch) {

  // convert to SI
  for (int id = 0; id < 6; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T, Te, ener_exch);

  // Calculate electron energy loss via elastic collisions
  amrex::Real ndens = 0.0;
  for(int sp=0; sp<NUM_SPECIES; sp++) ndens += C[sp]*N_A;
  amrex::Real emag = EN * ndens * 1.0e-21;
  amrex::Real mu = specMob(E_ID,Te,ndens,emag,T);
  amrex::Real nu = -ECHARGE/ME/mu;
  amrex::Real elec_elastic_coll_term= 3.0/2.0 * K_B * C[0] * N_A * (Te-T) * nu * (2.0*ME/(4.0*M_AMU));
  *ener_exch -= elec_elastic_coll_term;

  // convert to chemkin units
  for (int id = 0; id < 6; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[],
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[6]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 6; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 6; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 6; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[],
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[6]; // temporary storage
  amrex::Real PORT =
      1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 6; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 6; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[],
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[6]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 6; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 6; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[],
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[6];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 0.000549; // E
  XW += x[1] * 4.002602; // HE
  XW += x[2] * 4.002053; // HEp
  XW += x[3] * 4.002602; // HEm
  XW += x[4] * 8.004655; // HE2p
  XW += x[5] * 8.005204; // HE2m
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 6; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 6; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRG(int kcharge[]) {
  kcharge[0] = -1; // E
  kcharge[1] = 0;  // HE
  kcharge[2] = 1;  // HEp
  kcharge[3] = 0;  // HEm
  kcharge[4] = 1;  // HE2p
  kcharge[5] = 0;  // HE2m
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRGMASS(amrex::Real zk[]) {

  int kchrg[6];
  CKCHRG(kchrg);

  for (int id = 0; id < 6; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}
#endif

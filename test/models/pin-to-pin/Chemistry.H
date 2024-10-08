#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  E
1  N
*/

// Species
#define E_ID 0
#define N2_ID 1
#define N2p_ID 2

#define NUM_ELEMENTS 2
#define NUM_SPECIES 3
#define NUM_IONS 2
#define NUM_REACTIONS 1

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
  kk = 3;
  ii = 1;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[3] = {
    1822.8884868472639482, // E
    0.0356964374955379,    // N2
    0.0356964374955379,    // N2p
};
#endif
const amrex::Real h_global_imw[3] = {
    1822.8884868472639482, // E
    0.0356964374955379,    // N2
    0.0356964374955379,    // N2p
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[3] = {
    0.000549,  // E
    28.014000, // N2
    28.014000, // N2p
};
#endif
const amrex::Real h_global_mw[3] = {
    0.000549,  // E
    28.014000, // N2
    28.014000, // N2p
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_imw(amrex::Real *imw_new) {
  imw_new[0] = 1822.8884868472639482; // E
  imw_new[1] = 0.0356964374955379;    // N2
  imw_new[2] = 0.0356964374955379;    // N2p
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
  mw_new[0] = 0.000549;  // E
  mw_new[1] = 28.014000; // N2
  mw_new[2] = 28.014000; // N2p
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
    // species 1: N2
    species[1] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041202000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427000000000e-13 * tc[4];
    // species 2: N2p
    species[2] = -1.020899900000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * tc[0] - 7.041202000000000e-04 * tc[1] +
                 6.605369999999999e-07 * tc[2] - 4.701262500000001e-10 * tc[3] +
                 1.222427000000000e-13 * tc[4];
  } else {
    // species 0: E
    species[0] = -7.459784500000000e+02 * invT + 1.423710750000000e+01 -
                 2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 1: N2
    species[1] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439884000000000e-04 * tc[1] +
                 9.474600000000001e-08 * tc[2] - 8.414198333333333e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
    // species 2: N2p
    species[2] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * tc[0] - 7.439884000000000e-04 * tc[1] +
                 9.474600000000001e-08 * tc[2] - 8.414198333333333e-12 * tc[3] +
                 3.376675500000000e-16 * tc[4];
  }
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT(amrex::Real wt[]) {
  get_mw(wt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specMob(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{ 
    //mobility is scaled by charge
    amrex::Real mob=0.0;
    amrex::Real elecmob = -2.3987*std::pow(emag,-0.26);
    if(specid==E_ID)
    {
        mob=elecmob;
    }
    if(specid==N2p_ID)
    {
        mob=0.0;
    }
    if(specid==NUM_SPECIES)
    {
        mob=(5.0/3.0)*elecmob;
    }
    return(mob);
} 
     
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specDiff(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
    amrex::Real dcoeff=0.0; 
    amrex::Real elecdcoeff  = 4.3628e-3 * std::pow(emag,0.22);
    if(specid==E_ID)
    {   
        dcoeff=elecdcoeff;
    }
    if(specid==N2p_ID)
    {   
        dcoeff=0.0;
    }
    if(specid==NUM_SPECIES)
    {   
        dcoeff=(5.0/3.0)*elecdcoeff;
    }
    
    return(dcoeff);
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(amrex::Real *qf, amrex::Real *qr, const amrex::Real *sc,
          const amrex::Real * /*sc_qss*/, const amrex::Real *tc,
          const amrex::Real invT) {

  // reaction 0: N2 + E => E + N2p + E
  qf[0] = sc[0] * sc[1];
  qr[0] = 0.0;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 3; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;

  // reaction 0:  N2 + E => E + N2p + E
  k_f = 0;
  qf[0] *= k_f;
  qr[0] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[2])) * (refCinv);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real *wdot, const amrex::Real *sc, const amrex::Real T, const amrex::Real EN) {
  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 3; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 3; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // Convert E/N (Td) to |E|
  amrex::Real ndens = mixture * 6.02214076e23;
  amrex::Real efield_mag = EN * ndens * 1e-21;

  {
    // reaction 0:  N2 + E => E + N2p + E
    amrex::Real alpha = (1.1944e6 + 4.3666e26/std::pow(efield_mag,3.0))
    *std::exp(-2.73e7/efield_mag);
    amrex::Real eta = 340.75;
    amrex::Real alpha_bar = alpha-eta;
    amrex::Real mu_e = 2.3987 * std::pow(efield_mag,-0.26);
    const amrex::Real k_f = alpha_bar*mu_e*efield_mag; 
    // const amrex::Real qf = k_f * (sc[0] * sc[1]);
    // const amrex::Real qr = 0.0;
    // const amrex::Real qdot = qf - qr;
    wdot[0] += k_f*sc[E_ID];
    wdot[1] -= k_f*sc[E_ID];
    wdot[2] += k_f*sc[E_ID];
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch) {

  // convert to SI
  for (int id = 0; id < 3; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T, EN);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[]) {
  amrex::Real c[3]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 3; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  amrex::Real EN = 0.0;
  // convert to chemkin units
  productionRate(wdot, c, T, EN);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[]) {
  amrex::Real c[3]; // temporary storage
  amrex::Real PORT =
      1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * PORT;
  }

  amrex::Real EN = 0.0;
  // convert to chemkin units
  productionRate(wdot, c, T, EN);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[]) {
  amrex::Real c[3]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 3; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  amrex::Real EN = 0.0;
  // call productionRate
  productionRate(wdot, c, T, EN);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[]) {
  amrex::Real c[3];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 0.000549;  // E
  XW += x[1] * 28.014000; // N2
  XW += x[2] * 28.014000; // N2p
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * ROW;
  }

  amrex::Real EN = 0.0;
  // convert to chemkin units
  productionRate(wdot, c, T, EN);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRG(int kcharge[]) {
  kcharge[0] = -1; // E
  kcharge[1] = 0;  // N2
  kcharge[2] = 1;  // N2p
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRGMASS(amrex::Real zk[]) {

  int kchrg[3];
  CKCHRG(kchrg);

  for (int id = 0; id < 3; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}
#endif

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

#define NUM_GAS_ELEMENTS 2// Elements in the homogeneous phase
#define NUM_GAS_SPECIES 3// Species in the homogeneous phase
#define NUM_GAS_REACTIONS 1// Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00// mol/cm^2

#define NUM_SURFACE_ELEMENTS 0// Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0// Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0// Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 2

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp 
void atomicWeight(amrex::Real *  awt);
//  MISC 
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int * _rmap);
void CKINU(const int i, int &nspec, int * ki, int * nu);


// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
mm = 2;
kk = 3;
ii = 1;
nfit = -1; // Why do you need this anyway ? 
}

//  inverse molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[3]={
1822.8884868472639482,// E
0.0356964374955379,// N2
0.0356971365293034,// N2+
};
#endif
const amrex::Real h_global_imw[3]={
1822.8884868472639482,// E
0.0356964374955379,// N2
0.0356971365293034,// N2+
};

//  molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[3]={
0.000549,// E
28.014000,// N2
28.013451,// N2+
};
#endif
const amrex::Real h_global_mw[3]={
0.000549,// E
28.014000,// N2
28.013451,// N2+
};

//  inverse molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
imw_new[0] = 1822.8884868472639482;// E
imw_new[1] = 0.0356964374955379;// N2
imw_new[2] = 0.0356971365293034;// N2+
}

//  inverse molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real imw(const int n){
#if AMREX_DEVICE_COMPILE
return global_imw[n];
#else
return h_global_imw[n];
#endif
}
//  molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
mw_new[0] = 0.000549;// E
mw_new[1] = 28.014000;// N2
mw_new[2] = 28.013451;// N2+
}

//  molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mw(const int n){
#if AMREX_DEVICE_COMPILE
return global_mw[n];
#else
return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
 ru  = 8.31446261815324e+07; 
 ruc = 1.98721558317399615845; 
 pa  = 1.01325e+06; 
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with no change across T
// species 0: E
species[0] =
+1.50000000e+00 
;

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 1: N2
species[1] =
+2.29867700e+00 +1.40824040e-03 * T-3.96322200e-06  * T2+5.64151500e-09 * T3-2.44485400e-12 * T4
;
// species 2: N2+
species[2] =
+2.77540711e+00 -2.06459157e-03 * T+4.75752301e-06  * T2-3.15664228e-09 * T3+6.70509973e-13 * T4
;
}
else {
// species 1: N2
species[1] =
+1.92664000e+00 +1.48797680e-03 * T-5.68476000e-07  * T2+1.00970380e-10 * T3-6.75335100e-15 * T4
;
// species 2: N2+
species[2] =
+2.58661363e+00 +2.53071949e-04 * T+1.84778214e-07  * T2-4.55257223e-11 * T3+3.26818029e-15 * T4
;
}
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with no change across T
// species 0: E
species[0] =
+2.50000000e+00 
;

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 1: N2
species[1] =
+3.29867700e+00 +1.40824040e-03 * T-3.96322200e-06 * T2+5.64151500e-09 * T3-2.44485400e-12 * T4
;
// species 2: N2+
species[2] =
+3.77540711e+00 -2.06459157e-03 * T+4.75752301e-06 * T2-3.15664228e-09 * T3+6.70509973e-13 * T4
;
}
else {
// species 1: N2
species[1] =
+2.92664000e+00 +1.48797680e-03 * T-5.68476000e-07 * T2+1.00970380e-10 * T3-6.75335100e-15 * T4
;
// species 2: N2+
species[2] =
+3.58661363e+00 +2.53071949e-04 * T+1.84778214e-07 * T2-4.55257223e-11 * T3+3.26818029e-15 * T4
;
}
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);


// species with no change across T
// species 0: E
species[0] =
-7.453750000000000e+02 * invT+1.422081220000000e+01 -2.500000000000000e+00 * logT
;

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 1: N2
species[1] =
-1.020899900000000e+03 * invT-6.516950000000001e-01 -3.298677000000000e+00 * logT-7.041202000000000e-04 * T+6.605369999999999e-07 * T2-4.701262500000001e-10 * T3+1.222427000000000e-13 * T4
;
// species 2: N2+
species[2] =
+1.804811150000000e+05 * invT+1.082185330000000e+00 -3.775407110000000e+00 * logT+1.032295785000000e-03 * T-7.929205016666667e-07 * T2+2.630535233333333e-10 * T3-3.352549865000000e-14 * T4
;
}
else {
// species 1: N2
species[1] =
-9.227977000000000e+02 * invT-3.053888000000000e+00 -2.926640000000000e+00 * logT-7.439884000000000e-04 * T+9.474600000000001e-08 * T2-8.414198333333333e-12 * T3+3.376675500000000e-16 * T4
;
// species 2: N2+
species[2] =
+1.803909940000000e+05 * invT+4.907721999999999e-01 -3.586613630000000e+00 * logT-1.265359745000000e-04 * T-3.079636900000000e-08 * T2+3.793810191666666e-12 * T3-1.634090145000000e-16 * T4
;
}
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;


// species with no change across T
// species 0: E
species[0] =
+1.50000000e+00 -7.45375000e+02 * invT
;

// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 1: N2
species[1] =
+2.29867700e+00 +7.04120200e-04 * T-1.32107400e-06 * T2+1.41037875e-09 * T3-4.88970800e-13 * T4-1.02089990e+03 * invT
;
// species 2: N2+
species[2] =
+2.77540711e+00 -1.03229579e-03 * T+1.58584100e-06 * T2-7.89160570e-10 * T3+1.34101995e-13 * T4+1.80481115e+05 * invT
;
}
else {
// species 1: N2
species[1] =
+1.92664000e+00 +7.43988400e-04 * T-1.89492000e-07 * T2+2.52425950e-11 * T3-1.35067020e-15 * T4-9.22797700e+02 * invT
;
// species 2: N2+
species[2] =
+2.58661363e+00 +1.26535975e-04 * T+6.15927380e-08 * T2-1.13814306e-11 * T3+6.53636058e-16 * T4+1.80390994e+05 * invT
;
}
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real wt[])
{
get_mw(wt);
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
amrex::Real result = 0; 
amrex::Real cvor[3]; //  temporary storage
cv_R(cvor, T);

// perform dot product
for (int id = 0; id < 3; ++id) {
result += x[id]*cvor[id];
}

cvbl = result * 8.31446261815324e+07;
}

// get internal energy as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(const amrex::Real T, amrex::Real uml[])
{
amrex::Real RT = 8.31446261815324e+07*T; // R*T
speciesInternalEnergy(uml, T);

// convert to chemkin units
for (int id = 0; id < 3; ++id) {
uml[id] *= RT;
}
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real *wdot, const amrex::Real *sc, const amrex::Real T, const amrex::Real EN)
{
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
gibbs(g_RT, T);

  // Convert E/N (Td) to |E|
  amrex::Real ndens = mixture * 6.02214076e23;
  amrex::Real efield_mag = EN * ndens * 1e-21;

{
// reaction 0:  N2 + E => E + N2+ + E

amrex::Real alpha = (1.1944e6 + 4.3666e26/std::pow(efield_mag,3.0))
*std::exp(-2.73e7/efield_mag);
amrex::Real eta = 340.75;
amrex::Real alpha_bar = alpha-eta;    
amrex::Real mu_e = 2.3987 * std::pow(efield_mag,-0.26);
const amrex::Real k_f = alpha_bar*mu_e*efield_mag; 
//const amrex::Real qf = k_f * (sc[0]*sc[1]);
//const amrex::Real qr = 0.0;
//const amrex::Real qdot = qf - qr;
wdot[0] += k_f*sc[E_ID];
wdot[1] -= k_f*sc[E_ID];
wdot[2] += k_f*sc[E_ID];
}

}


// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{

productionRate(wdot, C, T, EN);
}

//  species unit charge number 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
kcharge[0] = -1;// E
kcharge[1] = 0;// N2
kcharge[2] = 1;// N2+
}

//  species charge per unit mass 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

int kchrg[3];
CKCHRG(kchrg);

for (int id = 0; id < 3; ++id) {
zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
}
}

#endif

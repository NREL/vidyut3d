#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <UnivConstants.H>
#include <VarDefines.H>

/* Elements
0  E
1  H
2  N
3  O
*/

// Species
#define E_ID 0
#define N2p_ID 1
#define O2p_ID 2
#define On_ID 3
#define NH3p_ID 4
#define H2Op_ID 5
#define O_1D_ID 6
#define N2_A3_ID 7
#define N2_B3_ID 8
#define N2_C3_ID 9
#define O_ID 10
#define NH2_ID 11
#define NH_ID 12
#define H_ID 13
#define H2_ID 14
#define OH_ID 15
#define H2O_ID 16
#define NH3_ID 17
#define O2_ID 18
#define N2_ID 19

#define NUM_GAS_ELEMENTS 4// Elements in the homogeneous phase
#define NUM_GAS_SPECIES 20// Species in the homogeneous phase
#define NUM_GAS_REACTIONS 15// Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00// mol/cm^2

#define NUM_SURFACE_ELEMENTS 0// Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0// Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0// Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 6

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
mm = 4;
kk = 20;
ii = 15;
nfit = -1; // Why do you need this anyway ? 
}

//  inverse molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[20]={
1822.8884868472639482,// E
0.0356971365293034,// N2+
0.0312524890458026,// O2+
0.0625017634094834,// O-
0.0587183495805802,// NH3+
0.0555109881883343,// H2O+
0.0625039064941559,// O(1D)
0.0356964374955379,// N2(A3)
0.0356964374955379,// N2(B3)
0.0356964374955379,// N2(C3)
0.0625039064941559,// O
0.0624102852150034,// NH2
0.0666000666000666,// NH
0.9920634920634921,// H
0.4960317460317460,// H2
0.0587993179279120,// OH
0.0555092978073827,// H2O
0.0587164582232400,// NH3
0.0312519532470779,// O2
0.0356964374955379,// N2
};
#endif
const amrex::Real h_global_imw[20]={
1822.8884868472639482,// E
0.0356971365293034,// N2+
0.0312524890458026,// O2+
0.0625017634094834,// O-
0.0587183495805802,// NH3+
0.0555109881883343,// H2O+
0.0625039064941559,// O(1D)
0.0356964374955379,// N2(A3)
0.0356964374955379,// N2(B3)
0.0356964374955379,// N2(C3)
0.0625039064941559,// O
0.0624102852150034,// NH2
0.0666000666000666,// NH
0.9920634920634921,// H
0.4960317460317460,// H2
0.0587993179279120,// OH
0.0555092978073827,// H2O
0.0587164582232400,// NH3
0.0312519532470779,// O2
0.0356964374955379,// N2
};

//  molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[20]={
0.000549,// E
28.013451,// N2+
31.997451,// O2+
15.999549,// O-
17.030451,// NH3+
18.014451,// H2O+
15.999000,// O(1D)
28.014000,// N2(A3)
28.014000,// N2(B3)
28.014000,// N2(C3)
15.999000,// O
16.023000,// NH2
15.015000,// NH
1.008000,// H
2.016000,// H2
17.007000,// OH
18.015000,// H2O
17.031000,// NH3
31.998000,// O2
28.014000,// N2
};
#endif
const amrex::Real h_global_mw[20]={
0.000549,// E
28.013451,// N2+
31.997451,// O2+
15.999549,// O-
17.030451,// NH3+
18.014451,// H2O+
15.999000,// O(1D)
28.014000,// N2(A3)
28.014000,// N2(B3)
28.014000,// N2(C3)
15.999000,// O
16.023000,// NH2
15.015000,// NH
1.008000,// H
2.016000,// H2
17.007000,// OH
18.015000,// H2O
17.031000,// NH3
31.998000,// O2
28.014000,// N2
};

//  inverse molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
imw_new[0] = 1822.8884868472639482;// E
imw_new[1] = 0.0356971365293034;// N2+
imw_new[2] = 0.0312524890458026;// O2+
imw_new[3] = 0.0625017634094834;// O-
imw_new[4] = 0.0587183495805802;// NH3+
imw_new[5] = 0.0555109881883343;// H2O+
imw_new[6] = 0.0625039064941559;// O(1D)
imw_new[7] = 0.0356964374955379;// N2(A3)
imw_new[8] = 0.0356964374955379;// N2(B3)
imw_new[9] = 0.0356964374955379;// N2(C3)
imw_new[10] = 0.0625039064941559;// O
imw_new[11] = 0.0624102852150034;// NH2
imw_new[12] = 0.0666000666000666;// NH
imw_new[13] = 0.9920634920634921;// H
imw_new[14] = 0.4960317460317460;// H2
imw_new[15] = 0.0587993179279120;// OH
imw_new[16] = 0.0555092978073827;// H2O
imw_new[17] = 0.0587164582232400;// NH3
imw_new[18] = 0.0312519532470779;// O2
imw_new[19] = 0.0356964374955379;// N2
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
mw_new[1] = 28.013451;// N2+
mw_new[2] = 31.997451;// O2+
mw_new[3] = 15.999549;// O-
mw_new[4] = 17.030451;// NH3+
mw_new[5] = 18.014451;// H2O+
mw_new[6] = 15.999000;// O(1D)
mw_new[7] = 28.014000;// N2(A3)
mw_new[8] = 28.014000;// N2(B3)
mw_new[9] = 28.014000;// N2(C3)
mw_new[10] = 15.999000;// O
mw_new[11] = 16.023000;// NH2
mw_new[12] = 15.015000;// NH
mw_new[13] = 1.008000;// H
mw_new[14] = 2.016000;// H2
mw_new[15] = 17.007000;// OH
mw_new[16] = 18.015000;// H2O
mw_new[17] = 17.031000;// NH3
mw_new[18] = 31.998000;// O2
mw_new[19] = 28.014000;// N2
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
// species 1: N2+
species[1] =
+1.804811150000000e+05 * invT+1.082185330000000e+00 -3.775407110000000e+00 * logT+1.032295785000000e-03 * T-7.929205016666667e-07 * T2+2.630535233333333e-10 * T3-3.352549865000000e-14 * T4
;
// species 2: O2+
species[2] =
+1.397422290000000e+05 * invT+4.811498610999999e+00 -4.610171670000000e+00 * logT+3.179759760000000e-03 * T-2.373760400000000e-06 * T2+1.008316025000000e-09 * T3-1.854784390000000e-13 * T4
;
// species 3: O-
species[3] =
+1.879352840000000e+05 * invT-1.893376760000000e+00 -2.500000000000000e+00 * logT
;
// species 4: NH3+
species[4] =
-6.525488000000000e+03 * invT-5.922786000000000e+00 -2.204352000000000e+00 * logT-5.057380000000000e-03 * T+2.442108333333333e-06 * T2-1.206029166666667e-09 * T3+2.664254500000000e-13 * T4
;
// species 5: H2O+
species[5] =
+1.171220000000000e+05 * invT+5.047672768000000e+00 -4.198640560000000e+00 * logT+1.018217050000000e-03 * T-1.086733685000000e-06 * T2+4.573308850000000e-10 * T3-8.859890850000000e-14 * T4
;
// species 6: O(1D)
species[6] =
+5.199679930000000e+04 * invT-2.116556980000000e+00 -2.500000000000000e+00 * logT+1.120545415000000e-13 * T-8.504059783333332e-17 * T2+4.061373533333333e-20 * T3-8.327925200000000e-24 * T4
;
// species 7: N2(A3)
species[7] =
+7.052112920000000e+04 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 8: N2(B3)
species[8] =
+8.420838440000000e+04 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 9: N2(C3)
species[9] =
+1.268940620000000e+05 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 10: O
species[10] =
+2.912225920000000e+04 * invT+1.116333640000000e+00 -3.168267100000000e+00 * logT+1.639659420000000e-03 * T-1.107177326666667e-06 * T2+5.106721866666666e-10 * T3-1.056329855000000e-13 * T4
;
// species 11: NH2
species[11] =
+2.177228000000000e+04 * invT+3.423820000000002e-01 -3.432493000000000e+00 * logT-1.649770000000000e-03 * T+1.102266666666667e-06 * T2-7.159122500000001e-10 * T3+1.786023500000000e-13 * T4
;
// species 12: NH
species[12] =
+4.185047000000000e+04 * invT+8.325769999999997e-01 -3.339758000000000e+00 * logT-6.265045000000000e-04 * T+5.819410000000001e-07 * T2-3.515676666666667e-10 * T3+7.788090000000001e-14 * T4
;
// species 13: H
species[13] =
+2.547365990000000e+04 * invT+2.946682853000000e+00 -2.500000000000000e+00 * logT-3.526664095000000e-13 * T+3.326532733333333e-16 * T2-1.917346933333333e-19 * T3+4.638661660000000e-23 * T4
;
// species 14: H2
species[14] =
-9.179351730000000e+02 * invT+1.661320882000000e+00 -2.344331120000000e+00 * logT-3.990260375000000e-03 * T+3.246358500000000e-06 * T2-1.679767450000000e-09 * T3+3.688058805000000e-13 * T4
;
// species 15: OH
species[15] =
+3.381538120000000e+03 * invT+4.815738570000000e+00 -4.125305610000000e+00 * logT+1.612724695000000e-03 * T-1.087941151666667e-06 * T2+4.832113691666666e-10 * T3-1.031186895000000e-13 * T4
;
// species 16: H2O
species[16] =
-3.029372670000000e+04 * invT+5.047672768000000e+00 -4.198640560000000e+00 * logT+1.018217050000000e-03 * T-1.086733685000000e-06 * T2+4.573308850000000e-10 * T3-8.859890850000000e-14 * T4
;
// species 17: NH3
species[17] =
-6.525488000000000e+03 * invT-5.922786000000000e+00 -2.204352000000000e+00 * logT-5.057380000000000e-03 * T+2.442108333333333e-06 * T2-1.206029166666667e-09 * T3+2.664254500000000e-13 * T4
;
// species 18: O2
species[18] =
-1.063943560000000e+03 * invT+1.247806300000001e-01 -3.782456360000000e+00 * logT+1.498367080000000e-03 * T-1.641217001666667e-06 * T2+8.067745908333334e-10 * T3-1.621864185000000e-13 * T4
;
// species 19: N2
species[19] =
-1.020899900000000e+03 * invT-6.516950000000001e-01 -3.298677000000000e+00 * logT-7.041202000000000e-04 * T+6.605369999999999e-07 * T2-4.701262500000001e-10 * T3+1.222427000000000e-13 * T4
;
}
else {
// species 1: N2+
species[1] =
+1.803909940000000e+05 * invT+4.907721999999999e-01 -3.586613630000000e+00 * logT-1.265359745000000e-04 * T-3.079636900000000e-08 * T2+3.793810191666666e-12 * T3-1.634090145000000e-16 * T4
;
// species 2: O2+
species[2] =
+1.398768230000000e+05 * invT-2.130505470000000e+00 -3.316759220000000e+00 * logT-5.576112200000001e-04 * T+6.391542600000000e-08 * T2-4.773205725000000e-12 * T3+1.388241905000000e-16 * T4
;
// species 3: O-
species[3] =
+1.879408740000000e+05 * invT-1.988834180000000e+00 -2.485420280000000e+00 * logT-1.284893475000000e-05 * T+2.147222966666667e-09 * T2-1.379379058333333e-13 * T3-5.496667200000000e-18 * T4
;
// species 4: NH3+
species[4] =
-6.493270000000000e+03 * invT-5.010192999999999e+00 -2.461904000000000e+00 * logT-3.029583000000000e-03 * T+3.341628333333333e-07 * T2-2.613335833333333e-11 * T3+9.691585000000000e-16 * T4
;
// species 5: H2O+
species[5] =
+1.171220000000000e+05 * invT-1.932777610000000e+00 -3.033992490000000e+00 * logT-1.088459020000000e-03 * T+2.734541966666666e-08 * T2+8.086832250000000e-12 * T3-8.410049600000000e-16 * T4
;
// species 6: O(1D)
species[6] =
+5.199679930000000e+04 * invT-2.116556980000000e+00 -2.500000000000000e+00 * logT-1.046348745000000e-15 * T+1.871095050000000e-19 * T2-2.106359225000000e-23 * T3+1.100479930000000e-27 * T4
;
// species 7: N2(A3)
species[7] =
+7.064415680000000e+04 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 8: N2(B3)
species[8] =
+8.433141200000000e+04 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 9: N2(C3)
species[9] =
+1.270170890000000e+05 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 10: O
species[10] =
+2.921757910000000e+04 * invT-2.214917859999999e+00 -2.569420780000000e+00 * logT+4.298705685000000e-05 * T-6.991409816666667e-09 * T2+8.348149916666666e-13 * T3-6.141684549999999e-17 * T4
;
// species 11: NH2
species[11] =
+2.191977000000000e+04 * invT-2.816567000000000e+00 -2.961311000000000e+00 * logT-1.466349500000000e-03 * T+1.510600000000000e-07 * T2-1.347714166666667e-11 * T3+6.021000000000000e-16 * T4
;
// species 12: NH
species[12] =
+4.207828000000000e+04 * invT-3.096950000000000e+00 -2.760249000000000e+00 * logT-6.876730000000000e-04 * T+7.419856666666667e-08 * T2-6.410660000000000e-12 * T3+2.508796000000000e-16 * T4
;
// species 13: H
species[13] =
+2.547365990000000e+04 * invT+2.946682924000000e+00 -2.500000010000000e+00 * logT+1.154214865000000e-11 * T-2.692699133333334e-15 * T2+3.945960291666667e-19 * T3-2.490986785000000e-23 * T4
;
// species 14: H2
species[14] =
-9.501589220000000e+02 * invT+6.542302510000000e+00 -3.337279200000000e+00 * logT+2.470123655000000e-05 * T-8.324279633333333e-08 * T2+1.496386616666667e-11 * T3-1.001276880000000e-15 * T4
;
// species 15: OH
species[15] =
+3.718857740000000e+03 * invT-2.836911870000000e+00 -2.864728860000000e+00 * logT-5.282522400000000e-04 * T+4.318045966666667e-08 * T2-2.543488950000000e-12 * T3+6.659793800000000e-17 * T4
;
// species 16: H2O
species[16] =
-3.000429710000000e+04 * invT-1.932777610000000e+00 -3.033992490000000e+00 * logT-1.088459020000000e-03 * T+2.734541966666666e-08 * T2+8.086832250000000e-12 * T3-8.410049600000000e-16 * T4
;
// species 17: NH3
species[17] =
-6.493270000000000e+03 * invT-5.010192999999999e+00 -2.461904000000000e+00 * logT-3.029583000000000e-03 * T+3.341628333333333e-07 * T2-2.613335833333333e-11 * T3+9.691585000000000e-16 * T4
;
// species 18: O2
species[18] =
-1.088457720000000e+03 * invT-2.170693450000000e+00 -3.282537840000000e+00 * logT-7.415437700000000e-04 * T+1.263277781666667e-07 * T2-1.745587958333333e-11 * T3+1.083588970000000e-15 * T4
;
// species 19: N2
species[19] =
-9.227977000000000e+02 * invT-3.053888000000000e+00 -2.926640000000000e+00 * logT-7.439884000000000e-04 * T+9.474600000000001e-08 * T2-8.414198333333333e-12 * T3+3.376675500000000e-16 * T4
;
}
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real wt[])
{
get_mw(wt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real specMob(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
amrex::Real mob = 0.0;
amrex::Real E_cs = 40e-20; //A^2 (cross section)
amrex::Real nu_elec = std::sqrt(8.0*K_B*Te/PI/ME)*ndens*E_cs;
amrex::Real elecmob = -2.3987*std::pow(emag,-0.26);

if (specid == E_ID)
{
    //elecmob = -1.0*std::min(0.3,std::exp(-1.35788649+(-1.92169813e-01)*std::log(Te) +  
    //(5.34119379e+03)/Te + (-6.01715433e+06)/std::pow(Te,2) + (2.46885421e+09)/std::pow(Te,3)));
    mob = elecmob;

    // amrex::Print() << "elecmob = " << mob << std::endl;
    //amrex::Print() << "Te = " << Te << std::endl;
    
}
else
{
    mob = 0.0;
}

if(specid==NUM_SPECIES)
{
    mob=fivebythree*elecmob;
}

return(mob);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real specDiff(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
amrex::Real dcoeff = 0.0;
amrex::Real Te_in_ev=Te/eV;
amrex::Real Tg_in_ev=T/eV;
amrex::Real elecmob=specMob(E_ID,Te,ndens,emag,T);
amrex::Real elecdcoeff  = 4.3628e-3 * std::pow(emag,0.22);

  if(specid==E_ID)
  {
      dcoeff=elecdcoeff;
  }
  else
  {
    dcoeff = 0.0;
  }
  
  if(specid==NUM_SPECIES)
  {
    dcoeff=(5.0/3.0)*elecdcoeff;
  }

return(dcoeff);
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
  if(rxntype == 4){
    *enerExch -= mee*NA*(qf - qr);
  }
  if(rxntype == 5){
    amrex::Real qq = sc[E_ID] * sc[elidx] * kf;
    *enerExch -= 3.0 * kB * (wt[E_ID] / wt[elidx]) * (Te - T) * qq * NA;
  }
  if(rxntype == 7){
    *enerExch += mee*NA*(qf - qr);
  }

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, const amrex::Real * sc, const amrex::Real T, const amrex::Real Te, amrex::Real EN, amrex::Real * enerExch)
{
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
const amrex::Real refC = 101325 / 8.31446 * invT;
const amrex::Real refCinv = 1 / refC;

for (int i = 0; i < 20; ++i) {
wdot[i] = 0.0;
}

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int i = 0; i < 20; ++i) {
mixture += sc[i];
}

// compute the Gibbs free energy
amrex::Real g_RT[20];
gibbs(g_RT, T);

// Precalculating values for electron energy exchange evaluation
amrex::Real ne = sc[E_ID] * 6.02214085774e23;
amrex::Real Ue = 1.5 * Te * ne * 1.380649e-23;

{
// reaction 3:  N2 + E => E + N2+ + E
const amrex::Real k_f = N_A*std::exp(-3.32079762e+01 + (2.65942596e-01)*std::log(Te) + 
                        (-2.63638556e+05)/Te + (-1.07235413e+09)/std::pow(Te,2) + 
                        (9.81595921e+11)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[19]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[1] += qdot;
wdot[19] -= qdot;
int rxntype = 2;
amrex::Real eexci = 15.58;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 4:  O2 + E => 2 O + E
const amrex::Real k_f = N_A*std::exp(-2.32649743e+01 + (-8.53413663e-01)*log(Te) + 
                        (-8.66160060e+04)/Te);
const amrex::Real qf = k_f * (sc[0]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] += 2.000000 * qdot;
wdot[18] -= qdot;
int rxntype = 1;
amrex::Real eexci = 6.0;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 5:  O2 + E => O + O(1D) + E
const amrex::Real k_f = N_A*std::exp(-2.50767376e+01 + (-4.59801638e-01)*std::log(Te) + 
                        (-1.68996936e+05)/Te + (1.67919627e+09)/std::pow(Te,2) + 
                        (-7.58597561e+12)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[6] += qdot;
wdot[10] += qdot;
wdot[18] -= qdot;
int rxntype = 1;
amrex::Real eexci = 8.4;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 6:  O2 + E => E + O2+ + E
const amrex::Real k_f = N_A*std::exp(-3.81089916e+01 + (6.55328174e-01)*log(Te) + 
                        (-2.28204343e+05)/Te);
const amrex::Real qf = k_f * (sc[0]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[2] += qdot;
wdot[18] -= qdot;
int rxntype = 2;
amrex::Real eexci = 12.07;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 9:  H2O + E => H + OH + E
const amrex::Real k_f = N_A*std::exp(-3.10013255e+01 + (-3.87482328e-02)*std::log(Te) + 
                        (-1.48745015e+05)/Te);
const amrex::Real qf = k_f * (sc[0]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[13] += qdot;
wdot[15] += qdot;
wdot[16] -= qdot;
int rxntype = 1;
amrex::Real eexci = 7.0;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 10:  H2O + E => H2 + O + E
const amrex::Real k_f = N_A*std::exp(-3.52569958e+01 + (2.89318843e-01)*std::log(Te) + 
                        (-2.00647811e+05)/Te + (-1.47217769e+09)/std::pow(Te,2) + 
                        (1.0)/std::pow(Te,3));                        
const amrex::Real qf = k_f * (sc[0]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] += qdot;
wdot[14] += qdot;
wdot[16] -= qdot;
int rxntype = 1;
amrex::Real eexci = 9.5;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 11:  H2O + E => E + H2O+ + E
const amrex::Real k_f = N_A*std::exp(-3.43664398e+01 + (3.33931971e-01)*std::log(Te) + 
                        (-2.30244103e+05)/Te + (-2.02477843e+09)/std::pow(Te,2) + 
                        (1.0)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[5] += qdot;
wdot[16] -= qdot;
int rxntype = 2;
amrex::Real eexci = 13.76;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 12:  NH3 + E => H + NH2 + E
const amrex::Real k_f = N_A*std::exp(-2.73349619e+01 + (-2.71647193e-01)*std::log(Te) + 
                        (-7.81463651e+04)/Te + (-5.30850417e+08)/std::pow(Te,2) + 
                        (1.00000000e+00)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[17]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[11] += qdot;
wdot[13] += qdot;
wdot[17] -= qdot;
int rxntype = 1;
amrex::Real eexci = 5.72;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 14:  NH3 + E => E + NH3+ + E
const amrex::Real k_f = N_A*std::exp(-3.39902168e+01 + (3.29753098e-01)*std::log(Te) + 
                        (-2.06638805e+05)/Te + (-5.98397173e+08)/std::pow(Te,2) + 
                        (1.00000000e+00)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[17]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[4] += qdot;
wdot[17] -= qdot;
int rxntype = 2;
amrex::Real eexci = 10.2;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 0:  E + N2 => E + N2(A3)
const amrex::Real k_f = N_A*std::exp((-1.75763409e+01) + (-1.43343753e+00)*std::log(Te) + 
                        (-2.57780053e+05)/Te + (2.85742316e+09)/std::pow(Te,2) + 
                        (-1.36013257e+13)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[19]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[7] += qdot;
wdot[19] -= qdot;
int rxntype = 1;
amrex::Real eexci = 6.17;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 1:  E + N2 => E + N2(B3)
const amrex::Real k_f = N_A*std::exp(-1.67067355e+01 + (-1.32208241e+00)*std::log(Te) + 
                        (-2.17286625e+05)/Te + (2.14505360e+09)/std::pow(Te,2) + 
                        (-9.57567162e+12)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[19]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[8] += qdot;
wdot[19] -= qdot;
int rxntype = 1;
amrex::Real eexci = 7.35;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 2:  E + N2 => E + N2(C3)
const amrex::Real k_f = N_A*std::exp(-1.16413995e+01 + (-1.66771820e+00)*std::log(Te) + 
                        (-2.70327743e+05)/Te + (1.64399908e+09)/std::pow(Te,2) + 
                        (-3.19064748e+12)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[19]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[9] += qdot;
wdot[19] -= qdot;
int rxntype = 1;
amrex::Real eexci = 11.03;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 7:  E + O2 => O + O-
const amrex::Real k_f = N_A*(8.8e-11)*std::exp(-51059.8/Te)*(1e-6);
const amrex::Real qf = k_f * (sc[0]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] += qdot;
wdot[10] += qdot;
wdot[18] -= qdot;
int rxntype = 4;
amrex::Real eexci = 0.0;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 8:  E + H2O => H2 + O-
const amrex::Real k_f = N_A*std::exp(-2.16440764e+01 + (-1.36692302e+00)*std::log(Te) + 
                        (-1.07819832e+05)/Te);
const amrex::Real qf = k_f * (sc[0]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] += qdot;
wdot[14] += qdot;
wdot[16] -= qdot;
int rxntype = 1; // Dissociative attachment
amrex::Real eexci = 9.5;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

{
// reaction 13:  E + NH3 => E + 2 H + NH
const amrex::Real k_f = N_A*std::exp(-2.78603532e+01 + (-3.60698526e-02)*log(Te) + 
                        (-1.96732080e+05)/Te + (-1.28507578e+08)/std::pow(Te,2) + 
                        (1.00000000e+00)/std::pow(Te,3));
const amrex::Real qf = k_f * (sc[0]*sc[17]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[12] += qdot;
wdot[13] += 2.000000 * qdot;
wdot[17] -= qdot;
int rxntype = 1;
amrex::Real eexci = 8.65;
int elidx = 0;
comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, T, Te);
}

// Calculate electron energy loss via elastic collisions
amrex::Real emag = 0.0;
amrex::Real ndens = 0.0;
for(int sp=0; sp<NUM_SPECIES; sp++) ndens += sc[sp]*N_A;
amrex::Real mu = specMob(E_ID,Te,ndens,emag,T);
amrex::Real nu = -ECHARGE/ME/mu;
amrex::Real elec_elastic_coll_term= 3.0/2.0 * K_B * sc[0] * N_A * (Te-T) * nu * (2.0*ME/(39.948*M_AMU));
*enerExch -= elec_elastic_coll_term;

}


// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{

// convert to SI
for (int id = 0; id < 20; ++id) {
C[id] *= 1.0e6;
}

// convert to chemkin units
productionRate(wdot, C, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 20; ++id) {
C[id] *= 1.0e-6;
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[20]; // temporary storage
amrex::Real YOW = 0; 
amrex::Real PWORT; 

// Compute inverse of mean molecular wt first
for (int i = 0; i < 20; i++)
{
YOW += y[i]*imw(i);
}
// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// multiply by 1e6 so c goes to SI
PWORT *= 1e6; 
// Now compute conversion (and go to SI)
for (int i = 0; i < 20; i++)
{
c[i] = PWORT * y[i]*imw(i);
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 20; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[20]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 20; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 20; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[20]; // temporary storage

// See Eq 8 with an extra 1e6 so c goes to SI
for (int i = 0; i < 20; i++)
{
c[i] = 1e6 * rho * y[i]*imw(i);
}

// call productionRate
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 20; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho, const amrex::Real T, const amrex::Real x[],  amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[20]; // temporary storage
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*0.000549; // E
XW += x[1]*28.013451; // N2+
XW += x[2]*31.997451; // O2+
XW += x[3]*15.999549; // O-
XW += x[4]*17.030451; // NH3+
XW += x[5]*18.014451; // H2O+
XW += x[6]*15.999000; // O(1D)
XW += x[7]*28.014000; // N2(A3)
XW += x[8]*28.014000; // N2(B3)
XW += x[9]*28.014000; // N2(C3)
XW += x[10]*15.999000; // O
XW += x[11]*16.023000; // NH2
XW += x[12]*15.015000; // NH
XW += x[13]*1.008000; // H
XW += x[14]*2.016000; // H2
XW += x[15]*17.007000; // OH
XW += x[16]*18.015000; // H2O
XW += x[17]*17.031000; // NH3
XW += x[18]*31.998000; // O2
XW += x[19]*28.014000; // N2
// Extra 1e6 factor to take c to SI
ROW = 1e6*rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 20; ++id) {
c[id] = x[id]*ROW;
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 20; ++id) {
wdot[id] *= 1.0e-6;
}
}

//  species unit charge number 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
kcharge[0] = -1;// E
kcharge[1] = 1;// N2+
kcharge[2] = 1;// O2+
kcharge[3] = -1;// O-
kcharge[4] = 1;// NH3+
kcharge[5] = 1;// H2O+
kcharge[6] = 0;// O(1D)
kcharge[7] = 0;// N2(A3)
kcharge[8] = 0;// N2(B3)
kcharge[9] = 0;// N2(C3)
kcharge[10] = 0;// O
kcharge[11] = 0;// NH2
kcharge[12] = 0;// NH
kcharge[13] = 0;// H
kcharge[14] = 0;// H2
kcharge[15] = 0;// OH
kcharge[16] = 0;// H2O
kcharge[17] = 0;// NH3
kcharge[18] = 0;// O2
kcharge[19] = 0;// N2
}

//  species charge per unit mass 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

int kchrg[20];
CKCHRG(kchrg);

for (int id = 0; id < 20; ++id) {
zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
}
}

#endif

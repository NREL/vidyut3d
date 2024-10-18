#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  E
1  H
2  N
3  O
*/

// Species
#define E_ID 0
#define Op_ID 1
#define O2p_ID 2
#define O4p_ID 3
#define On_ID 4
#define O2n_ID 5
#define N4p_ID 6
#define N2p_ID 7
#define NH3p_ID 8
#define H2Op_ID 9
#define O_1D_ID 10
#define O_ID 11
#define N_ID 12
#define N2_A1_ID 13
#define N2_B3_ID 14
#define N2_C3_ID 15
#define NH2_ID 16
#define NH_ID 17
#define H_ID 18
#define H2_ID 19
#define OH_ID 20
#define H2NO_ID 21
#define NOH_ID 22
#define H2O_ID 23
#define NH3_ID 24
#define O2_ID 25
#define N2_ID 26

#define NUM_GAS_ELEMENTS 4// Elements in the homogeneous phase
#define NUM_GAS_SPECIES 27// Species in the homogeneous phase
#define NUM_GAS_REACTIONS 81// Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00// mol/cm^2

#define NUM_SURFACE_ELEMENTS 0// Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0// Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0// Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 10

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
kk = 27;
ii = 81;
nfit = -1; // Why do you need this anyway ? 
}

//  inverse molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[27]={
1822.8884868472639482,// E
0.0625060497257993,// O+
0.0312524890458026,// O2+
0.0156261105720719,// O4+
0.0625017634094834,// O-
0.0312514174667246,// O2-
0.0178483935044992,// N4+
0.0356971365293034,// N2+
0.0587183495805802,// NH3+
0.0555109881883343,// H2O+
0.0625039064941559,// O(1D)
0.0625039064941559,// O
0.0713928749910759,// N
0.0356964374955379,// N2(A1)
0.0356964374955379,// N2(B3)
0.0356964374955379,// N2(C3)
0.0624102852150034,// NH2
0.0666000666000666,// NH
0.9920634920634921,// H
0.4960317460317460,// H2
0.0587993179279120,// OH
0.0312285303853601,// H2NO
0.0322435029341588,// NOH
0.0555092978073827,// H2O
0.0587164582232400,// NH3
0.0312519532470779,// O2
0.0356964374955379,// N2
};
#endif
const amrex::Real h_global_imw[27]={
1822.8884868472639482,// E
0.0625060497257993,// O+
0.0312524890458026,// O2+
0.0156261105720719,// O4+
0.0625017634094834,// O-
0.0312514174667246,// O2-
0.0178483935044992,// N4+
0.0356971365293034,// N2+
0.0587183495805802,// NH3+
0.0555109881883343,// H2O+
0.0625039064941559,// O(1D)
0.0625039064941559,// O
0.0713928749910759,// N
0.0356964374955379,// N2(A1)
0.0356964374955379,// N2(B3)
0.0356964374955379,// N2(C3)
0.0624102852150034,// NH2
0.0666000666000666,// NH
0.9920634920634921,// H
0.4960317460317460,// H2
0.0587993179279120,// OH
0.0312285303853601,// H2NO
0.0322435029341588,// NOH
0.0555092978073827,// H2O
0.0587164582232400,// NH3
0.0312519532470779,// O2
0.0356964374955379,// N2
};

//  molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[27]={
0.000549,// E
15.998451,// O+
31.997451,// O2+
63.995451,// O4+
15.999549,// O-
31.998549,// O2-
56.027451,// N4+
28.013451,// N2+
17.030451,// NH3+
18.014451,// H2O+
15.999000,// O(1D)
15.999000,// O
14.007000,// N
28.014000,// N2(A1)
28.014000,// N2(B3)
28.014000,// N2(C3)
16.023000,// NH2
15.015000,// NH
1.008000,// H
2.016000,// H2
17.007000,// OH
32.022000,// H2NO
31.014000,// NOH
18.015000,// H2O
17.031000,// NH3
31.998000,// O2
28.014000,// N2
};
#endif
const amrex::Real h_global_mw[27]={
0.000549,// E
15.998451,// O+
31.997451,// O2+
63.995451,// O4+
15.999549,// O-
31.998549,// O2-
56.027451,// N4+
28.013451,// N2+
17.030451,// NH3+
18.014451,// H2O+
15.999000,// O(1D)
15.999000,// O
14.007000,// N
28.014000,// N2(A1)
28.014000,// N2(B3)
28.014000,// N2(C3)
16.023000,// NH2
15.015000,// NH
1.008000,// H
2.016000,// H2
17.007000,// OH
32.022000,// H2NO
31.014000,// NOH
18.015000,// H2O
17.031000,// NH3
31.998000,// O2
28.014000,// N2
};

//  inverse molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
imw_new[0] = 1822.8884868472639482;// E
imw_new[1] = 0.0625060497257993;// O+
imw_new[2] = 0.0312524890458026;// O2+
imw_new[3] = 0.0156261105720719;// O4+
imw_new[4] = 0.0625017634094834;// O-
imw_new[5] = 0.0312514174667246;// O2-
imw_new[6] = 0.0178483935044992;// N4+
imw_new[7] = 0.0356971365293034;// N2+
imw_new[8] = 0.0587183495805802;// NH3+
imw_new[9] = 0.0555109881883343;// H2O+
imw_new[10] = 0.0625039064941559;// O(1D)
imw_new[11] = 0.0625039064941559;// O
imw_new[12] = 0.0713928749910759;// N
imw_new[13] = 0.0356964374955379;// N2(A1)
imw_new[14] = 0.0356964374955379;// N2(B3)
imw_new[15] = 0.0356964374955379;// N2(C3)
imw_new[16] = 0.0624102852150034;// NH2
imw_new[17] = 0.0666000666000666;// NH
imw_new[18] = 0.9920634920634921;// H
imw_new[19] = 0.4960317460317460;// H2
imw_new[20] = 0.0587993179279120;// OH
imw_new[21] = 0.0312285303853601;// H2NO
imw_new[22] = 0.0322435029341588;// NOH
imw_new[23] = 0.0555092978073827;// H2O
imw_new[24] = 0.0587164582232400;// NH3
imw_new[25] = 0.0312519532470779;// O2
imw_new[26] = 0.0356964374955379;// N2
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
mw_new[1] = 15.998451;// O+
mw_new[2] = 31.997451;// O2+
mw_new[3] = 63.995451;// O4+
mw_new[4] = 15.999549;// O-
mw_new[5] = 31.998549;// O2-
mw_new[6] = 56.027451;// N4+
mw_new[7] = 28.013451;// N2+
mw_new[8] = 17.030451;// NH3+
mw_new[9] = 18.014451;// H2O+
mw_new[10] = 15.999000;// O(1D)
mw_new[11] = 15.999000;// O
mw_new[12] = 14.007000;// N
mw_new[13] = 28.014000;// N2(A1)
mw_new[14] = 28.014000;// N2(B3)
mw_new[15] = 28.014000;// N2(C3)
mw_new[16] = 16.023000;// NH2
mw_new[17] = 15.015000;// NH
mw_new[18] = 1.008000;// H
mw_new[19] = 2.016000;// H2
mw_new[20] = 17.007000;// OH
mw_new[21] = 32.022000;// H2NO
mw_new[22] = 31.014000;// NOH
mw_new[23] = 18.015000;// H2O
mw_new[24] = 17.031000;// NH3
mw_new[25] = 31.998000;// O2
mw_new[26] = 28.014000;// N2
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
// species 1: O+
species[1] =
+1.879352840000000e+05 * invT-1.893376760000000e+00 -2.500000000000000e+00 * logT
;
// species 2: O2+
species[2] =
+1.397422290000000e+05 * invT+4.811498610999999e+00 -4.610171670000000e+00 * logT+3.179759760000000e-03 * T-2.373760400000000e-06 * T2+1.008316025000000e-09 * T3-1.854784390000000e-13 * T4
;
// species 3: O4+
species[3] =
+1.725300970000000e+05 * invT-1.699199678000000e+01 -1.151078020000000e+00 * logT-1.077851640000000e-02 * T+3.504515650000000e-06 * T2-6.842840075000000e-10 * T3+3.248804900000000e-14 * T4
;
// species 4: O-
species[4] =
+1.879352840000000e+05 * invT-1.893376760000000e+00 -2.500000000000000e+00 * logT
;
// species 5: O2-
species[5] =
+1.397422290000000e+05 * invT+4.811498610999999e+00 -4.610171670000000e+00 * logT+3.179759760000000e-03 * T-2.373760400000000e-06 * T2+1.008316025000000e-09 * T3-1.854784390000000e-13 * T4
;
// species 6: N4+
species[6] =
+2.155531830000000e+05 * invT-4.496941250000000e+00 -3.325965150000000e+00 * logT-7.195459050000000e-03 * T+2.613287550000000e-06 * T2-7.935732066666667e-10 * T3+1.267833240000000e-13 * T4
;
// species 7: N2+
species[7] =
+1.804811150000000e+05 * invT+1.082185330000000e+00 -3.775407110000000e+00 * logT+1.032295785000000e-03 * T-7.929205016666667e-07 * T2+2.630535233333333e-10 * T3-3.352549865000000e-14 * T4
;
// species 8: NH3+
species[8] =
-6.525488000000000e+03 * invT-5.922786000000000e+00 -2.204352000000000e+00 * logT-5.057380000000000e-03 * T+2.442108333333333e-06 * T2-1.206029166666667e-09 * T3+2.664254500000000e-13 * T4
;
// species 9: H2O+
species[9] =
+1.171220000000000e+05 * invT+5.047672768000000e+00 -4.198640560000000e+00 * logT+1.018217050000000e-03 * T-1.086733685000000e-06 * T2+4.573308850000000e-10 * T3-8.859890850000000e-14 * T4
;
// species 10: O(1D)
species[10] =
+5.199679930000000e+04 * invT-2.116556980000000e+00 -2.500000000000000e+00 * logT+1.120545415000000e-13 * T-8.504059783333332e-17 * T2+4.061373533333333e-20 * T3-8.327925200000000e-24 * T4
;
// species 11: O
species[11] =
+2.912225920000000e+04 * invT+1.116333640000000e+00 -3.168267100000000e+00 * logT+1.639659420000000e-03 * T-1.107177326666667e-06 * T2+5.106721866666666e-10 * T3-1.056329855000000e-13 * T4
;
// species 12: N
species[12] =
+5.609890000000000e+04 * invT-1.664495000000000e+00 -2.503071000000000e+00 * logT+1.090009000000000e-05 * T-9.034215000000001e-09 * T2+4.706300000000000e-12 * T3-1.049952000000000e-15 * T4
;
// species 13: N2(A1)
species[13] =
+7.052112920000000e+04 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 14: N2(B3)
species[14] =
+8.420838440000000e+04 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 15: N2(C3)
species[15] =
+1.268940620000000e+05 * invT+5.635349000000001e-01 -3.531005280000000e+00 * logT+6.183049400000000e-05 * T+8.383323883333334e-08 * T2-2.029421766666667e-10 * T3+7.044061750000001e-14 * T4
;
// species 16: NH2
species[16] =
+2.177228000000000e+04 * invT+3.423820000000002e-01 -3.432493000000000e+00 * logT-1.649770000000000e-03 * T+1.102266666666667e-06 * T2-7.159122500000001e-10 * T3+1.786023500000000e-13 * T4
;
// species 17: NH
species[17] =
+4.185047000000000e+04 * invT+8.325769999999997e-01 -3.339758000000000e+00 * logT-6.265045000000000e-04 * T+5.819410000000001e-07 * T2-3.515676666666667e-10 * T3+7.788090000000001e-14 * T4
;
// species 18: H
species[18] =
+2.547365990000000e+04 * invT+2.946682853000000e+00 -2.500000000000000e+00 * logT-3.526664095000000e-13 * T+3.326532733333333e-16 * T2-1.917346933333333e-19 * T3+4.638661660000000e-23 * T4
;
// species 19: H2
species[19] =
-9.179351730000000e+02 * invT+1.661320882000000e+00 -2.344331120000000e+00 * logT-3.990260375000000e-03 * T+3.246358500000000e-06 * T2-1.679767450000000e-09 * T3+3.688058805000000e-13 * T4
;
// species 20: OH
species[20] =
+3.381538120000000e+03 * invT+4.815738570000000e+00 -4.125305610000000e+00 * logT+1.612724695000000e-03 * T-1.087941151666667e-06 * T2+4.832113691666666e-10 * T3-1.031186895000000e-13 * T4
;
// species 22: NOH
species[22] =
+2.458501410000000e+04 * invT+3.778653500000000e-01 -4.307960970000000e+00 * logT+2.507692435000000e-03 * T-3.784242283333334e-06 * T2+2.213969758333333e-09 * T3-5.149934950000000e-13 * T4
;
// species 23: H2O
species[23] =
-3.029372670000000e+04 * invT+5.047672768000000e+00 -4.198640560000000e+00 * logT+1.018217050000000e-03 * T-1.086733685000000e-06 * T2+4.573308850000000e-10 * T3-8.859890850000000e-14 * T4
;
// species 24: NH3
species[24] =
-6.525488000000000e+03 * invT-5.922786000000000e+00 -2.204352000000000e+00 * logT-5.057380000000000e-03 * T+2.442108333333333e-06 * T2-1.206029166666667e-09 * T3+2.664254500000000e-13 * T4
;
// species 25: O2
species[25] =
-1.063943560000000e+03 * invT+1.247806300000001e-01 -3.782456360000000e+00 * logT+1.498367080000000e-03 * T-1.641217001666667e-06 * T2+8.067745908333334e-10 * T3-1.621864185000000e-13 * T4
;
// species 26: N2
species[26] =
-1.020899900000000e+03 * invT-6.516950000000001e-01 -3.298677000000000e+00 * logT-7.041202000000000e-04 * T+6.605369999999999e-07 * T2-4.701262500000001e-10 * T3+1.222427000000000e-13 * T4
;
}
else {
// species 1: O+
species[1] =
+1.879408740000000e+05 * invT-1.988834180000000e+00 -2.485420280000000e+00 * logT-1.284893475000000e-05 * T+2.147222966666667e-09 * T2-1.379379058333333e-13 * T3-5.496667200000000e-18 * T4
;
// species 2: O2+
species[2] =
+1.398768230000000e+05 * invT-2.130505470000000e+00 -3.316759220000000e+00 * logT-5.576112200000001e-04 * T+6.391542600000000e-08 * T2-4.773205725000000e-12 * T3+1.388241905000000e-16 * T4
;
// species 3: O4+
species[3] =
+1.708397930000000e+05 * invT+2.225873356000000e+01 -7.577843460000000e+00 * logT-1.240260260000000e-03 * T+1.618555716666667e-07 * T2-1.381152091666667e-11 * T3+5.152387000000000e-16 * T4
;
// species 4: O-
species[4] =
+1.879408740000000e+05 * invT-1.988834180000000e+00 -2.485420280000000e+00 * logT-1.284893475000000e-05 * T+2.147222966666667e-09 * T2-1.379379058333333e-13 * T3-5.496667200000000e-18 * T4
;
// species 5: O2-
species[5] =
+1.398768230000000e+05 * invT-2.130505470000000e+00 -3.316759220000000e+00 * logT-5.576112200000001e-04 * T+6.391542600000000e-08 * T2-4.773205725000000e-12 * T3+1.388241905000000e-16 * T4
;
// species 6: N4+
species[6] =
+2.145405450000000e+05 * invT+1.829579346000000e+01 -7.052858160000000e+00 * logT-1.455472345000000e-03 * T+1.859566650000000e-07 * T2-1.565699533333333e-11 * T3+5.788450650000000e-16 * T4
;
// species 7: N2+
species[7] =
+1.803909940000000e+05 * invT+4.907721999999999e-01 -3.586613630000000e+00 * logT-1.265359745000000e-04 * T-3.079636900000000e-08 * T2+3.793810191666666e-12 * T3-1.634090145000000e-16 * T4
;
// species 8: NH3+
species[8] =
-6.493270000000000e+03 * invT-5.010192999999999e+00 -2.461904000000000e+00 * logT-3.029583000000000e-03 * T+3.341628333333333e-07 * T2-2.613335833333333e-11 * T3+9.691585000000000e-16 * T4
;
// species 9: H2O+
species[9] =
+1.171220000000000e+05 * invT-1.932777610000000e+00 -3.033992490000000e+00 * logT-1.088459020000000e-03 * T+2.734541966666666e-08 * T2+8.086832250000000e-12 * T3-8.410049600000000e-16 * T4
;
// species 10: O(1D)
species[10] =
+5.199679930000000e+04 * invT-2.116556980000000e+00 -2.500000000000000e+00 * logT-1.046348745000000e-15 * T+1.871095050000000e-19 * T2-2.106359225000000e-23 * T3+1.100479930000000e-27 * T4
;
// species 11: O
species[11] =
+2.921757910000000e+04 * invT-2.214917859999999e+00 -2.569420780000000e+00 * logT+4.298705685000000e-05 * T-6.991409816666667e-09 * T2+8.348149916666666e-13 * T3-6.141684549999999e-17 * T4
;
// species 12: N
species[12] =
+5.611604000000000e+04 * invT-1.998490000000000e+00 -2.450268000000000e+00 * logT-5.330730000000000e-05 * T+1.244222833333333e-08 * T2-1.566376666666667e-12 * T3+5.129920000000000e-17 * T4
;
// species 13: N2(A1)
species[13] =
+7.064415680000000e+04 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 14: N2(B3)
species[14] =
+8.433141200000000e+04 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 15: N2(C3)
species[15] =
+1.270170890000000e+05 * invT-2.919311250000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 16: NH2
species[16] =
+2.191977000000000e+04 * invT-2.816567000000000e+00 -2.961311000000000e+00 * logT-1.466349500000000e-03 * T+1.510600000000000e-07 * T2-1.347714166666667e-11 * T3+6.021000000000000e-16 * T4
;
// species 17: NH
species[17] =
+4.207828000000000e+04 * invT-3.096950000000000e+00 -2.760249000000000e+00 * logT-6.876730000000000e-04 * T+7.419856666666667e-08 * T2-6.410660000000000e-12 * T3+2.508796000000000e-16 * T4
;
// species 18: H
species[18] =
+2.547365990000000e+04 * invT+2.946682924000000e+00 -2.500000010000000e+00 * logT+1.154214865000000e-11 * T-2.692699133333334e-15 * T2+3.945960291666667e-19 * T3-2.490986785000000e-23 * T4
;
// species 19: H2
species[19] =
-9.501589220000000e+02 * invT+6.542302510000000e+00 -3.337279200000000e+00 * logT+2.470123655000000e-05 * T-8.324279633333333e-08 * T2+1.496386616666667e-11 * T3-1.001276880000000e-15 * T4
;
// species 20: OH
species[20] =
+3.718857740000000e+03 * invT-2.836911870000000e+00 -2.864728860000000e+00 * logT-5.282522400000000e-04 * T+4.318045966666667e-08 * T2-2.543488950000000e-12 * T3+6.659793800000000e-17 * T4
;
// species 22: NOH
species[22] =
+2.439655860000000e+04 * invT-3.455839999999988e-03 -4.007457490000000e+00 * logT-1.209123330000000e-03 * T+1.356223335000000e-07 * T2-1.046442200000000e-11 * T3+3.640313010000000e-16 * T4
;
// species 23: H2O
species[23] =
-3.000429710000000e+04 * invT-1.932777610000000e+00 -3.033992490000000e+00 * logT-1.088459020000000e-03 * T+2.734541966666666e-08 * T2+8.086832250000000e-12 * T3-8.410049600000000e-16 * T4
;
// species 24: NH3
species[24] =
-6.493270000000000e+03 * invT-5.010192999999999e+00 -2.461904000000000e+00 * logT-3.029583000000000e-03 * T+3.341628333333333e-07 * T2-2.613335833333333e-11 * T3+9.691585000000000e-16 * T4
;
// species 25: O2
species[25] =
-1.088457720000000e+03 * invT-2.170693450000000e+00 -3.282537840000000e+00 * logT-7.415437700000000e-04 * T+1.263277781666667e-07 * T2-1.745587958333333e-11 * T3+1.083588970000000e-15 * T4
;
// species 26: N2
species[26] =
-9.227977000000000e+02 * invT-3.053888000000000e+00 -2.926640000000000e+00 * logT-7.439884000000000e-04 * T+9.474600000000001e-08 * T2-8.414198333333333e-12 * T3+3.376675500000000e-16 * T4
;
}

// species with midpoint at T=1500 kelvin
if (T < 1500) {
// species 21: H2NO
species[21] =
+6.868030000000000e+03 * invT-8.735920000000000e+00 -2.530590000000000e+00 * logT-4.298017500000000e-03 * T+9.118383333333333e-07 * T2-1.896874166666667e-10 * T3+2.324036500000000e-14 * T4
;
}
else {
// species 21: H2NO
species[21] =
+5.569325000000000e+03 * invT+1.182688600000000e+01 -5.673346000000000e+00 * logT-1.149418500000000e-03 * T+2.957410000000000e-08 * T2+9.195683333333332e-12 * T3-9.298810000000001e-16 * T4
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
return(mob);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real specDiff(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
amrex::Real diff = 0.0;
return(diff);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, const amrex::Real * sc, const amrex::Real T, const amrex::Real Te, amrex::Real EN, amrex::Real * enerExch)
{
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
const amrex::Real refC = 101325 / 8.31446 * invT;
const amrex::Real refCinv = 1 / refC;

for (int i = 0; i < 27; ++i) {
wdot[i] = 0.0;
}

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int i = 0; i < 27; ++i) {
mixture += sc[i];
}

// compute the Gibbs free energy
amrex::Real g_RT[27];
gibbs(g_RT, T);

// Precalculating values for electron energy exchange evaluation
amrex::Real ne = sc[E_ID] * 6.02214085774e23;
amrex::Real Ue = 1.5 * Te * ne * 1.380649e-23;

{
// reaction 0:  O2 + E => E + O2+ + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[2] += qdot;
wdot[25] -= qdot;
}

{
// reaction 1:  O2 + E => 2 O + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[11] += 2.000000 * qdot;
wdot[25] -= qdot;
}

{
// reaction 4:  O- + E => E + O + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[4]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[4] -= qdot;
wdot[11] += qdot;
}

{
// reaction 5:  O2 + O2+ + O2 => O4+ + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*(sc[25] * sc[25]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[3] += qdot;
wdot[25] -= qdot;
}

{
// reaction 7:  E + O2 + O2 => O2- + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*(sc[25] * sc[25]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] += qdot;
wdot[25] -= qdot;
}

{
// reaction 10:  O2- + O4+ + O2 => 3 O2 + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[3]*sc[5]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[5] -= qdot;
wdot[25] += 3.000000 * qdot;
}

{
// reaction 11:  O2+ + O2- + O2 => 2 O2 + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*sc[5]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[5] -= qdot;
wdot[25] += 2.000000 * qdot;
}

{
// reaction 12:  O2 + E => O + O(1D) + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] += qdot;
wdot[11] += qdot;
wdot[25] -= qdot;
}

{
// reaction 14:  O2- + O => E + O2 + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[5]*sc[11]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[5] -= qdot;
wdot[25] += qdot;
}

{
// reaction 16:  N2 + E => E + N2+ + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[7] += qdot;
wdot[26] -= qdot;
}

{
// reaction 17:  N2 + N2+ + N2 => N4+ + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[7]*(sc[26] * sc[26]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[6] += qdot;
wdot[7] -= qdot;
wdot[26] -= qdot;
}

{
// reaction 26:  N2 + N2+ + O2 => N4+ + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[7]*sc[25]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[6] += qdot;
wdot[7] -= qdot;
wdot[26] -= qdot;
}

{
// reaction 29:  O2 + O2+ + N2 => O4+ + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*sc[25]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[3] += qdot;
wdot[25] -= qdot;
}

{
// reaction 30:  O2- + O4+ + N2 => 3 O2 + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[3]*sc[5]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[5] -= qdot;
wdot[25] += 3.000000 * qdot;
}

{
// reaction 31:  O2+ + O2- + N2 => 2 O2 + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*sc[5]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[5] -= qdot;
wdot[25] += 2.000000 * qdot;
}

{
// reaction 37:  NH3 + E => H + NH2 + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[16] += qdot;
wdot[18] += qdot;
wdot[24] -= qdot;
}

{
// reaction 39:  NH3 + E => E + NH3+ + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[8] += qdot;
wdot[24] -= qdot;
}

{
// reaction 64:  H2O + E => H + OH + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[18] += qdot;
wdot[20] += qdot;
wdot[23] -= qdot;
}

{
// reaction 65:  H2O + E => H2 + O + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[11] += qdot;
wdot[19] += qdot;
wdot[23] -= qdot;
}

{
// reaction 66:  H2O + E => E + H2O+ + E
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[9] += qdot;
wdot[23] -= qdot;
}

{
// reaction 2:  E + O2 => O + O-
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[4] += qdot;
wdot[11] += qdot;
wdot[25] -= qdot;
}

{
// reaction 3:  O- + O2+ => O + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*sc[4]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[4] -= qdot;
wdot[11] += qdot;
wdot[25] += qdot;
}

{
// reaction 6:  E + O4+ => 2 O + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[3]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] -= qdot;
wdot[11] += 2.000000 * qdot;
wdot[25] += qdot;
}

{
// reaction 8:  E + O2+ => 2 O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[2]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[2] -= qdot;
wdot[11] += 2.000000 * qdot;
}

{
// reaction 9:  O2- + O4+ => 3 O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[3]*sc[5]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[5] -= qdot;
wdot[25] += 3.000000 * qdot;
}

{
// reaction 13:  O + O- => E + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[4]*sc[11]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[4] -= qdot;
wdot[11] -= qdot;
wdot[25] += qdot;
}

{
// reaction 15:  O(1D) + O2 => O + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[11] += qdot;
wdot[25] -= qdot;
wdot[25] += qdot;
}

{
// reaction 18:  E + N4+ => N2 + N2(C3)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[6]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[6] -= qdot;
wdot[15] += qdot;
wdot[26] += qdot;
}

{
// reaction 19:  E + N2+ => 2 N
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[7]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[7] -= qdot;
wdot[12] += 2.000000 * qdot;
}

{
// reaction 20:  E + N2 => E + N2(C3)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[15] += qdot;
wdot[26] -= qdot;
}

{
// reaction 21:  N2(C3) => N2(B3)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[15]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[14] += qdot;
wdot[15] -= qdot;
}

{
// reaction 22:  N2 + N2(C3) => N2 + N2(B3)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[15]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[14] += qdot;
wdot[15] -= qdot;
wdot[26] -= qdot;
wdot[26] += qdot;
}

{
// reaction 23:  E + N2 => E + N2(A1)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[13] += qdot;
wdot[26] -= qdot;
}

{
// reaction 24:  E + N2 => E + N2(B3)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[14] += qdot;
wdot[26] -= qdot;
}

{
// reaction 25:  N2 + N2(B3) => N2 + N2(A1)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[14]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[13] += qdot;
wdot[14] -= qdot;
wdot[26] -= qdot;
wdot[26] += qdot;
}

{
// reaction 27:  N4+ + O2 => 2 N2 + O2+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[6]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] += qdot;
wdot[6] -= qdot;
wdot[25] -= qdot;
wdot[26] += 2.000000 * qdot;
}

{
// reaction 28:  N2+ + O2 => N2 + O2+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[7]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] += qdot;
wdot[7] -= qdot;
wdot[25] -= qdot;
wdot[26] += qdot;
}

{
// reaction 32:  N2+ + O- => 2 N + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[4]*sc[7]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[4] -= qdot;
wdot[7] -= qdot;
wdot[11] += qdot;
wdot[12] += 2.000000 * qdot;
}

{
// reaction 33:  N2(A1) + O2 => N2 + 2 O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[13]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[11] += 2.000000 * qdot;
wdot[13] -= qdot;
wdot[25] -= qdot;
wdot[26] += qdot;
}

{
// reaction 34:  N2(B3) + O2 => N2 + 2 O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[14]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[11] += 2.000000 * qdot;
wdot[14] -= qdot;
wdot[25] -= qdot;
wdot[26] += qdot;
}

{
// reaction 35:  N2(C3) + O2 => N2 + O + O(1D)
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[15]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] += qdot;
wdot[11] += qdot;
wdot[15] -= qdot;
wdot[25] -= qdot;
wdot[26] += qdot;
}

{
// reaction 36:  N2 + O(1D) => N2 + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[26]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[11] += qdot;
wdot[26] -= qdot;
wdot[26] += qdot;
}

{
// reaction 38:  E + NH3 => E + 2 H + NH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[17] += qdot;
wdot[18] += 2.000000 * qdot;
wdot[24] -= qdot;
}

{
// reaction 40:  H + NH3 => H2 + NH2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[18]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[16] += qdot;
wdot[18] -= qdot;
wdot[19] += qdot;
wdot[24] -= qdot;
}

{
// reaction 41:  E + NH3+ => 2 H + NH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[8]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[8] -= qdot;
wdot[17] += qdot;
wdot[18] += 2.000000 * qdot;
}

{
// reaction 42:  E + NH3+ => H + NH2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[8]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[8] -= qdot;
wdot[16] += qdot;
wdot[18] += qdot;
}

{
// reaction 43:  NH3 + O(1D) => NH2 + OH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[16] += qdot;
wdot[20] += qdot;
wdot[24] -= qdot;
}

{
// reaction 44:  NH3 + O(1D) => H2O + NH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[17] += qdot;
wdot[23] += qdot;
wdot[24] -= qdot;
}

{
// reaction 45:  NH3 + O(1D) => H2 + NOH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[19] += qdot;
wdot[22] += qdot;
wdot[24] -= qdot;
}

{
// reaction 46:  NH3 + O(1D) => H + H2NO
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[18] += qdot;
wdot[21] += qdot;
wdot[24] -= qdot;
}

{
// reaction 47:  NH3+ + O2- => NH3 + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[5]*sc[8]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[5] -= qdot;
wdot[8] -= qdot;
wdot[24] += qdot;
wdot[25] += qdot;
}

{
// reaction 48:  NH3+ + O- => NH3 + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[4]*sc[8]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[4] -= qdot;
wdot[8] -= qdot;
wdot[11] += qdot;
wdot[24] += qdot;
}

{
// reaction 49:  NH3 + O2+ => NH3+ + O2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[2]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[8] += qdot;
wdot[24] -= qdot;
wdot[25] += qdot;
}

{
// reaction 50:  N2+ + NH3 => N2 + NH3+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[7]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[7] -= qdot;
wdot[8] += qdot;
wdot[24] -= qdot;
wdot[26] += qdot;
}

{
// reaction 51:  N + NH => H + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[12]*sc[17]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] -= qdot;
wdot[17] -= qdot;
wdot[18] += qdot;
wdot[26] += qdot;
}

{
// reaction 52:  H + NH => H2 + N
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[17]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] += qdot;
wdot[17] -= qdot;
wdot[18] -= qdot;
wdot[19] += qdot;
}

{
// reaction 53:  2 NH => H2 + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * ((sc[17] * sc[17]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[17] -= 2.000000 * qdot;
wdot[19] += qdot;
wdot[26] += qdot;
}

{
// reaction 54:  2 NH => N + NH2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * ((sc[17] * sc[17]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] += qdot;
wdot[16] += qdot;
wdot[17] -= 2.000000 * qdot;
}

{
// reaction 55:  2 NH => 2 H + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * ((sc[17] * sc[17]));
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[17] -= 2.000000 * qdot;
wdot[18] += 2.000000 * qdot;
wdot[26] += qdot;
}

{
// reaction 56:  H + NH2 => H2 + NH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[16]*sc[18]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[16] -= qdot;
wdot[17] += qdot;
wdot[18] -= qdot;
wdot[19] += qdot;
}

{
// reaction 57:  N + NH2 => 2 H + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[12]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] -= qdot;
wdot[16] -= qdot;
wdot[18] += 2.000000 * qdot;
wdot[26] += qdot;
}

{
// reaction 58:  N + NH2 => H2 + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[12]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] -= qdot;
wdot[16] -= qdot;
wdot[19] += qdot;
wdot[26] += qdot;
}

{
// reaction 59:  NH + NH2 => N + NH3
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[16]*sc[17]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[12] += qdot;
wdot[16] -= qdot;
wdot[17] -= qdot;
wdot[24] += qdot;
}

{
// reaction 60:  H2 + NH2 => H + NH3
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[16]*sc[19]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[16] -= qdot;
wdot[18] += qdot;
wdot[19] -= qdot;
wdot[24] += qdot;
}

{
// reaction 61:  N2(B3) + NH3 => H + N2 + NH2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[14]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[14] -= qdot;
wdot[16] += qdot;
wdot[18] += qdot;
wdot[24] -= qdot;
wdot[26] += qdot;
}

{
// reaction 62:  N2(A1) + NH3 => H + N2 + NH2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[13]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[13] -= qdot;
wdot[16] += qdot;
wdot[18] += qdot;
wdot[24] -= qdot;
wdot[26] += qdot;
}

{
// reaction 63:  N2(C3) + NH3 => H2 + N2 + NH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[15]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[15] -= qdot;
wdot[17] += qdot;
wdot[19] += qdot;
wdot[24] -= qdot;
wdot[26] += qdot;
}

{
// reaction 67:  E + H2O => H2 + O-
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[4] += qdot;
wdot[19] += qdot;
wdot[23] -= qdot;
}

{
// reaction 68:  E + H2O+ => H2 + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[9]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[9] -= qdot;
wdot[11] += qdot;
wdot[19] += qdot;
}

{
// reaction 69:  E + H2O+ => 2 H + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[9]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[9] -= qdot;
wdot[11] += qdot;
wdot[18] += 2.000000 * qdot;
}

{
// reaction 70:  E + H2O+ => H + OH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[0]*sc[9]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[9] -= qdot;
wdot[18] += qdot;
wdot[20] += qdot;
}

{
// reaction 71:  H2O + N2+ => H2O+ + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[7]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[7] -= qdot;
wdot[9] += qdot;
wdot[23] -= qdot;
wdot[26] += qdot;
}

{
// reaction 72:  H2O+ + NH3 => H2O + NH3+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[9]*sc[24]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[8] += qdot;
wdot[9] -= qdot;
wdot[23] += qdot;
wdot[24] -= qdot;
}

{
// reaction 73:  H2O+ + NH2 => NH3+ + OH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[9]*sc[16]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[8] += qdot;
wdot[9] -= qdot;
wdot[16] -= qdot;
wdot[20] += qdot;
}

{
// reaction 74:  H2O+ + O => H2 + O2+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[9]*sc[11]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] += qdot;
wdot[9] -= qdot;
wdot[11] -= qdot;
wdot[19] += qdot;
}

{
// reaction 75:  H2O+ + O2 => H2O + O2+
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[9]*sc[25]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[2] += qdot;
wdot[9] -= qdot;
wdot[23] += qdot;
wdot[25] -= qdot;
}

{
// reaction 76:  H2O + O(1D) => 2 OH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[20] += 2.000000 * qdot;
wdot[23] -= qdot;
}

{
// reaction 77:  H2O + O(1D) => H2O + O
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[10]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[10] -= qdot;
wdot[11] += qdot;
wdot[23] -= qdot;
wdot[23] += qdot;
}

{
// reaction 78:  H2O + N2(A1) => H + N2 + OH
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[13]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[13] -= qdot;
wdot[18] += qdot;
wdot[20] += qdot;
wdot[23] -= qdot;
wdot[26] += qdot;
}

{
// reaction 79:  H2O + N2(B3) => H2O + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[14]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[14] -= qdot;
wdot[23] -= qdot;
wdot[23] += qdot;
wdot[26] += qdot;
}

{
// reaction 80:  H2O + N2(C3) => H2O + N2
const amrex::Real k_f = 0
           ;
const amrex::Real qf = k_f * (sc[15]*sc[23]);
const amrex::Real qr = 0.0;
const amrex::Real qdot = qf - qr;
wdot[15] -= qdot;
wdot[23] -= qdot;
wdot[23] += qdot;
wdot[26] += qdot;
}


}


// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{

// convert to SI
for (int id = 0; id < 27; ++id) {
C[id] *= 1.0e6;
}

// convert to chemkin units
productionRate(wdot, C, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 27; ++id) {
C[id] *= 1.0e-6;
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[27]; // temporary storage
amrex::Real YOW = 0; 
amrex::Real PWORT; 

// Compute inverse of mean molecular wt first
for (int i = 0; i < 27; i++)
{
YOW += y[i]*imw(i);
}
// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// multiply by 1e6 so c goes to SI
PWORT *= 1e6; 
// Now compute conversion (and go to SI)
for (int i = 0; i < 27; i++)
{
c[i] = PWORT * y[i]*imw(i);
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 27; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[27]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 27; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 27; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[27]; // temporary storage

// See Eq 8 with an extra 1e6 so c goes to SI
for (int i = 0; i < 27; i++)
{
c[i] = 1e6 * rho * y[i]*imw(i);
}

// call productionRate
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 27; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho, const amrex::Real T, const amrex::Real x[],  amrex::Real wdot[],const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch)
{
amrex::Real c[27]; // temporary storage
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*0.000549; // E
XW += x[1]*15.998451; // O+
XW += x[2]*31.997451; // O2+
XW += x[3]*63.995451; // O4+
XW += x[4]*15.999549; // O-
XW += x[5]*31.998549; // O2-
XW += x[6]*56.027451; // N4+
XW += x[7]*28.013451; // N2+
XW += x[8]*17.030451; // NH3+
XW += x[9]*18.014451; // H2O+
XW += x[10]*15.999000; // O(1D)
XW += x[11]*15.999000; // O
XW += x[12]*14.007000; // N
XW += x[13]*28.014000; // N2(A1)
XW += x[14]*28.014000; // N2(B3)
XW += x[15]*28.014000; // N2(C3)
XW += x[16]*16.023000; // NH2
XW += x[17]*15.015000; // NH
XW += x[18]*1.008000; // H
XW += x[19]*2.016000; // H2
XW += x[20]*17.007000; // OH
XW += x[21]*32.022000; // H2NO
XW += x[22]*31.014000; // NOH
XW += x[23]*18.015000; // H2O
XW += x[24]*17.031000; // NH3
XW += x[25]*31.998000; // O2
XW += x[26]*28.014000; // N2
// Extra 1e6 factor to take c to SI
ROW = 1e6*rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 27; ++id) {
c[id] = x[id]*ROW;
}

// convert to chemkin units
productionRate(wdot, c, T, Te, EN, ener_exch);

// convert to chemkin units
for (int id = 0; id < 27; ++id) {
wdot[id] *= 1.0e-6;
}
}

//  species unit charge number 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
kcharge[0] = -1;// E
kcharge[1] = 1;// O+
kcharge[2] = 1;// O2+
kcharge[3] = 1;// O4+
kcharge[4] = -1;// O-
kcharge[5] = -1;// O2-
kcharge[6] = 1;// N4+
kcharge[7] = 1;// N2+
kcharge[8] = 1;// NH3+
kcharge[9] = 1;// H2O+
kcharge[10] = 0;// O(1D)
kcharge[11] = 0;// O
kcharge[12] = 0;// N
kcharge[13] = 0;// N2(A1)
kcharge[14] = 0;// N2(B3)
kcharge[15] = 0;// N2(C3)
kcharge[16] = 0;// NH2
kcharge[17] = 0;// NH
kcharge[18] = 0;// H
kcharge[19] = 0;// H2
kcharge[20] = 0;// OH
kcharge[21] = 0;// H2NO
kcharge[22] = 0;// NOH
kcharge[23] = 0;// H2O
kcharge[24] = 0;// NH3
kcharge[25] = 0;// O2
kcharge[26] = 0;// N2
}

//  species charge per unit mass 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

int kchrg[27];
CKCHRG(kchrg);

for (int id = 0; id < 27; ++id) {
zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
}
}

#endif

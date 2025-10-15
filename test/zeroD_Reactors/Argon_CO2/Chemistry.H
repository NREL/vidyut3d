#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <UnivConstants.H>
#include <VarDefines.H>

/* Elements
0  E
1  O
2  C
3  Ar
*/

// Species
#define E_ID 0
#define O_ID 1
#define O2_ID 2
#define AR_ID 3
#define CO2_ID 4
#define CO_ID 5
#define CO2v1_ID 6
#define CO2v2_ID 7
#define CO2v3_ID 8
#define CO2v4_ID 9
#define CO2p_ID 10
#define Om_ID 11
#define Op_ID 12
#define COp_ID 13
#define ARe_ID 14
#define ARp_ID 15

#define NUM_ELEMENTS 4
#define NUM_SPECIES 16

// FIX these and make species of same type continuous
#define FIRST_ION 10
#define NUM_IONS 5
#define FIRST_NEUTRAL 1
#define NUM_NEUTRALS 10
//==================================================

#define NUM_REACTIONS 68

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real* awt);
//  MISC
void CKAWT(amrex::Real* awt);
void CKNCF(int* ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int* _rmap);
void CKINU(const int i, int& nspec, int* ki, int* nu);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
    mm = 4;
    kk = 16;
    ii = 68;
    nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[16] = {
    1822.8884868472639482, // E
    0.0625039064941559,    // O
    0.0312519532470779,    // O2
    0.0250312891113892,    // AR
    0.0227226249176305,    // CO2
    0.0357015351660121,    // CO
    0.0227226249176305,    // CO2v1
    0.0227226249176305,    // CO2v2
    0.0227226249176305,    // CO2v3
    0.0227226249176305,    // CO2v4
    0.0227229081626688,    // CO2p
    0.0625017634094834,    // Om
    0.0625060497257993,    // Op
    0.0357022343994464,    // COp
    0.0250312891113892,    // ARe
    0.0250316328373182,    // ARp
};
#endif
const amrex::Real h_global_imw[16] = {
    1822.8884868472639482, // E
    0.0625039064941559,    // O
    0.0312519532470779,    // O2
    0.0250312891113892,    // AR
    0.0227226249176305,    // CO2
    0.0357015351660121,    // CO
    0.0227226249176305,    // CO2v1
    0.0227226249176305,    // CO2v2
    0.0227226249176305,    // CO2v3
    0.0227226249176305,    // CO2v4
    0.0227229081626688,    // CO2p
    0.0625017634094834,    // Om
    0.0625060497257993,    // Op
    0.0357022343994464,    // COp
    0.0250312891113892,    // ARe
    0.0250316328373182,    // ARp
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[16] = {
    0.000549,  // E
    15.999000, // O
    31.998000, // O2
    39.950000, // AR
    44.009000, // CO2
    28.010000, // CO
    44.009000, // CO2v1
    44.009000, // CO2v2
    44.009000, // CO2v3
    44.009000, // CO2v4
    44.008451, // CO2p
    15.999549, // Om
    15.998451, // Op
    28.009451, // COp
    39.950000, // ARe
    39.949451, // ARp
};
#endif
const amrex::Real h_global_mw[16] = {
    0.000549,  // E
    15.999000, // O
    31.998000, // O2
    39.950000, // AR
    44.009000, // CO2
    28.010000, // CO
    44.009000, // CO2v1
    44.009000, // CO2v2
    44.009000, // CO2v3
    44.009000, // CO2v4
    44.008451, // CO2p
    15.999549, // Om
    15.998451, // Op
    28.009451, // COp
    39.950000, // ARe
    39.949451, // ARp
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_imw(amrex::Real* imw_new)
{
    imw_new[0] = 1822.8884868472639482; // E
    imw_new[1] = 0.0625039064941559;    // O
    imw_new[2] = 0.0312519532470779;    // O2
    imw_new[3] = 0.0250312891113892;    // AR
    imw_new[4] = 0.0227226249176305;    // CO2
    imw_new[5] = 0.0357015351660121;    // CO
    imw_new[6] = 0.0227226249176305;    // CO2v1
    imw_new[7] = 0.0227226249176305;    // CO2v2
    imw_new[8] = 0.0227226249176305;    // CO2v3
    imw_new[9] = 0.0227226249176305;    // CO2v4
    imw_new[10] = 0.0227229081626688;   // CO2p
    imw_new[11] = 0.0625017634094834;   // Om
    imw_new[12] = 0.0625060497257993;   // Op
    imw_new[13] = 0.0357022343994464;   // COp
    imw_new[14] = 0.0250312891113892;   // ARe
    imw_new[15] = 0.0250316328373182;   // ARp
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real imw(const int n)
{
#if AMREX_DEVICE_COMPILE
    return global_imw[n];
#else
    return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_mw(amrex::Real* mw_new)
{
    mw_new[0] = 0.000549;   // E
    mw_new[1] = 15.999000;  // O
    mw_new[2] = 31.998000;  // O2
    mw_new[3] = 39.950000;  // AR
    mw_new[4] = 44.009000;  // CO2
    mw_new[5] = 28.010000;  // CO
    mw_new[6] = 44.009000;  // CO2v1
    mw_new[7] = 44.009000;  // CO2v2
    mw_new[8] = 44.009000;  // CO2v3
    mw_new[9] = 44.009000;  // CO2v4
    mw_new[10] = 44.008451; // CO2p
    mw_new[11] = 15.999549; // Om
    mw_new[12] = 15.998451; // Op
    mw_new[13] = 28.009451; // COp
    mw_new[14] = 39.950000; // ARe
    mw_new[15] = 39.949451; // ARp
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real mw(const int n)
{
#if AMREX_DEVICE_COMPILE
    return global_mw[n];
#else
    return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
    ru = 8.31446261815324e+07;
    ruc = 1.98721558317399615845;
    pa = 1.01325e+06;
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real* tc)
{

    // temperature
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1.0 / T;

    // species with midpoint at T=1000 kelvin
    if (T < 1000)
    {
        // species 0: E
        species[0] =
            -7.459783900000000e+02 * invT + 1.423710750000000e+01 -
            2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
            0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
            0.000000000000000e+00 * tc[4];
        // species 1: O
        species[1] =
            +2.912225920000000e+04 * invT + 1.116333640000000e+00 -
            3.168267100000000e+00 * tc[0] + 1.639659420000000e-03 * tc[1] -
            1.107177326666667e-06 * tc[2] + 5.106721866666666e-10 * tc[3] -
            1.056329855000000e-13 * tc[4];
        // species 2: O2
        species[2] =
            -1.063943560000000e+03 * invT + 1.247806300000001e-01 -
            3.782456360000000e+00 * tc[0] + 1.498367080000000e-03 * tc[1] -
            1.641217001666667e-06 * tc[2] + 8.067745908333334e-10 * tc[3] -
            1.621864185000000e-13 * tc[4];
        // species 4: CO2
        species[4] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 5: CO
        species[5] =
            -1.434408600000000e+04 * invT + 7.112418999999992e-02 -
            3.579533470000000e+00 * tc[0] + 3.051768400000000e-04 * tc[1] -
            1.694690550000000e-07 * tc[2] - 7.558382366666667e-11 * tc[3] +
            4.522122495000000e-14 * tc[4];
        // species 6: CO2v1
        species[6] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 7: CO2v2
        species[7] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 8: CO2v3
        species[8] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 9: CO2v4
        species[9] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 10: CO2p
        species[10] =
            -4.837196970000000e+04 * invT - 7.544278700000000e+00 -
            2.356773520000000e+00 * tc[0] - 4.492298385000000e-03 * tc[1] +
            1.187260448333333e-06 * tc[2] - 2.049325183333333e-10 * tc[3] +
            7.184977399999999e-15 * tc[4];
        // species 11: Om
        species[11] =
            +2.912225920000000e+04 * invT + 1.116333640000000e+00 -
            3.168267100000000e+00 * tc[0] + 1.639659420000000e-03 * tc[1] -
            1.107177326666667e-06 * tc[2] + 5.106721866666666e-10 * tc[3] -
            1.056329855000000e-13 * tc[4];
        // species 12: Op
        species[12] =
            +2.912225920000000e+04 * invT + 1.116333640000000e+00 -
            3.168267100000000e+00 * tc[0] + 1.639659420000000e-03 * tc[1] -
            1.107177326666667e-06 * tc[2] + 5.106721866666666e-10 * tc[3] -
            1.056329855000000e-13 * tc[4];
        // species 13: COp
        species[13] =
            -1.434408600000000e+04 * invT + 7.112418999999992e-02 -
            3.579533470000000e+00 * tc[0] + 3.051768400000000e-04 * tc[1] -
            1.694690550000000e-07 * tc[2] - 7.558382366666667e-11 * tc[3] +
            4.522122495000000e-14 * tc[4];
    } else
    {
        // species 0: E
        species[0] =
            -7.459784500000000e+02 * invT + 1.423710750000000e+01 -
            2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
            0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
            0.000000000000000e+00 * tc[4];
        // species 1: O
        species[1] =
            +2.921757910000000e+04 * invT - 2.214917859999999e+00 -
            2.569420780000000e+00 * tc[0] + 4.298705685000000e-05 * tc[1] -
            6.991409816666667e-09 * tc[2] + 8.348149916666666e-13 * tc[3] -
            6.141684549999999e-17 * tc[4];
        // species 2: O2
        species[2] =
            -1.088457720000000e+03 * invT - 2.170693450000000e+00 -
            3.282537840000000e+00 * tc[0] - 7.415437700000000e-04 * tc[1] +
            1.263277781666667e-07 * tc[2] - 1.745587958333333e-11 * tc[3] +
            1.083588970000000e-15 * tc[4];
        // species 4: CO2
        species[4] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 5: CO
        species[5] =
            -1.415187240000000e+04 * invT - 5.103502110000000e+00 -
            2.715185610000000e+00 * tc[0] - 1.031263715000000e-03 * tc[1] +
            1.664709618333334e-07 * tc[2] - 1.917108400000000e-11 * tc[3] +
            1.018238580000000e-15 * tc[4];
        // species 6: CO2v1
        species[6] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 7: CO2v2
        species[7] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 8: CO2v3
        species[8] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 9: CO2v4
        species[9] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 10: CO2p
        species[10] =
            -4.875916600000000e+04 * invT + 1.585822230000000e+00 -
            3.857460290000000e+00 * tc[0] - 2.207185130000000e-03 * tc[1] +
            3.691356733333334e-07 * tc[2] - 4.362418233333334e-11 * tc[3] +
            2.360420820000000e-15 * tc[4];
        // species 11: Om
        species[11] =
            +2.921757910000000e+04 * invT - 2.214917859999999e+00 -
            2.569420780000000e+00 * tc[0] + 4.298705685000000e-05 * tc[1] -
            6.991409816666667e-09 * tc[2] + 8.348149916666666e-13 * tc[3] -
            6.141684549999999e-17 * tc[4];
        // species 12: Op
        species[12] =
            +2.921757910000000e+04 * invT - 2.214917859999999e+00 -
            2.569420780000000e+00 * tc[0] + 4.298705685000000e-05 * tc[1] -
            6.991409816666667e-09 * tc[2] + 8.348149916666666e-13 * tc[3] -
            6.141684549999999e-17 * tc[4];
        // species 13: COp
        species[13] =
            -1.415187240000000e+04 * invT - 5.103502110000000e+00 -
            2.715185610000000e+00 * tc[0] - 1.031263715000000e-03 * tc[1] +
            1.664709618333334e-07 * tc[2] - 1.917108400000000e-11 * tc[3] +
            1.018238580000000e-15 * tc[4];
    }

    // species with no change at a midpoint T
    // species 3: AR
    species[3] = -7.453750000000000e+02 * invT - 1.866000000000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 14: ARe
    species[14] = -7.453750000000000e+02 * invT - 1.866000000000000e+00 -
                  2.500000000000000e+00 * tc[0] -
                  0.000000000000000e+00 * tc[1] -
                  0.000000000000000e+00 * tc[2] -
                  0.000000000000000e+00 * tc[3] - 0.000000000000000e+00 * tc[4];
    // species 15: ARp
    species[15] = -7.453750000000000e+02 * invT - 1.866000000000000e+00 -
                  2.500000000000000e+00 * tc[0] -
                  0.000000000000000e+00 * tc[1] -
                  0.000000000000000e+00 * tc[2] -
                  0.000000000000000e+00 * tc[3] - 0.000000000000000e+00 * tc[4];
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT(amrex::Real wt[])
{
    get_mw(wt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real specMob(
    const int specid,
    const amrex::Real Te,
    const amrex::Real ndens,
    const amrex::Real emag,
    const amrex::Real T)
{

    amrex::Real Arp_cs = 120e-20; // A^2 (cross section)
    amrex::Real E_cs = 40e-20;    // A^2 (cross section)
    amrex::Real nu_elec = std::sqrt(8.0 * K_B * Te / PI / ME) * ndens * E_cs;
    amrex::Real elecmob = -ECHARGE / ME / nu_elec;

    amrex::Real mob = 0.0;
    if (specid == E_ID)
    {
        mob = elecmob;
    }

    return (mob);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real specDiff(
    const int specid,
    const amrex::Real Te,
    const amrex::Real ndens,
    const amrex::Real emag,
    const amrex::Real T)
{

    return 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_ener_exch(
    amrex::Real qf,
    amrex::Real qr,
    const amrex::Real* sc,
    amrex::Real kf,
    int rxntype,
    amrex::Real eexci,
    int elidx,
    amrex::Real* enerExch,
    amrex::Real Ue,
    amrex::Real T,
    const amrex::Real Te)
{
    // (0)='no electron'
    // (1)='excitation/de-excitation'
    // (2)='ionization'
    // (3)='recombination'
    // (4)='attachment'
    // (5)='elastic exchange'
    // (6)='chemi-ioniz'
    // (7)='detachment'

    amrex::Real eV2J = 1.602176620800000e-19; // eV -> J
    amrex::Real NA = 6.02214085774e23;        // 1/mol
    amrex::Real kB = 1.380649e-23;            // J/K
    amrex::Real ne = sc[E_ID] * NA;           // 1/m3
    amrex::Real mee = Ue / ne;                // J

    // Molecular weight used to calculate elastic collision energy exchange
    amrex::Real wt[NUM_SPECIES];
    get_mw(wt);

    // Superelastic factor (keeps Te from dropping below gas temperarture)
    amrex::Real sefact = (tanh((Te - T) / 0.5) + 1.0) / 2.0;

    // *enerExch = 0.0;
    if (rxntype == 1 || rxntype == 2)
    {
        *enerExch -= (eexci * eV2J) * NA * (qf - qr) * sefact;
    }
    if (rxntype == 3)
    {
        *enerExch -= mee * NA * (qf - qr);
    }
    if (rxntype == 4)
    {
        *enerExch -= mee * NA * (qf - qr);
    }
    if (rxntype == 5)
    {
        amrex::Real qq = sc[E_ID] * sc[elidx] * kf;
        *enerExch -= 3.0 * kB * (wt[E_ID] / wt[elidx]) * (Te - T) * qq * NA;
    }
    if (rxntype == 7)
    {
        *enerExch += mee * NA * (qf - qr);
    }

    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void productionRate(
    amrex::Real* wdot,
    const amrex::Real* sc,
    const amrex::Real T,
    const amrex::Real Te,
    amrex::Real* enerExch)
{
    const amrex::Real tc[5] = {
        log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];

    // reference concentration: P_atm / (RT) in inverse mol/m^3
    const amrex::Real refC = 101325 / 8.31446 * invT;
    const amrex::Real refCinv = 1 / refC;

    for (int i = 0; i < 16; ++i)
    {
        wdot[i] = 0.0;
    }

    // compute the mixture concentration
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 16; ++i)
    {
        mixture += sc[i];
    }

    // compute the Gibbs free energy
    amrex::Real g_RT[16];
    gibbs(g_RT, tc);

    // Precalculating values for electron energy exchange evaluation
    amrex::Real ne = sc[E_ID] * 6.02214085774e23;
    amrex::Real Ue = 1.5 * Te * ne * 1.380649e-23;

    // NOTE: units of JANEV fits for electron impact rxns are cm3/s and must be
    // convered to m3/mol-s Precalculating values
    double Janev_sum;
    double invTe = (Te == 0) ? 1.0 : 1.0 / Te;
    double TeeV = Te / 11595.0;
    double logTe =
        log(TeeV); // Fits are performed assuming Te is eV rather than K
    double invTeeV = (Te == 0) ? 1.0 : 1.0 / (TeeV);
    double Te_pow[] = {
        1.0,
        logTe,
        amrex::Math::powi<2>(logTe),
        amrex::Math::powi<3>(logTe),
        amrex::Math::powi<4>(logTe),
        amrex::Math::powi<5>(logTe),
        amrex::Math::powi<6>(logTe),
        amrex::Math::powi<7>(logTe),
        amrex::Math::powi<8>(logTe)};
    double invTe_pow[] = {
        invTeeV, amrex::Math::powi<2>(invTeeV), amrex::Math::powi<3>(invTeeV),
        amrex::Math::powi<4>(invTeeV)};

    std::vector<double> Jfit_coefs = {0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0};
    std::vector<double> Ffit_coefs = {0.0, 0.0, 0.0, 0.0};

    {
        // reaction 5:  CO2 + E => CO + O + E
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            28.3950215483559, -119.510009432235, 160.437496467960,
            -74.1425574458809};
        double Ffit_A = 1.50741518933862e-16;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
    }

    {
        // reaction 6:  CO2 + E => CO2+ + E + E
        // OLD FIT USING PURE CO2 CROSS SECTION DATA
        // Janev_sum = 0.0;
        // amrex::Real k_f;
        // if(TeeV < 8.0){
        //     Jfit_coefs = {-73.2131885836362, 288.352725017751,
        //     -674.533620464049, 949.480977796531, -825.714286521495,
        //     451.196428429444, -151.538614955178, 28.6691180230888,
        //     -2.34367647830683}; double Jfit_A = 3.99500000000000e-16; for(int
        //     j = 0; j<9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j]; k_f =
        //     Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        // } else {
        //     Ffit_coefs = {-41.5052482777841, 19.3784093595350,
        //     -60.4073644704944, 74.5295616459192}; double Ffit_A
        //     = 4.15042653485891e-14; for(int j = 0; j<4; j++) Janev_sum +=
        //     Ffit_coefs[j] * invTe_pow[j]; k_f = Ffit_A * exp(Janev_sum)
        //     * 6.02214085774e23;
        // }
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            65.5733184943082, -511.020275942352, 911.053057794656,
            -540.822599138513};
        double Ffit_A = 1.12781185638050e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[10] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 13.8;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 10:  CO + E => CO+ + E + E
        Janev_sum = 0.0;
        amrex::Real k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -2290.33055644355, 11806.3672878286,  -26408.0009783117,
                33193.6096823637,  -25629.3277592247, 12465.3156663098,
                -3735.92751665149, 631.802846717600,  -46.2207594603064};
            double Jfit_A = 1.05300000000000e-14;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -61.3900740209655, 436.492020460250, -2057.32080961874,
                2715.83171383976};
            double Ffit_A = 3.64027197545915e-13;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const amrex::Real qf = k_f * (sc[0] * sc[5]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[5] -= qdot;
        wdot[13] += qdot;
    }

    {
        // reaction 13:  AR + E => AR+ + E + E
        Janev_sum = 0.0;
        amrex::Real k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1951.40116682516, 9125.28938191463,  -19064.2751720459,
                22763.7019090326,  -16855.4371729404, 7904.47881474249,
                -2291.72587558317, 375.740475221758,  -26.6926001068518};
            double Jfit_A = 1.42900000000000e-14;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -77.9510931269820, 533.060249822779, -2091.33857389100,
                1893.83145115372};
            double Ffit_A = 1.13390537981533e-12;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const amrex::Real qf = k_f * (sc[0] * sc[3]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[3] -= qdot;
        wdot[15] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 15.759;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 15:  E + O+ + M => O + M
        const amrex::Real k_f = 3626600000;
        const amrex::Real Corr = mixture;
        const amrex::Real qf = Corr * k_f * (sc[0] * sc[12]);
        const amrex::Real qr =
            Corr * k_f * exp(-(g_RT[0] - g_RT[1] + g_RT[12])) * (refC) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[1] += qdot;
        wdot[12] -= qdot;
    }

    {
        // reaction 17:  E + O + M => O- + M
        const amrex::Real k_f = 36266;
        const amrex::Real Corr = mixture;
        const amrex::Real qf = Corr * k_f * (sc[0] * sc[1]);
        const amrex::Real qr =
            Corr * k_f * exp(-(g_RT[0] + g_RT[1] - g_RT[11])) * (refC) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[1] -= qdot;
        wdot[11] += qdot;
    }

    {
        // reaction 63:  O+ + O- + M => O2 + M
        const amrex::Real k_f = 72530000000;
        const amrex::Real Corr = mixture;
        const amrex::Real qf = Corr * k_f * (sc[11] * sc[12]);
        const amrex::Real qr = Corr * k_f *
                               exp(-(-g_RT[2] + g_RT[11] + g_RT[12])) * (refC) *
                               (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[2] += qdot;
        wdot[11] -= qdot;
        wdot[12] -= qdot;
    }

    {
        // reaction 65:  CO + O + M => CO2 + M
        const amrex::Real k_f = 297.38267 * exp(-1510.0 / tc[1]);
        const amrex::Real Corr = mixture;
        const amrex::Real qf = Corr * k_f * (sc[1] * sc[5]);
        const amrex::Real qr =
            Corr * k_f * exp(-(g_RT[1] - g_RT[4] + g_RT[5])) * (refC) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[4] += qdot;
        wdot[5] -= qdot;
    }

    {
        // reaction 67:  2 O + M => O2 + M
        const amrex::Real k_f =
            4605.804775 * (300.0 / tc[1]) * exp(-170.0 / tc[1]);
        const amrex::Real Corr = mixture;
        const amrex::Real qf = Corr * k_f * ((sc[1] * sc[1]));
        const amrex::Real qr =
            Corr * k_f * exp(-(2.000000 * g_RT[1] - g_RT[2])) * (refC) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[1] -= 2.000000 * qdot;
        wdot[2] += qdot;
    }

    {
        // reaction 0:  CO2 + E => CO2 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.74567823745404,   0.489892009971784,    0.475497844841483,
            -0.0828111482837000, -0.0869447627283896,  0.00364885418266187,
            0.00738742037568563, 0.000815768871163232, -1.80446057296271e-05};
        double Jfit_A = 3.07100000000000e-13;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        // Capping rate for now since fit diverges at large Te
        const amrex::Real k_f =
            std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 1.0e12);
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        int rxntype = 5;
        amrex::Real eexci = 0.0;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 0:  CO2 + E => CO2 + E
        // NOTE: excitational losses for CO2(e1)
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            13.9944258667849, -44.2367004233860, 35.3755884340306,
            -12.2196503543953};
        double Ffit_A = 4.66040157646352e-16;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 7.0;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 0:  CO2 + E => CO2 + E
        // NOTE: excitational losses for CO2(e2)
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            -12.6852242013963, -39.8289588750467, 75.9272218485778,
            -44.8968368578706};
        double Ffit_A = 1.41992889938866e-13;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 10.5;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 1:  CO2 + E => CO2(v1) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.502839945828532, 0.694016395482428,    0.190052560381054,
            -0.336722426572758, -0.112995969891858,   0.0625050243076499,
            0.0180812366062628, -0.00406798152564209, -0.00117479961936418};
        double Jfit_A = 9.13200000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[6] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.083;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 2:  CO2 + E => CO2(v2) + E
        // NOTE: CO2(020)
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.08379968039393,  1.78557074369457,    0.0320873806354341,
            -0.785497270396418, -0.0745881449840463, 0.173672656884839,
            0.0242491592490644, -0.0137085592778002, -0.00283607126964260};
        double Jfit_A = 4.31900000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[7] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.17;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 2:  CO2 + E => CO2(v2) + E
        // NOTE: CO2(100)
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.74710336851205,   5.28866024293583,    -3.06722524904491,
            -1.26813972992265,   0.765286490697251,   0.220019094021629,
            -0.0660599223620025, -0.0174513395327540, -0.000100232840708316};
        double Jfit_A = 2.10800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[7] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.17;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.26399041657649,  4.01593315343776,  -3.86564987870520,
            0.630843859313518,  0.918455783525699, -0.466408614313122,
            -0.122041546730277, 0.133320282875472, -0.0249335094184505};
        double Jfit_A = 1.22800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.252;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 3:  CO2 + E => CO2(v3) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.0161587259013659, 0.192761580775350,    -0.233095798484862,
            -0.0504138326883387, 0.0729053506151256,   0.0268158476130099,
            -0.0130233553843901, -0.00234706991825410, 0.000208569444277222};
        double Jfit_A = 5.87800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[8] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.291;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.430350115379283,  1.88554119808260,   -2.10311501560593,
            0.205835857388545,   0.349603651488484,  -0.113940765190612,
            -0.0356282206715671, 0.0288224520101671, -0.00505168357993776};
        double Jfit_A = 6.65200000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.339;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }
    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.26866503084565,  4.16951503938324,  -4.27011214236418,
            0.528779192020846,  1.67719403845492,  -0.767308824575880,
            -0.423310590818580, 0.369354928743402, -0.0691523136873366};
        double Jfit_A = 1.02300000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.339;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }
    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.16752755372413,  3.97264512328113,  -4.02945519654257,
            0.673891717708656,  1.00483069462086,  -0.508781573502289,
            -0.133736912597309, 0.144832888470466, -0.0269944163072462};
        double Jfit_A = 3.90400000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.422;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }
    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.57218058544877,  4.53507634411035,  -3.90997647678875,
            0.335368207124721,  1.05569484438475,  -0.410009273039837,
            -0.159898555151304, 0.129957825662046, -0.0219563391933853};
        double Jfit_A = 2.89300000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }
    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.803462147557306,  3.31431359636670,   -3.80327983310172,
            0.791648348030451,   0.692211498451406,  -0.342472582689845,
            -0.0688634709416767, 0.0708789327184151, -0.0117253943723200};
        double Jfit_A = 1.68100000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 0.505;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }
    {
        // reaction 4:  CO2 + E => CO2(v4) + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.13341741015961,  3.84235032466761,  -3.89843865763089,
            0.695260507371384,  0.918161189362339, -0.486298345254785,
            -0.117448184645462, 0.135925105891503, -0.0258281763494971};
        double Jfit_A = 5.46700000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[4] -= qdot;
        wdot[9] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 2.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 7:  CO2 + E => CO + 2 E + O+
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            -77.4751667545915, 368.702583061800, -1027.52886179464,
            -875.218481131602};
        double Ffit_A = 1.67574146792253e-13;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += 2.000000 * qdot;
        wdot[0] -= qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[12] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 19.1;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 8:  CO2 + E => CO+ + 2 E + O
        Janev_sum = 0.0;
        amrex::Real k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -16374.5761168493, 70345.6825911028,  -132493.614520093,
                142424.193083429,  -95424.9718989068, 40779.8217137456,
                -10852.6008636845, 1644.38446103563,  -108.620096064776};
            double Jfit_A = 1.10800000000000e-15;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -117.796225708354, 928.699594978788, -3999.36334971610,
                4171.87939848509};
            double Ffit_A = 4.84703792532310e-13;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += 2.000000 * qdot;
        wdot[0] -= qdot;
        wdot[1] += qdot;
        wdot[4] -= qdot;
        wdot[13] += qdot;
        int rxntype = 2;
        amrex::Real eexci = 19.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 9:  CO2 + E => CO + O-
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            5.68537810168761, -15.1912011431635, 6.86995641002507,
            -1.21745244042248};
        double Ffit_A = 4.23133262631561e-18;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[11] += qdot;
        int rxntype = 4;
        amrex::Real eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 11:  AR + E => AR + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -2.62130096848181,    1.55888667906962,     0.0350711485300335,
            0.00439515704875221,  -0.00138205864597430, -0.0293026895863627,
            -0.00284458054782676, 0.00241786232297547,  0.000390470745241025};
        double Jfit_A = 2.33300000000000e-13;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const amrex::Real k_f = std::min(
            Jfit_A * exp(Janev_sum) * 6.02214085774e23,
            3.0e-13 * 6.02214085774e23);
        const amrex::Real qf = k_f * (sc[0] * sc[3]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        int rxntype = 5;
        amrex::Real eexci = 0.0;
        int elidx = AR_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 12:  AR + E => AR* + E
        Janev_sum = 0.0;
        amrex::Real k_f;
        Ffit_coefs = {
            -17.8244315740958, -55.1739383730274, 102.978958246859,
            -55.9966048149755};
        double Ffit_A = 9.18800465289008e-14;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const amrex::Real qf = k_f * (sc[0] * sc[3]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[0] += qdot;
        wdot[3] -= qdot;
        wdot[14] += qdot;
        int rxntype = 1;
        amrex::Real eexci = 11.55;
        int elidx = AR_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 14:  CO2+ + E => CO + O
        const amrex::Real k_f =
            12044281520000 * (1.0 / tc[1]) * pow(TeeV, -0.5);
        const amrex::Real qf = k_f * (sc[0] * sc[10]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[1] += qdot;
        wdot[5] += qdot;
        wdot[10] -= qdot;
    }

    {
        // reaction 16:  CO+ + E => CO
        const amrex::Real k_f = 2.217954441908e10 * pow(TeeV, -0.55);
        const amrex::Real qf = k_f * (sc[0] * sc[13]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[5] += qdot;
        wdot[13] -= qdot;
    }

    {
        // reaction 18:  AR* + CO2 => AR + CO + O
        const amrex::Real k_f = 319200000;
        const amrex::Real qf = k_f * (sc[4] * sc[14]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[3] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[14] -= qdot;
    }

    {
        // reaction 19:  AR+ + CO2 => AR + CO + O+
        const amrex::Real k_f = 0;
        const amrex::Real qf = k_f * (sc[4] * sc[15]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[12] += qdot;
        wdot[15] -= qdot;
    }

    {
        // reaction 20:  AR+ + CO2 => AR + CO+ + O
        const amrex::Real k_f = 0;
        const amrex::Real qf = k_f * (sc[4] * sc[15]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[3] += qdot;
        wdot[4] -= qdot;
        wdot[13] += qdot;
        wdot[15] -= qdot;
    }

    {
        // reaction 21:  AR+ + CO2 => AR + CO2+
        const amrex::Real k_f = 30110000;
        const amrex::Real qf = k_f * (sc[4] * sc[15]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] += qdot;
        wdot[4] -= qdot;
        wdot[10] += qdot;
        wdot[15] -= qdot;
    }

    {
        // reaction 22:  CO2 + CO2(v1) => 2 CO2
        const amrex::Real k_f = 6444;
        const amrex::Real qf = k_f * (sc[4] * sc[6]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += 2.000000 * qdot;
        wdot[6] -= qdot;
    }

    {
        // reaction 23:  CO + CO2(v1) => CO + CO2
        const amrex::Real k_f = 4505;
        const amrex::Real qf = k_f * (sc[5] * sc[6]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] += qdot;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[6] -= qdot;
    }

    {
        // reaction 24:  CO2(v1) + O2 => CO2 + O2
        const amrex::Real k_f = 4505;
        const amrex::Real qf = k_f * (sc[2] * sc[6]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[4] += qdot;
        wdot[6] -= qdot;
    }

    {
        // reaction 25:  AR + CO2(v1) => AR + CO2
        const amrex::Real k_f = 4505;
        const amrex::Real qf = k_f * (sc[3] * sc[6]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[4] += qdot;
        wdot[6] -= qdot;
    }

    {
        // reaction 26:  CO2 + CO2(v2) => 2 CO2
        const amrex::Real k_f = 5.42;
        const amrex::Real qf = k_f * (sc[4] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += 2.000000 * qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 27:  CO + CO2(v2) => CO + CO2
        const amrex::Real k_f = 16.8;
        const amrex::Real qf = k_f * (sc[5] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] += qdot;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 28:  CO2(v2) + O2 => CO2 + O2
        const amrex::Real k_f = 16.8;
        const amrex::Real qf = k_f * (sc[2] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[4] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 29:  AR + CO2(v2) => AR + CO2
        const amrex::Real k_f = 16.8;
        const amrex::Real qf = k_f * (sc[3] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[4] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 30:  CO2 + CO2(v2) => CO2 + CO2(v1)
        const amrex::Real k_f = 17460;
        const amrex::Real qf = k_f * (sc[4] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 31:  CO + CO2(v2) => CO + CO2(v1)
        const amrex::Real k_f = 12220;
        const amrex::Real qf = k_f * (sc[5] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 32:  CO2(v2) + O2 => CO2(v1) + O2
        const amrex::Real k_f = 12220;
        const amrex::Real qf = k_f * (sc[2] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 33:  AR + CO2(v2) => AR + CO2(v1)
        const amrex::Real k_f = 12220;
        const amrex::Real qf = k_f * (sc[3] * sc[7]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // reaction 34:  CO2 + CO2(v3) => CO2 + CO2(v2)
        const amrex::Real k_f = 464.9;
        const amrex::Real qf = k_f * (sc[4] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[7] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 35:  CO + CO2(v3) => CO + CO2(v2)
        const amrex::Real k_f = 13.97;
        const amrex::Real qf = k_f * (sc[5] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[7] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 36:  CO2(v3) + O2 => CO2(v2) + O2
        const amrex::Real k_f = 186.1;
        const amrex::Real qf = k_f * (sc[2] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[7] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 37:  AR + CO2(v3) => AR + CO2(v2)
        const amrex::Real k_f = 186.1;
        const amrex::Real qf = k_f * (sc[3] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[7] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 38:  CO2 + CO2(v3) => CO2 + CO2(v4)
        const amrex::Real k_f = 3643;
        const amrex::Real qf = k_f * (sc[4] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[8] -= qdot;
        wdot[9] += qdot;
    }

    {
        // reaction 39:  CO + CO2(v3) => CO + CO2(v4)
        const amrex::Real k_f = 1090;
        const amrex::Real qf = k_f * (sc[5] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[8] -= qdot;
        wdot[9] += qdot;
    }

    {
        // reaction 40:  CO2(v3) + O2 => CO2(v4) + O2
        const amrex::Real k_f = 1457;
        const amrex::Real qf = k_f * (sc[2] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[8] -= qdot;
        wdot[9] += qdot;
    }

    {
        // reaction 41:  AR + CO2(v3) => AR + CO2(v4)
        const amrex::Real k_f = 1457;
        const amrex::Real qf = k_f * (sc[3] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[8] -= qdot;
        wdot[9] += qdot;
    }

    {
        // reaction 42:  CO2 + CO2(v3) => CO2(v1) + CO2(v2)
        const amrex::Real k_f = 1457;
        const amrex::Real qf = k_f * (sc[4] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[6] += qdot;
        wdot[7] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 43:  CO2 + CO2(v3) => CO2 + CO2(v1)
        const amrex::Real k_f = 1.024;
        const amrex::Real qf = k_f * (sc[4] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[6] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 44:  CO + CO2(v3) => CO + CO2(v1)
        const amrex::Real k_f = 0.3071;
        const amrex::Real qf = k_f * (sc[5] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[6] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 45:  CO2(v3) + O2 => CO2(v1) + O2
        const amrex::Real k_f = 0.4095;
        const amrex::Real qf = k_f * (sc[2] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[6] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 46:  AR + CO2(v3) => AR + CO2(v1)
        const amrex::Real k_f = 0.4095;
        const amrex::Real qf = k_f * (sc[3] * sc[8]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[6] += qdot;
        wdot[8] -= qdot;
    }

    {
        // reaction 47:  CO2 + CO2(v4) => CO2 + CO2(v2)
        const amrex::Real k_f = 26080;
        const amrex::Real qf = k_f * (sc[4] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[7] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 48:  CO + CO2(v4) => CO + CO2(v2)
        const amrex::Real k_f = 18250;
        const amrex::Real qf = k_f * (sc[5] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[7] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 49:  CO2(v4) + O2 => CO2(v2) + O2
        const amrex::Real k_f = 18250;
        const amrex::Real qf = k_f * (sc[2] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[7] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 50:  AR + CO2(v4) => AR + CO2(v2)
        const amrex::Real k_f = 18250;
        const amrex::Real qf = k_f * (sc[3] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[7] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 51:  CO2 + CO2(v4) => CO2 + CO2(v1)
        const amrex::Real k_f = 5.468;
        const amrex::Real qf = k_f * (sc[4] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[4] += qdot;
        wdot[6] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 52:  CO + CO2(v4) => CO + CO2(v1)
        const amrex::Real k_f = 3722;
        const amrex::Real qf = k_f * (sc[5] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[5] -= qdot;
        wdot[5] += qdot;
        wdot[6] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 53:  CO2(v4) + O2 => CO2(v1) + O2
        const amrex::Real k_f = 3722;
        const amrex::Real qf = k_f * (sc[2] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[2] += qdot;
        wdot[6] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 54:  AR + CO2(v4) => AR + CO2(v1)
        const amrex::Real k_f = 3722;
        const amrex::Real qf = k_f * (sc[3] * sc[9]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[3] += qdot;
        wdot[6] += qdot;
        wdot[9] -= qdot;
    }

    {
        // reaction 55:  CO2 + O+ => CO2+ + O
        const amrex::Real k_f = 271000000;
        const amrex::Real qf = k_f * (sc[4] * sc[12]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[4] -= qdot;
        wdot[10] += qdot;
        wdot[12] -= qdot;
    }

    {
        // reaction 56:  CO+ + CO2 => CO + CO2+
        const amrex::Real k_f = 602200000;
        const amrex::Real qf = k_f * (sc[4] * sc[13]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[10] += qdot;
        wdot[13] -= qdot;
    }

    {
        // reaction 57:  CO + O+ => CO+ + O
        const amrex::Real k_f =
            2950848.9724 * pow(tc[1], 0.5) * exp(-4580.0 / tc[1]);
        const amrex::Real qf = k_f * (sc[5] * sc[12]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[5] -= qdot;
        wdot[12] -= qdot;
        wdot[13] += qdot;
    }

    {
        // reaction 58:  CO2+ + O => CO2 + O+
        const amrex::Real k_f = 57930000;
        const amrex::Real qf = k_f * (sc[1] * sc[10]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[4] += qdot;
        wdot[10] -= qdot;
        wdot[12] += qdot;
    }

    {
        // reaction 59:  CO + O- => CO2 + E
        const amrex::Real k_f = 331200000;
        const amrex::Real qf = k_f * (sc[5] * sc[11]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[4] += qdot;
        wdot[5] -= qdot;
        wdot[11] -= qdot;
    }

    {
        // reaction 60:  CO+ + O => CO + O+
        const amrex::Real k_f = 84310000;
        const amrex::Real qf = k_f * (sc[1] * sc[13]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[5] += qdot;
        wdot[12] += qdot;
        wdot[13] -= qdot;
    }

    {
        // reaction 61:  O + O- => E + O2
        const amrex::Real k_f = 138500000;
        const amrex::Real qf = k_f * (sc[1] * sc[11]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[1] -= qdot;
        wdot[2] += qdot;
        wdot[11] -= qdot;
    }

    {
        // reaction 62:  O+ + O- => 2 O
        const amrex::Real k_f = 24090000000;
        const amrex::Real qf = k_f * (sc[11] * sc[12]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += 2.000000 * qdot;
        wdot[11] -= qdot;
        wdot[12] -= qdot;
    }

    {
        // reaction 64:  CO2 + O => CO + O2
        const amrex::Real k_f = 16861994.128 * exp(-26500.0 / tc[1]);
        const amrex::Real qf = k_f * (sc[1] * sc[4]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[2] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
    }

    {
        // reaction 66:  CO + O2 => CO2 + O
        const amrex::Real k_f = 2529299.1192 * exp(-24000.0 / tc[1]);
        const amrex::Real qf = k_f * (sc[2] * sc[5]);
        const amrex::Real qr = 0.0;
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[2] -= qdot;
        wdot[4] += qdot;
        wdot[5] -= qdot;
    }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(
    const amrex::Real T,
    amrex::Real C[],
    amrex::Real wdot[],
    const amrex::Real Te,
    const amrex::Real EN,
    amrex::Real* ener_exch)
{

    productionRate(wdot, C, T, Te, ener_exch);
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRG(int kcharge[])
{
    kcharge[0] = -1;  // E
    kcharge[1] = 0;   // O
    kcharge[2] = 0;   // O2
    kcharge[3] = 0;   // AR
    kcharge[4] = 0;   // CO2
    kcharge[5] = 0;   // CO
    kcharge[6] = 0;   // CO2v1
    kcharge[7] = 0;   // CO2v2
    kcharge[8] = 0;   // CO2v3
    kcharge[9] = 0;   // CO2v4
    kcharge[10] = 1;  // CO2p
    kcharge[11] = -1; // Om
    kcharge[12] = 1;  // Op
    kcharge[13] = 1;  // COp
    kcharge[14] = 0;  // ARe
    kcharge[15] = 1;  // ARp
}

// compute electron temperature from EbyN using local field approximation
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
etempLFA(const amrex::Real Tg, amrex::Real C[], const amrex::Real EN)
{
    // 0.5 eV
    return (5802.25);
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRGMASS(amrex::Real zk[])
{

    int kchrg[16];
    CKCHRG(kchrg);

    for (int id = 0; id < 16; ++id)
    {
        zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
    }
}

#endif

#include "Chemistry.H"
const int rmap[68] = {5,  6,  10, 13, 15, 17, 63, 65, 67, 0,  1,  2,  3,  4,
                      7,  8,  9,  11, 12, 14, 16, 18, 19, 20, 21, 22, 23, 24,
                      25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
                      39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
                      53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 66};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 68; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[68] = {4, 4, 4, 4, 4, 3, 3, 5, 5, 4, 3, 4, 4, 3, 4, 3, 3,
                      3, 5, 5, 5, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4,
                      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 2};
  const int kiv[340] = {
    4,  0,  4,  0,  0, 4,  0,  6,  0,  0, 4,  0,  7, 0,  0, 4,  0,  8,  0, 0,
    4,  0,  9,  0,  0, 4,  5,  1,  0,  0, 4,  10, 0, 0,  0, 4,  0,  5,  0, 12,
    4,  0,  13, 0,  1, 4,  0,  5,  11, 0, 5,  13, 0, 0,  0, 3,  0,  3,  0, 0,
    3,  0,  14, 0,  0, 3,  15, 0,  0,  0, 10, 0,  5, 1,  0, 0,  12, 1,  0, 0,
    13, 0,  5,  0,  0, 0,  1,  11, 0,  0, 14, 4,  3, 5,  1, 15, 4,  3,  5, 12,
    15, 4,  3,  13, 1, 15, 4,  3,  10, 0, 4,  6,  4, 0,  0, 5,  6,  5,  4, 0,
    6,  2,  4,  2,  0, 3,  6,  3,  4,  0, 4,  7,  4, 0,  0, 5,  7,  5,  4, 0,
    7,  2,  4,  2,  0, 3,  7,  3,  4,  0, 4,  7,  4, 6,  0, 5,  7,  5,  6, 0,
    7,  2,  6,  2,  0, 3,  7,  3,  6,  0, 4,  8,  4, 7,  0, 5,  8,  5,  7, 0,
    8,  2,  7,  2,  0, 3,  8,  3,  7,  0, 4,  8,  4, 9,  0, 5,  8,  5,  9, 0,
    8,  2,  9,  2,  0, 3,  8,  3,  9,  0, 4,  8,  6, 7,  0, 4,  8,  4,  6, 0,
    5,  8,  5,  6,  0, 8,  2,  6,  2,  0, 3,  8,  3, 6,  0, 4,  9,  4,  7, 0,
    5,  9,  5,  7,  0, 9,  2,  7,  2,  0, 3,  9,  3, 7,  0, 4,  9,  4,  6, 0,
    5,  9,  5,  6,  0, 9,  2,  6,  2,  0, 3,  9,  3, 6,  0, 4,  12, 10, 1, 0,
    4,  13, 5,  10, 0, 5,  12, 13, 1,  0, 10, 1,  4, 12, 0, 5,  11, 4,  0, 0,
    13, 1,  5,  12, 0, 1,  11, 0,  2,  0, 11, 12, 1, 0,  0, 11, 12, 2,  0, 0,
    4,  1,  5,  2,  0, 5,  1,  4,  0,  0, 5,  2,  4, 1,  0, 1,  2,  0,  0, 0};
  const int nuv[340] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 2, 1,
    -1, -1, 1, 2, 1, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 68) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[16]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 16; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 68; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  // compute the Gibbs free energy
  amrex::Real g_RT[16];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 0.000549;  // E
  awt[1] = 15.999000; // O
  awt[2] = 12.011000; // C
  awt[3] = 39.950000; // Ar
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 4;
  // Zero ncf
  for (int id = 0; id < kd * 16; ++id) {
    ncf[id] = 0;
  }

  // E
  ncf[0 * kd + 0] = 1; // E

  // O
  ncf[1 * kd + 1] = 1; // O

  // O2
  ncf[2 * kd + 1] = 2; // O

  // AR
  ncf[3 * kd + 3] = 1; // Ar

  // CO2
  ncf[4 * kd + 2] = 1; // C
  ncf[4 * kd + 1] = 2; // O

  // CO
  ncf[5 * kd + 2] = 1; // C
  ncf[5 * kd + 1] = 1; // O

  // CO2v1
  ncf[6 * kd + 2] = 1; // C
  ncf[6 * kd + 1] = 2; // O

  // CO2v2
  ncf[7 * kd + 2] = 1; // C
  ncf[7 * kd + 1] = 2; // O

  // CO2v3
  ncf[8 * kd + 2] = 1; // C
  ncf[8 * kd + 1] = 2; // O

  // CO2v4
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 2; // O

  // CO2p
  ncf[10 * kd + 2] = 1;  // C
  ncf[10 * kd + 0] = -1; // E
  ncf[10 * kd + 1] = 2;  // O

  // Om
  ncf[11 * kd + 0] = 1; // E
  ncf[11 * kd + 1] = 1; // O

  // Op
  ncf[12 * kd + 0] = -1; // E
  ncf[12 * kd + 1] = 1;  // O

  // COp
  ncf[13 * kd + 2] = 1;  // C
  ncf[13 * kd + 0] = -1; // E
  ncf[13 * kd + 1] = 1;  // O

  // ARe
  ncf[14 * kd + 3] = 1; // Ar

  // ARp
  ncf[15 * kd + 3] = 1;  // Ar
  ncf[15 * kd + 0] = -1; // E
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "E";
  ename[1] = "O";
  ename[2] = "C";
  ename[3] = "Ar";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(16);
  kname[0] = "E";
  kname[1] = "O";
  kname[2] = "O2";
  kname[3] = "AR";
  kname[4] = "CO2";
  kname[5] = "CO";
  kname[6] = "CO2v1";
  kname[7] = "CO2v2";
  kname[8] = "CO2v3";
  kname[9] = "CO2v4";
  kname[10] = "CO2p";
  kname[11] = "Om";
  kname[12] = "Op";
  kname[13] = "COp";
  kname[14] = "ARe";
  kname[15] = "ARp";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 17; k++) {
    for (int l = 0; l < 17; l++) {
      if (Jac[17 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 17; k++) {
    for (int l = 0; l < 17; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[17 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 17; k++) {
    for (int l = 0; l < 17; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[17 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 17;
    int offset_col = nc * 17;
    for (int k = 0; k < 17; k++) {
      for (int l = 0; l < 17; l++) {
        if (Jac[17 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 17;
      for (int l = 0; l < 17; l++) {
        for (int k = 0; k < 17; k++) {
          if (Jac[17 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 17;
      for (int l = 0; l < 17; l++) {
        for (int k = 0; k < 17; k++) {
          if (Jac[17 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 17;
      for (int l = 0; l < 17; l++) {
        for (int k = 0; k < 17; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[17 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 17;
      for (int l = 0; l < 17; l++) {
        for (int k = 0; k < 17; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[17 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 17; k++) {
    for (int l = 0; l < 17; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 17 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[17 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 17 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 289> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 16> conc = {0.0};
  for (int n = 0; n < 16; n++) {
    conc[n] = 1.0 / 16.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 17; l++) {
      for (int k = 0; k < 17; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[17 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 17; l++) {
      for (int k = 0; k < 17; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[17 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}

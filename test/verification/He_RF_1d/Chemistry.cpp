#include "Chemistry.H"
const int rmap[9] = {2, 3, 4, 7, 8, 0, 1, 5, 6};

// Returns 0-based map of reaction order
void GET_RMAP(int *_rmap) {
  for (int j = 0; j < 9; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int &nspec, int ki[], int nu[]) {
  const int ns[9] = {4, 4, 3, 3, 3, 4, 4, 3, 3};
  const int kiv[36] = {0, 1, 0, 3, 0, 3, 0, 1, 1, 0, 2, 0, 3, 0, 2, 0, 5, 0,
                       4, 0, 0, 4, 1, 3, 3, 0, 1, 2, 1, 3, 5, 0, 1, 2, 4, 0};
  const int nuv[36] = {-1, -1, 1, 1, -1, -1, 1, 1, -1, 1,  1, 0,
                       -1, 1,  1, 0, -1, 1,  1, 0, -1, -1, 1, 1,
                       -2, 1,  1, 1, -1, -1, 1, 0, -1, -1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 4;
  } else {
    if (i > 9) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 4 + j] + 1;
        nu[j] = nuv[(i - 1) * 4 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[],
            amrex::Real q_f[], amrex::Real q_r[]) {
  amrex::Real c[6]; // temporary storage
  amrex::Real PORT =
      1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 6; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 9; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *q_f, amrex::Real *q_r, amrex::Real *sc,
                    amrex::Real T) {
  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  // compute the Gibbs free energy
  amrex::Real g_RT[6];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void atomicWeight(amrex::Real *awt) {
  awt[0] = 0.000549; // E
  awt[1] = 4.002602; // He
}

// get atomic weight for all elements
void CKAWT(amrex::Real *awt) { atomicWeight(awt); }

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void CKNCF(int *ncf) {
  int kd = 2;
  // Zero ncf
  for (int id = 0; id < kd * 6; ++id) {
    ncf[id] = 0;
  }

  // E
  ncf[0 * kd + 0] = 1; // E

  // HE
  ncf[1 * kd + 1] = 1; // He

  // HEp
  ncf[2 * kd + 0] = -1; // E
  ncf[2 * kd + 1] = 1;  // He

  // HEm
  ncf[3 * kd + 1] = 1; // He

  // HE2p
  ncf[4 * kd + 0] = -1; // E
  ncf[4 * kd + 1] = 2;  // He

  // HE2m
  ncf[5 * kd + 1] = 2; // He
}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string> &ename) {
  ename.resize(2);
  ename[0] = "E";
  ename[1] = "He";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string> &kname) {
  kname.resize(6);
  kname[0] = "E";
  kname[1] = "HE";
  kname[2] = "HEp";
  kname[3] = "HEm";
  kname[4] = "HE2p";
  kname[5] = "HE2m";
}

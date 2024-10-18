#include "Chemistry.H"
const int rmap[NUM_REACTIONS] = {3,5,6,7,10,11,12,13,15,0,1,2,4,8,9,14};

// Returns 0-based map of reaction order
void GET_RMAP
(int * _rmap)
{
for (int j=0; j<NUM_REACTIONS; ++j)
{
_rmap[j] = rmap[j];
}
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int& nspec, int ki[], int nu[])
{
const int ns[NUM_GAS_REACTIONS] =
     {4,4,4,3,4,2,3,3,4,4,3,3,3,3,5,3};
const int kiv[NUM_GAS_REACTIONS*5] =
     {0,20,0,8,0,0,20,0,9,0,0,20,0,10,0,20,0,1,0,0,0,19,0,6,0,19,11,0,0,0,19,11,7,0,0,19,0,2,0,0,0,19,11,3,0,0,17,15,3,0,17,14,16,0,0,17,15,11,0,0,17,0,5,0,0,18,14,12,0,0,0,18,0,14,13,18,0,4,0,0};
const int nuv[NUM_GAS_REACTIONS*5] =
     {-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,2,0,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,2,1,-1,1,1,0,0};
if (i < 1) {
// Return max num species per reaction
nspec = 5;
} else {
if (i > NUM_GAS_REACTIONS) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*5 + j] + 1;
nu[j] = nuv[(i-1)*5 + j];
}
}
}
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 0.000549; // E
awt[1] = 1.008000; // H
awt[2] = 14.007000; // N
awt[3] = 15.999000; // O
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int kd = 4; 
// Zero ncf
for (int id = 0; id < kd * 21; ++ id) {
 ncf[id] = 0; 
}

// E
ncf[ 0 * kd + 0 ] = 1; // E

// N2+
ncf[ 1 * kd + 0 ] = -1; // E
ncf[ 1 * kd + 2 ] = 2; // N

// O2+
ncf[ 2 * kd + 0 ] = -1; // E
ncf[ 2 * kd + 3 ] = 2; // O

// O-
ncf[ 3 * kd + 0 ] = 1; // E
ncf[ 3 * kd + 3 ] = 1; // O

// NH3+
ncf[ 4 * kd + 0 ] = -1; // E
ncf[ 4 * kd + 1 ] = 3; // H
ncf[ 4 * kd + 2 ] = 1; // N

// H2O+
ncf[ 5 * kd + 0 ] = -1; // E
ncf[ 5 * kd + 1 ] = 2; // H
ncf[ 5 * kd + 3 ] = 1; // O

// O2(B1)
ncf[ 6 * kd + 3 ] = 2; // O

// O(1D)
ncf[ 7 * kd + 3 ] = 1; // O

// N2(A3)
ncf[ 8 * kd + 2 ] = 2; // N

// N2(B3)
ncf[ 9 * kd + 2 ] = 2; // N

// N2(C3)
ncf[ 10 * kd + 2 ] = 2; // N

// O
ncf[ 11 * kd + 3 ] = 1; // O

// NH2
ncf[ 12 * kd + 1 ] = 2; // H
ncf[ 12 * kd + 2 ] = 1; // N

// NH
ncf[ 13 * kd + 1 ] = 1; // H
ncf[ 13 * kd + 2 ] = 1; // N

// H
ncf[ 14 * kd + 1 ] = 1; // H

// H2
ncf[ 15 * kd + 1 ] = 2; // H

// OH
ncf[ 16 * kd + 1 ] = 1; // H
ncf[ 16 * kd + 3 ] = 1; // O

// H2O
ncf[ 17 * kd + 1 ] = 2; // H
ncf[ 17 * kd + 3 ] = 1; // O

// NH3
ncf[ 18 * kd + 1 ] = 3; // H
ncf[ 18 * kd + 2 ] = 1; // N

// O2
ncf[ 19 * kd + 3 ] = 2; // O

// N2
ncf[ 20 * kd + 2 ] = 2; // N

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(4);
ename[0] = "E";
ename[1] = "H";
ename[2] = "N";
ename[3] = "O";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(21);
kname[0] = "E";
kname[1] = "N2+";
kname[2] = "O2+";
kname[3] = "O-";
kname[4] = "NH3+";
kname[5] = "H2O+";
kname[6] = "O2(B1)";
kname[7] = "O(1D)";
kname[8] = "N2(A3)";
kname[9] = "N2(B3)";
kname[10] = "N2(C3)";
kname[11] = "O";
kname[12] = "NH2";
kname[13] = "NH";
kname[14] = "H";
kname[15] = "H2";
kname[16] = "OH";
kname[17] = "H2O";
kname[18] = "NH3";
kname[19] = "O2";
kname[20] = "N2";
}

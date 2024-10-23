#include "Chemistry.H"
const int rmap[NUM_REACTIONS] = {0,1,4,5,7,10,11,12,14,16,17,26,29,30,31,37,39,64,65,66,2,3,6,8,9,13,15,18,19,20,21,22,23,24,25,27,28,32,33,34,35,36,38,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,67,68,69,70,71,72,73,74,75,76,77,78,79,80};

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
     {3,2,4,4,3,3,4,3,3,3,3,3,3,4,3,4,3,3,4,3,4,2,4,4,4,4,3,4,4,3,3,3,4,4,4,5,4,3,5,3,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,4,4,4,4,4,5,5,5,3,3,3,4,4,4,4,4,4,4,4,4,3,4,5,4,4};
const int kiv[NUM_GAS_REACTIONS*5] =
     {25,0,2,0,0,25,11,0,0,0,0,25,11,4,0,4,2,11,25,0,4,0,11,0,0,25,2,3,0,0,0,3,11,25,0,0,25,5,0,0,0,2,11,0,0,5,3,25,0,0,5,3,25,0,0,2,5,25,0,0,25,11,10,0,0,11,4,0,25,0,5,0,25,0,0,10,25,11,25,0,26,0,7,0,0,26,7,6,0,0,0,6,26,15,0,0,7,12,0,0,0,26,0,15,0,15,14,0,0,0,26,15,26,14,0,0,26,0,13,0,0,26,0,14,0,26,14,26,13,0,26,7,6,0,0,6,25,26,2,0,7,25,26,2,0,25,2,3,0,0,5,3,25,0,0,2,5,25,0,0,7,4,12,11,0,13,25,26,11,0,14,25,26,11,0,15,25,26,11,10,26,10,26,11,0,24,18,16,0,0,0,24,0,18,17,24,0,8,0,0,18,24,19,16,0,0,8,18,17,0,0,8,18,16,0,24,10,16,20,0,24,10,23,17,0,24,10,19,22,0,24,10,18,21,0,8,5,24,25,0,8,4,24,11,0,24,2,8,25,0,7,24,26,8,0,12,17,18,26,0,18,17,19,12,0,17,19,26,0,0,17,12,16,0,0,17,18,26,0,0,18,16,19,17,0,12,16,18,26,0,12,16,19,26,0,17,16,12,24,0,19,16,18,24,0,14,24,18,26,16,13,24,18,26,16,15,24,19,26,17,23,18,20,0,0,23,19,11,0,0,23,0,9,0,0,0,23,19,4,0,0,9,19,11,0,0,9,18,11,0,0,9,18,20,0,23,7,9,26,0,9,24,23,8,0,9,16,8,20,0,9,11,19,2,0,9,25,23,2,0,23,10,20,0,0,23,10,23,11,0,23,13,18,26,20,23,14,23,26,0,23,15,23,26,0};
const int nuv[NUM_GAS_REACTIONS*5] =
     {-1,1,1,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,2,1,0,-1,-1,1,0,0,-1,-1,2,0,0,-1,-1,3,0,0,-1,-1,3,0,0,-1,-1,2,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,3,0,0,-1,-1,2,0,0,-1,-1,2,1,0,-1,-1,1,2,0,-1,-1,1,2,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,2,1,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,1,0,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0};
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
for (int id = 0; id < kd * 27; ++ id) {
 ncf[id] = 0; 
}

// E
ncf[ 0 * kd + 0 ] = 1; // E

// O+
ncf[ 1 * kd + 0 ] = -1; // E
ncf[ 1 * kd + 3 ] = 1; // O

// O2+
ncf[ 2 * kd + 0 ] = -1; // E
ncf[ 2 * kd + 3 ] = 2; // O

// O4+
ncf[ 3 * kd + 0 ] = -1; // E
ncf[ 3 * kd + 3 ] = 4; // O

// O-
ncf[ 4 * kd + 0 ] = 1; // E
ncf[ 4 * kd + 3 ] = 1; // O

// O2-
ncf[ 5 * kd + 0 ] = 1; // E
ncf[ 5 * kd + 3 ] = 2; // O

// N4+
ncf[ 6 * kd + 0 ] = -1; // E
ncf[ 6 * kd + 2 ] = 4; // N

// N2+
ncf[ 7 * kd + 0 ] = -1; // E
ncf[ 7 * kd + 2 ] = 2; // N

// NH3+
ncf[ 8 * kd + 0 ] = -1; // E
ncf[ 8 * kd + 1 ] = 3; // H
ncf[ 8 * kd + 2 ] = 1; // N

// H2O+
ncf[ 9 * kd + 0 ] = -1; // E
ncf[ 9 * kd + 1 ] = 2; // H
ncf[ 9 * kd + 3 ] = 1; // O

// O(1D)
ncf[ 10 * kd + 3 ] = 1; // O

// O
ncf[ 11 * kd + 3 ] = 1; // O

// N
ncf[ 12 * kd + 2 ] = 1; // N

// N2(A1)
ncf[ 13 * kd + 2 ] = 2; // N

// N2(B3)
ncf[ 14 * kd + 2 ] = 2; // N

// N2(C3)
ncf[ 15 * kd + 2 ] = 2; // N

// NH2
ncf[ 16 * kd + 1 ] = 2; // H
ncf[ 16 * kd + 2 ] = 1; // N

// NH
ncf[ 17 * kd + 1 ] = 1; // H
ncf[ 17 * kd + 2 ] = 1; // N

// H
ncf[ 18 * kd + 1 ] = 1; // H

// H2
ncf[ 19 * kd + 1 ] = 2; // H

// OH
ncf[ 20 * kd + 1 ] = 1; // H
ncf[ 20 * kd + 3 ] = 1; // O

// H2NO
ncf[ 21 * kd + 1 ] = 2; // H
ncf[ 21 * kd + 2 ] = 1; // N
ncf[ 21 * kd + 3 ] = 1; // O

// NOH
ncf[ 22 * kd + 1 ] = 1; // H
ncf[ 22 * kd + 2 ] = 1; // N
ncf[ 22 * kd + 3 ] = 1; // O

// H2O
ncf[ 23 * kd + 1 ] = 2; // H
ncf[ 23 * kd + 3 ] = 1; // O

// NH3
ncf[ 24 * kd + 1 ] = 3; // H
ncf[ 24 * kd + 2 ] = 1; // N

// O2
ncf[ 25 * kd + 3 ] = 2; // O

// N2
ncf[ 26 * kd + 2 ] = 2; // N

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
kname.resize(27);
kname[0] = "E";
kname[1] = "O+";
kname[2] = "O2+";
kname[3] = "O4+";
kname[4] = "O-";
kname[5] = "O2-";
kname[6] = "N4+";
kname[7] = "N2+";
kname[8] = "NH3+";
kname[9] = "H2O+";
kname[10] = "O(1D)";
kname[11] = "O";
kname[12] = "N";
kname[13] = "N2(A1)";
kname[14] = "N2(B3)";
kname[15] = "N2(C3)";
kname[16] = "NH2";
kname[17] = "NH";
kname[18] = "H";
kname[19] = "H2";
kname[20] = "OH";
kname[21] = "H2NO";
kname[22] = "NOH";
kname[23] = "H2O";
kname[24] = "NH3";
kname[25] = "O2";
kname[26] = "N2";
}
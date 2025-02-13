#include "Chemistry.H"
const int rmap[NUM_REACTIONS] = {0};

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
     {3};
const int kiv[NUM_GAS_REACTIONS*3] =
     {1,0,2};
const int nuv[NUM_GAS_REACTIONS*3] =
     {-1,1,1};
if (i < 1) {
// Return max num species per reaction
nspec = 3;
} else {
if (i > NUM_GAS_REACTIONS) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*3 + j] + 1;
nu[j] = nuv[(i-1)*3 + j];
}
}
}
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 0.000549; // E
awt[1] = 14.007000; // N
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
int kd = 2; 
// Zero ncf
for (int id = 0; id < kd * 3; ++ id) {
 ncf[id] = 0; 
}

// E
ncf[ 0 * kd + 0 ] = 1; // E

// N2
ncf[ 1 * kd + 1 ] = 2; // N

// N2+
ncf[ 2 * kd + 0 ] = -1; // E
ncf[ 2 * kd + 1 ] = 2; // N

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(2);
ename[0] = "E";
ename[1] = "N";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(3);
kname[0] = "E";
kname[1] = "N2";
kname[2] = "N2+";
}

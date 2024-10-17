#include "Chemistry.H"
const int rmap[NUM_REACTIONS] = {1,2,0,3,4,5};

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
     {4,3,3,4,4,3};
const int kiv[NUM_GAS_REACTIONS*4] =
     {1,0,2,0,1,3,0,0,2,3,0,0,2,1,3,0,2,0,1,0,1,2,1,0};
const int nuv[NUM_GAS_REACTIONS*4] =
     {-1,-1,1,1,-1,1,1,0,-1,1,1,0,-2,1,1,1,-1,-1,1,1,-1,-1,2,0};
if (i < 1) {
// Return max num species per reproductionRateaction
nspec = 4;
} else {
if (i > NUM_GAS_REACTIONS) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*4 + j] + 1;
nu[j] = nuv[(i-1)*4 + j];
}
}
}
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 0.000549; // E
awt[1] = 39.950000; // Ar
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
for (int id = 0; id < kd * 4; ++ id) {
 ncf[id] = 0; 
}

// E
ncf[ 0 * kd + 0 ] = 1; // E

// AR
ncf[ 1 * kd + 1 ] = 1; // Ar

// ARm
ncf[ 2 * kd + 1 ] = 1; // Ar

// ARp
ncf[ 3 * kd + 1 ] = 1; // Ar
ncf[ 3 * kd + 0 ] = -1; // E

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(2);
ename[0] = "E";
ename[1] = "Ar";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(4);
kname[0] = "E";
kname[1] = "AR";
kname[2] = "ARm";
kname[3] = "ARp";
}

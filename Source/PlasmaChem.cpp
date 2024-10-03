#include <Prob.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <PlasmaChem.H>
#include <UnivConstants.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED int spec_chrg[NUM_SPECIES]={0};

    void init()
    {
        // Get vector with species names
        CKSYMS_STR(specnames);

        // Get species charge numbers
        CKCHRG(spec_chrg);
    }    

    void close(){}

    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }

    AMREX_GPU_HOST_DEVICE int get_charge(int specid)
    {
        return spec_chrg[specid];
    }
}

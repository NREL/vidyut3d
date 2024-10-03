#ifndef _REACTION_SOURCE_H_
#define _REACTION_SOURCE_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

    AMREX_GPU_DEVICE AMREX_INLINE
    void  compute_react_source(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& rxnsrc,
                              amrex::Real Tg, amrex::Real Pg)
    {
        amrex::Real rateconsts[NUM_REACTIONS]={0.0};
        amrex::Real specden[NUM_ALL_SPECIES]={0.0};
        amrex::Real concmult,reactrate;

        for(int sp=0;sp<(NUM_ALL_SPECIES);sp++)
        {
           specden[sp]=phi(i,j,k,sp); //no offset needed
        }
        amrex::Real efieldmag =std::sqrt(std::pow(phi(i,j,k,EFX_ID),2.0)+
                               std::pow(phi(i,j,k,EFY_ID),2.0)+
                               std::pow(phi(i,j,k,EFZ_ID),2.0));
        
        plasmachem::get_reaction_rateconstants(phi(i,j,k,ETEMP_ID), Tg, Pg, efieldmag, 
                                   specden,rateconsts);
        
        //include electron energy
        for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
        {
           rxnsrc(i,j,k,sp)=0.0;
        }

        for(int r=0;r<NUM_REACTIONS;r++)
        {
            concmult=1.0;
            for(int sp=0;sp<NUM_ALL_SPECIES;sp++)
            {
                if(plasmachem::rct_rxnarray[r][sp] > 0.0)
                {
                    concmult*=std::pow(phi(i,j,k,sp),plasmachem::rct_rxnarray[r][sp]);
                }
            }
            reactrate=rateconsts[r]*concmult;

            //update production/removal
            for(int sp=0;sp<NUM_ALL_SPECIES;sp++)
            {
                if(plasmachem::rct_rxnarray[r][sp] > 0.0)
                {
                    rxnsrc(i,j,k,sp) -= plasmachem::rct_rxnarray[r][sp]*reactrate;
                }
                if(plasmachem::pdt_rxnarray[r][sp] > 0.0)
                {
                    rxnsrc(i,j,k,sp) += plasmachem::pdt_rxnarray[r][sp]*reactrate;
                }
            }
            
            //for electron energy (EEN_ID is same as NUM_ALL_SPECIES)
            rxnsrc(i,j,k,EEN_ID) += plasmachem::reaction_elecenergy[r]*reactrate;
        }
    }
#endif

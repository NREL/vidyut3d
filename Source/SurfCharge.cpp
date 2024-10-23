#include <Vidyut.H>
#include <UserSources.H>
#include <Chemistry.H>
#include <UserFunctions.H>
#include <PlasmaChem.H>
#include <ProbParm.H>
#include <AMReX_MLABecLaplacian.H>
#include <HelperFuncs.H>

void Vidyut::update_surf_charge(Vector<MultiFab>& Sborder,
                                Real current_time, Real dt)
{
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> phi_arr = Sborder[ilev].array(mfi);
            Real time = current_time; // for GPU capture

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                //so the ghost cell index at left side is i-1 while it is i on the right
                if (bx.smallEnd(idim) == domain.smallEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryLo(bx, idim), 
                                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                                           dielectricflag=user_transport::is_dielectric(i, j, k, idim, -1, 
                                                                                        phi_arr,
                                                                                        prob_lo, prob_hi, dx, 
                                                                                        time); 
                                       }
                                       if(dielectricflag)
                                       {
                                           //find electron flux
                                           //find ion flux
                                           //do euler update of charge
                                       } 
                                       });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim))
            {
                amrex::ParallelFor(amrex::bdryHi(bx, idim), 
                                   [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                                       dielectricflag=user_transport::is_dielectric(i, j, k, idim, -1, 
                                                                                    phi_arr,
                                                                                    prob_lo, prob_hi, dx, 
                                                                                    time); 
                                       if(dielectricflag)
                                       {
                                           //find electron flux
                                           //find ion flux
                                           //do euler update of charge
                                       } 
                                   });
            }
        }
    }

}

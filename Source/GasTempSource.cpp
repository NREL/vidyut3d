#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <UserFunctions.H>
#include <compute_explicit_flux.H>
#include <AMReX_MLABecLaplacian.H>

void Vidyut::compute_gastemp_source(int lev, 
                            MultiFab& Sborder, 
                            MultiFab& rxnsrc, 
                            Array<MultiFab,AMREX_SPACEDIM>& efield, 
                            Array<MultiFab,AMREX_SPACEDIM>& gradne, 
                            MultiFab& dsdt,
                            Real time, Real dt, int floor_jh)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    amrex::Real gas_molar_density = (gas_pressure) / (Ru*gas_temperature); // density assumed to remain constant (thus, use intial value of temperature)

    
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();
        
    GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
    GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);

        Array4<Real> sborder_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);
        Array4<Real> phi_arr = phi_new[lev].array(mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            amrex::Real Cv_gas = 0.0;
            amrex::Real captured_gastemp = phi_arr(i,j,k,GASTEMP_ID);
            amrex::Real spec_C[NUM_SPECIES];

            for(int sp=0; sp<NUM_SPECIES; sp++) spec_C[sp] = sborder_arr(i,j,k,sp) / N_A;

            calcBulkGasCv(captured_gastemp, spec_C, &Cv_gas);
            dsdt_arr(i,j,k) = 0.3*phi_arr(i,j,k,EJH_ID)/(gas_molar_density*Cv_gas);

        });
    }
}

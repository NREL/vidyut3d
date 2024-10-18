#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>
#define PI 3.14159265358979323846264338327950288

// compute the integral of f dV, where f is one of the fields in the plotfile
// this understands axisymmetric geometry.

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    int dir=0;
    amrex::Real n0=1e6;
    amrex::Real fintime=0.1;

    int farg = 1;
    if (pltfile.empty() && farg <= narg) 
    {
        pltfile = amrex::get_command_argument(farg);
    }

    PlotFileData pf(pltfile);
    const Vector<std::string>& var_names_pf = pf.varNames();

    Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    Array<Real,AMREX_SPACEDIM> probhi = pf.probHi();
    const int dim = pf.spaceDim();
    int fine_level = pf.finestLevel();
    Vector<Real> pos;
    Real errorsum_e = 0.0;
    Real electronsum=0.0;
    int coord = pf.coordSys();

    Real domvol=AMREX_D_TERM((probhi[0]-problo[0]), 
                             *(probhi[1]-problo[1]), 
                             *(probhi[2]-problo[2]));

    for (int ilev = 0; ilev <= fine_level; ++ilev) 
    {

        Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(ilev);

        Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

        auto xlo = problo[dir];
        auto dxdir = dx[dir];

        if (ilev < fine_level) 
        {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) 
            {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);
            const MultiFab& mfe = pf.get(ilev, "E");
            auto const& ima = mask.const_arrays();
            auto const& mae = mfe.const_arrays();
            
            electronsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  return((ima[bno](i,j,k) == 0) ? mae[bno](i,j,k)*vol : 0._rt);
                              });
            

            errorsum_e += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real xlen=probhi[dir]-problo[dir];
                                  Real exactsoln=n0*(1.0+std::sin(PI*x/xlen))*std::exp(-5.0*fintime);
                                  return((ima[bno](i,j,k) == 0) ?  std::pow(mae[bno](i,j,k)-exactsoln,2.0)
                                         *vol : 0._rt);
                              });  
            
        } 
        else 
        {
            const MultiFab& mfe = pf.get(ilev, "E");
            const MultiFab& mfphi = pf.get(ilev, "Potential");
            auto const& mae = mfe.const_arrays();
            auto const& maphi = mfphi.const_arrays();

            electronsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  return(mae[bno](i,j,k)*vol);
                              });
            
            errorsum_e += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real xlen=probhi[dir]-problo[dir];
                                  Real exactsoln=n0*(1.0+std::sin(PI*x/xlen))*std::exp(-5.0*fintime);
                                  return(std::pow(mae[bno](i,j,k)-exactsoln,2.0)*vol);
                              });
            
        }
    }

    ParallelDescriptor::ReduceRealSum(errorsum_e);
    ParallelDescriptor::ReduceRealSum(electronsum);
    amrex::Real errnorm_e=std::sqrt(errorsum_e/domvol);
    amrex::Real electronavg=electronsum/domvol;
    amrex::Real errfrac_e=errnorm_e/electronavg;

    amrex::Print() << "L2 error average of E  = " << errnorm_e << '\n';
    amrex::Print() << "error frac = " << errfrac_e << '\n';
    amrex::Print() << "avg Electron = " << electronavg << '\n';

    if(errfrac_e > 1e-3)
    {
        amrex::Abort("Too much error in MMS3");
    }
    else
    {
        amrex::Print() << "MMS error less than 1e-3, error="<<errfrac_e<<"\n";
        amrex::Print()<<"Good to go!!\n";
    }

}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}

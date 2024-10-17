#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>

// compute the integral of f dV, where f is one of the fields in the plotfile
// this understands axisymmetric geometry.

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    int dir;
    amrex::Real n0=1e6;
    amrex::Real alpha=1.60217662e-19/8.854187817e-12;

    int farg = 1;
    while (farg <= narg) 
    {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-dir") 
        {
            dir = std::stoi(amrex::get_command_argument(++farg));
        } 
        else 
        {
            break;
        }
        ++farg;
    }

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
    Real potsum=0.0;
    Real errorsum_pot=0.0;
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
            const MultiFab& mfphi = pf.get(ilev, "Potential");
            auto const& ima = mask.const_arrays();
            auto const& maphi = mfphi.const_arrays();
            auto const& mae = mfe.const_arrays();
            
            electronsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  return((ima[bno](i,j,k) == 0) ? mae[bno](i,j,k)*vol : 0._rt);
                              });
            
            potsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  return((ima[bno](i,j,k) == 0) ? maphi[bno](i,j,k)*vol : 0._rt);
                              });

            errorsum_e += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real exactsoln=x*x*x/alpha+n0;
                                  return((ima[bno](i,j,k) == 0) ?  std::pow(mae[bno](i,j,k)-exactsoln,2.0)
                                         *vol : 0._rt);
                              });  
            
            errorsum_pot += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real exactsoln=(std::pow(x,5.0)-x)/40.0;
                                  return((ima[bno](i,j,k) == 0) ?  std::pow(maphi[bno](i,j,k)-exactsoln,2.0)
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
            
            potsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  return(maphi[bno](i,j,k)*vol);
                              });

            errorsum_e += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real exactsoln=x*x*x/alpha+n0;
                                  return(std::pow(mae[bno](i,j,k)-exactsoln,2.0)*vol);
                              });
            
            errorsum_pot += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mfe,
                              [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                              -> GpuTuple<Real> 
                              {
                                  IntVect iv(i,j,k);
                                  Real x=xlo+(iv[dir]+0.5)*dxdir;
                                  Real exactsoln=(std::pow(x,5.0)-x)/40.0;
                                  return(std::pow(maphi[bno](i,j,k)-exactsoln,2.0)*vol);
                              });
        }
    }

    ParallelDescriptor::ReduceRealSum(errorsum_e);
    ParallelDescriptor::ReduceRealSum(errorsum_pot);
    ParallelDescriptor::ReduceRealSum(electronsum);
    ParallelDescriptor::ReduceRealSum(potsum);
    amrex::Real errnorm_e=std::sqrt(errorsum_e/domvol);
    amrex::Real errnorm_pot=std::sqrt(errorsum_pot/domvol);
    amrex::Real electronavg=electronsum/domvol;
    amrex::Real potavg=potsum/domvol;
    amrex::Real errfrac_e=errnorm_e/electronavg;
    amrex::Real errfrac_pot=amrex::Math::abs(errnorm_pot)/amrex::Math::abs(potavg);

    amrex::Print() << "L2 error average of E, pot  = " << errnorm_e << '\t' << errnorm_pot <<"\n";
    amrex::Print() << "error frac = " << errfrac_e << '\t' << errfrac_pot<<"\n";
    amrex::Print() << "Electron/pot avg = " << electronavg << '\t' << potavg << "\n";

    if(errfrac_e > 1e-2)
    {
        amrex::Abort("Too much error");
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

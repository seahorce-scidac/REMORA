#include "prob.H"
#include "REMORA_prob_common.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

void Problem::init_analytic_vmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);

      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      Gpu::streamSynchronize();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 1.0e-5_rt;

        Akt(i,j,k,Temp_comp) = 1.0e-6_rt;
        Akt(i,j,k,Salt_comp) = 1.0e-6_rt;
        Akt(i,j,k,Tracer_comp) = 0.0_rt;
      });
    }
}

void Problem::init_analytic_hmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_visc2_p,
        MultiFab& mf_visc2_r,
        MultiFab& mf_diff2)
{
    for ( MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& visc2_p_arr = mf_visc2_p.array(mfi);
      Array4<Real> const& visc2_r_arr = mf_visc2_r.array(mfi);
      Array4<Real> const& diff2_arr   = mf_diff2.array(mfi);

      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      int ncomp = mf_diff2.nComp(); // temperature and salt and tracer
      Gpu::streamSynchronize();

      amrex::ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
          visc2_p_arr(i,j,0) = 5.0e0;
          visc2_r_arr(i,j,0) = 5.0e0;

          for (int n = 0; n < ncomp; n++) {
              diff2_arr(i,j,0,n) = 5.0e3;
          }
      });
    } // mfi
}

void Problem::init_analytic_smflux(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    mf_sustr.setVal(0.0);
    mf_svstr.setVal(0.0);
}

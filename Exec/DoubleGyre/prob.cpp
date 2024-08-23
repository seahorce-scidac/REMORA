#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"
#include "DepthStretchTransform.H"

using namespace amrex;

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void
init_custom_bathymetry (int /*lev*/, const Geometry& /*geom*/,
                        MultiFab& mf_h,
                        const SolverChoice& /*m_solverChoice*/,
                        int /*rrx*/, int /*rry*/)
{
    mf_h.setVal(500.0_rt);
}

/**
 * \brief Initializes coriolis factor
 */
void
init_custom_coriolis    (const Geometry& /*geom*/,
                         MultiFab& /*mf_fcor*/,
                         const SolverChoice& /*m_solverChoice*/) {}

/**
 * \brief Initializes custom sea surface height
 */
void
init_custom_zeta (const Geometry& geom,
                      MultiFab& mf_zeta,
                      const SolverChoice& m_solverChoice)
{
    mf_zeta.setVal(0.0_rt);
}

void
init_custom_prob(
        const Box& bx,
        Array4<Real      > const& state,
        Array4<Real      > const& x_vel,
        Array4<Real      > const& y_vel,
        Array4<Real      > const& z_vel,
        Array4<Real const> const& /*z_w*/,
        Array4<Real const> const& z_r,
        Array4<Real const> const& /*Hz*/,
        Array4<Real const> const& /*h*/,
        Array4<Real const> const& /*Zt_avg1*/,
        GeometryData const& geomdata,
        const SolverChoice& m_solverChoice)
{
    bool l_use_salt = m_solverChoice.use_salt;

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    auto T0 = m_solverChoice.T0;
    Real val1 = (44.69_rt / 39.382_rt) * (44.69_rt / 39.382_rt);
    Real val2 = val1 * (m_solverChoice.rho0 * 100.0_rt/m_solverChoice.g) * (5.0e-5_rt/((42.689_rt/44.69_rt) * (42.689_rt/44.69_rt)));
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo         = geomdata.ProbLo();
        const auto prob_hi         = geomdata.ProbHi();
        const auto dx              = geomdata.CellSize();

        const Real y  = prob_lo[1] + (j + 0.5_rt) * dx[1];
        const Real z = z_r(i,j,k);
        const Real yextent = prob_hi[1] - prob_lo[1];

        const Real val3 = T0 + val2 * std::exp(z/100.0_rt) * (10.0_rt - 0.4_rt * std::tanh(z / 100.0_rt));
        const Real val4 = y / yextent;

        state(i,j,k,Temp_comp)=val3 - 3.0_rt * val4;
        if (l_use_salt) {
            state(i,j,k,Salt_comp)=34.5_rt - 0.001_rt * z - val4;
        }

        // Set scalar = 0 everywhere
        state(i, j, k, Scalar_comp) = 0.0_rt;
    });

  const Box& xbx = surroundingNodes(bx,0);
  const Box& ybx = surroundingNodes(bx,1);
  const Box& zbx = surroundingNodes(bx,2);

  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i,j,k) = 0.0_rt;
  });
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0_rt;
  });

  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0_rt;
  });

  Gpu::streamSynchronize();
}

void
init_custom_vmix(const Geometry& /*geom*/, MultiFab& mf_Akv, MultiFab& mf_Akt,
                 MultiFab& mf_z_w, const SolverChoice& /*m_solverChoice*/)
{
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Array4<Real> const& z_w = (mf_z_w).array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 1.0_rt; //2.0e-03_rt+8.0e-03_rt*std::exp(z_w(i,j,k)/150.0_rt);

        Akt(i,j,k,Temp_comp) = 1.0_rt;
        Akt(i,j,k,Salt_comp) = 1.0_rt;
        Akt(i,j,k,Scalar_comp) = 0.0_rt;
      });
    }
}

void
init_custom_hmix(const Geometry& /*geom*/, MultiFab& mf_visc2_p, MultiFab& mf_visc2_r,
                 MultiFab& mf_diff2, const SolverChoice& /*m_solverChoice*/)
{
    for ( MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& visc2_p = (mf_visc2_p).array(mfi);
      Array4<Real> const& visc2_r = (mf_visc2_r).array(mfi);
      Array4<Real> const& diff2   = mf_diff2.array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();

      int ncomp = mf_diff2.nComp();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        visc2_p(i,j,k) = 1280.0_rt;
        visc2_r(i,j,k) = 1280.0_rt;

        for (int n = 0; n < ncomp; n++) {
            diff2(i,j,k,n) = 1280.0_rt;
        }
      });
    }
}

void
init_custom_smflux(const Geometry& geom, const Real time, MultiFab& mf_sustr, MultiFab& mf_svstr,
                   const SolverChoice& m_solverChoice)
{
    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    const auto prob_lo         = geomdata.ProbLo();
    const auto prob_hi         = geomdata.ProbHi();

    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359_rt;
    Real tdays=time/Real(24.0*60.0*60.0);
    Real dstart=0.0_rt;
    Real air_rho=1.0_rt;

    const Real yextent = prob_hi[1] - prob_lo[1];
    const Real windamp = -0.05_rt / m_solverChoice.rho0;
    const Real val1  = 2.0_rt * pi / yextent;
    for ( MFIter mfi((mf_sustr), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        const Box& xbx = surroundingNodes(bx,0);
        const Box& xbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));

        Array4<Real> const& sustr = mf_sustr.array(mfi);
        ParallelFor(xbx2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Create bounding box for x and y to make spatially-dependent T and S
            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            const Real y  = prob_lo[1] + (j + 0.5) * dx[1];// - ycent;

            sustr(i,j,0) = windamp * std::cos(val1 * y);
        });
    }
    mf_svstr.setVal(0.0_rt);
}

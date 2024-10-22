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
init_custom_bathymetry (int /*lev*/, const Geometry& geom,
                        MultiFab& mf_h,
                        const SolverChoice& m_solverChoice,
                        int /*rrx*/, int /*rry*/)
{
    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    // mf_h.setVal(geomdata.ProbHi(2));
    mf_h.setVal(50.0_rt);

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& h  = (mf_h).array(mfi);

        Box bx = mfi.tilebox();

        Box bxD = bx;
        bxD.makeSlab(2,0);

        Gpu::streamSynchronize();

        ParallelForRNG(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int , auto engine)
        {
            h(i,j,0) = 50.0_rt;
        });

    } // mfi
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
        Array4<Real const> const& h,
        Array4<Real const> const& /*Zt_avg1*/,
        GeometryData const& geomdata,
        const SolverChoice& m_solverChoice)
{
    bool l_use_salt = m_solverChoice.use_salt;

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    bool EWPeriodic = geomdata.isPeriodic(1);
    bool NSPeriodic = geomdata.isPeriodic(0);

    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;
    ParallelFor(grow(grow(bx,0,1),1,1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo         = geomdata.ProbLo();
        const auto prob_hi         = geomdata.ProbHi();
        const auto dx              = geomdata.CellSize();

        const Real z = z_r(i,j,k);

        // Create bounding box for x and y to make spatially-dependent T and S
        const Real xcent = 0.5*(prob_lo[0] + prob_hi[0]);
        const Real ycent = 0.5*(prob_lo[1] + prob_hi[1]);

        const Real x  = prob_lo[0] + (i + 0.5) * dx[0] - xcent;
        const Real y  = prob_lo[1] + (j + 0.5) * dx[1];

        state(i,j,k,Temp_comp)=T0 + z / (9.8_rt * 1.7_rt);
        if (l_use_salt) {
            state(i,j,k,Salt_comp)=S0 + y / (9.8_rt * 760.0_rt);
        }

        // Set scalar = 0 everywhere
        state(i, j, k, Scalar_comp) = 0.0_rt;
    });

    const Box& xbx = surroundingNodes(bx,0);
    const Box& ybx = surroundingNodes(bx,1);
    const Box& zbx = surroundingNodes(bx,2);

    ParallelFor(grow(grow(xbx,1,1),0,1), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Create bounding box for x and y to make spatially-dependent T and S
        const auto prob_lo         = geomdata.ProbLo();
        const auto prob_hi         = geomdata.ProbHi();
        const auto dx              = geomdata.CellSize();

        const Real xcent = 0.5*(prob_lo[0] + prob_hi[0]);
        const Real ycent = 0.5*(prob_lo[1] + prob_hi[1]);

        const Real z = z_r(i,j,k);
        // h(i,j,0) = -geomdata.ProbLo(2);
        const Real depth = 50.0_rt;

        const Real x  = prob_lo[0] + (i + 0.5) * dx[0];// - xcent;
        const Real y  = prob_lo[1] + (j + 0.5) * dx[1];// - ycent;

        const Real val1  = x * 6.28318530718_rt * 1.0e-5_rt;
        const Real val2  = y * 6.28318530718_rt * 1.0e-5_rt;
        x_vel(i, j, k) = (depth + z) * 0.01_rt + 0.05_rt * std::sin(val1) * std::sin(val2);
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
            Akv(i,j,k) = 2.0e-03_rt+8.0e-03_rt*std::exp(z_w(i,j,k)/150.0_rt);

            Akt(i,j,k,Temp_comp) = 1.0e-6_rt;
            Akt(i,j,k,Salt_comp) = 1.0e-6_rt;
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
          visc2_p(i,j,k) = 2.0_rt;
          visc2_r(i,j,k) = 2.0_rt;

          for (int n = 0; n < ncomp; n++) {
              diff2(i,j,k,n) = 0.0_rt;
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

    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359_rt;
    Real tdays=time/Real(24.0*60.0*60.0);
    Real dstart=0.0_rt;
    Real windamp;
    Real userval=0.01_rt;
    Real mf_fcor=1.0e-4_rt;
    Real air_rho=1.0_rt;
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

    // Flow in x-direction (EW):
    if (NSPeriodic) {
        mf_sustr.setVal(0.0_rt);
    }
    else if (EWPeriodic) {
        windamp=userval*std::sin(0.9_rt*mf_fcor*86400.0_rt*(tdays-dstart))/m_solverChoice.rho0;
        mf_sustr.setVal(windamp);
    }

    // Flow in y-direction (NS):
    if (NSPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1_rt*Real(sin(pi*(tdays-dstart)/4.0_rt))/Real(m_solverChoice.rho0);
        else
            windamp=-0.1_rt/m_solverChoice.rho0;
        mf_svstr.setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr.setVal(0.0_rt);
    }
}

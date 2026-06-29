#include "prob.H"
#include "REMORA_prob_common.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_h)
{
    mf_h.setVal(500.0_rt);
}

/**
 * \brief Initializes custom sea surface height
 */
void Problem::init_analytic_zeta (
        int /*lev*/, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_zeta)
{
    mf_zeta.setVal(0.0_rt);
}

void Problem::init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel)
{
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    auto T0 = m_solverChoice.T0;
    Real val1 = (44.69_rt / 39.382_rt) * (44.69_rt / 39.382_rt);
    Real val2 = val1 * (m_solverChoice.rho0 * 100.0_rt/g) * (5.0e-5_rt/((42.689_rt/44.69_rt) * (42.689_rt/44.69_rt)));

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

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
            state(i, j, k, Tracer_comp) = 0.0_rt;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i,j,k) = 0.0_rt;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = 0.0_rt;
        });
    }
    Gpu::streamSynchronize();
}

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
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 1.0_rt;

        Akt(i,j,k,Temp_comp) = 1.0_rt;
        Akt(i,j,k,Salt_comp) = 1.0_rt;
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
      Array4<Real> const& visc2_p = (mf_visc2_p).array(mfi);
      Array4<Real> const& visc2_r = (mf_visc2_r).array(mfi);
      Array4<Real> const& diff2   = mf_diff2.array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();

      int ncomp = mf_diff2.nComp();

      amrex::ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
        visc2_p(i,j,0) = 1280.0_rt;
        visc2_r(i,j,0) = 1280.0_rt;

        for (int n = 0; n < ncomp; n++) {
            diff2(i,j,0,n) = 1280.0_rt;
        }
      });
    }
}

void Problem::init_analytic_smflux(
        int /*lev*/,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    auto geomdata = geom.data();
    [[maybe_unused]] bool EWPeriodic = geomdata.isPeriodic(0);
    [[maybe_unused]] bool NSPeriodic = geomdata.isPeriodic(1);

    const auto prob_lo         = geomdata.ProbLo();
    const auto prob_hi         = geomdata.ProbHi();

    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359_rt;

    const Real yextent = prob_hi[1] - prob_lo[1];
    const Real windamp = -0.05_rt / m_solverChoice.rho0;
    const Real val1  = 2.0_rt * pi / yextent;
    for ( MFIter mfi((mf_sustr), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        [[maybe_unused]] const Box& xbx = surroundingNodes(bx,0);
        const Box& xbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));

        Array4<Real> const& sustr = mf_sustr.array(mfi);
        ParallelFor(xbx2, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
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

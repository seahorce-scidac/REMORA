#include "prob.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h)
{
    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    const int Lm = geomdata.Domain().size()[0];
    const int Mm = geomdata.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);
      Array4<const Real> const& xr  = remora.vec_xr[lev]->array(mfi);
      Array4<const Real> const& yr  = remora.vec_yr[lev]->array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW,NGROW,0));

      Box gbx2D = gbx2;
      gbx2D.makeSlab(2,0);

      Gpu::streamSynchronize();

      ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
          Real val1 = (xr(i,j,0) + 500.0_rt) / 15000.0_rt;
          h(i,j,0) = 14.0_rt + 25.0_rt * (1.0_rt - std::exp(-PI * xr(i,j,0) * 1.0e-5_rt)) - 8.0_rt * std::exp(-val1 * val1);
      });

    } // mfi
}

void Problem::init_analytic_grid_scale (
        int /*lev*/, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_pm, amrex::MultiFab& mf_pn)
{
    auto geomdata = geom.data();
    const int Lm = geomdata.Domain().size()[0];
    const auto dxi = geom.InvCellSize();
    for ( MFIter mfi(mf_pm, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& pm = (mf_pm).array(mfi);

        Box bx = mfi.tilebox();
        Box gbx2 = bx;
        gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Gpu::streamSynchronize();
        ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real dx = 0.5_rt * (4000.0_rt / (Lm + 1)) * (i+1) + 675.0_rt;
            pm(i,j,0) = 1.0_rt / dx;
        });
    }
    mf_pn.setVal(dxi[1]);
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
        amrex::MultiFab& mf_yvel,
        amrex::MultiFab& mf_zvel)
{
    bool l_use_salt = m_solverChoice.use_salt;
    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    auto geomdata = geom.data();

    const int khi = geomdata.Domain().bigEnd()[2];

    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);
        Array4<      Real> const& z_vel = mf_zvel.array(mfi);

        Array4<const Real> const& z_r   = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(grow(bx,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real z = z_r(i,j,k);
            Real val1 = std::tanh(1.1_rt * z + 11.0_rt);
            state(i,j,k,Temp_comp)=T0+6.25_rt*val1;
            state(i,j,k,Salt_comp)=S0-0.75_rt*val1;

            state(i, j, k, Scalar_comp) = 0.0_rt;
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        const Box& zbx = surroundingNodes(bx,2);

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            x_vel(i, j, k) = 0.0_rt;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = 0.0_rt;
        });

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel(i, j, k) = 0.0_rt;
        });
    }
    Gpu::streamSynchronize();
}

void Problem::init_analytic_vmix(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);

      Array4<const Real> const& z_w = remora.vec_z_w[lev]->const_array(mfi);
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

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        visc2_p(i,j,k) = 5.0_rt;
        visc2_r(i,j,k) = 5.0_rt;

        for (int n = 0; n < ncomp; n++) {
            diff2(i,j,k,n) = 0.0_rt;
        }
      });
    }
}

void Problem::init_analytic_wind(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_Uwind, MultiFab& mf_Vwind)
{
    Real time = remora.get_t_old(lev);
    Real tdays=time/Real(24.0*60.0*60.0);
    Real dstart=0.0_rt;

    mf_Uwind.setVal(0.0_rt);
    if (tdays - dstart <= 6.0_rt) {
        mf_Vwind.setVal(10.0_rt);
    } else {
        mf_Vwind.setVal(0.0_rt);
    }
}

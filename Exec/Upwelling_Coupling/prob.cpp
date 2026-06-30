#include "prob.H"

#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "REMORA_IndexDefines.H"
#include "REMORA_DepthStretchTransform.H"

using namespace amrex;

ProbParm parms;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{
    // Problem-specific parameters can be parsed here if needed later.
}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int /*lev*/, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_h)
{
    auto geomdata = geom.data();

    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));
    Print() << "Initializing coastal bathymetry" << std::endl;
    for (MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& h = mf_h.array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = bx;
        gbx.grow(mf_h.nGrowVect());

        Box gbxD = gbx;
        gbxD.makeSlab(2,0);

        Gpu::streamSynchronize();

        ParallelFor(gbxD, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const auto prob_lo = geomdata.ProbLo();
            const auto dx = geomdata.CellSize();

            const Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const Real x0 = prob_lo[0];
            const Real x1 = geomdata.ProbHi(0);
            const Real h0 = 300.0_rt;
            const Real h1 = 10.0_rt;

            Real depth = h0 + (h1 - h0) * (x - x0) / (x1 - x0);
            depth = amrex::min(amrex::max(depth, h1), h0);
            h(i,j,0) = depth;
        });
    }
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
    const bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    [[maybe_unused]] bool EWPeriodic = geomdata.isPeriodic(0);
    [[maybe_unused]] bool NSPeriodic = geomdata.isPeriodic(1);

    auto S0 = m_solverChoice.S0;
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
    	amrex::Print() << "Box is " << bx << "comparing to khi+1 " << khi+1 << std::endl;
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<Real> const& state = mf_cons.array(mfi);
        Array4<Real> const& x_vel = mf_xvel.array(mfi);
        Array4<Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Convert model z coordinate to positive depth below the surface.
            const Real depth = amrex::max(0.0_rt, -z_r(i, j, k));

            Real temp = 6.5_rt;
            if (depth <= 20.0_rt) {
                temp = 14.0_rt;
            } else if (depth <= 50.0_rt) {
                const Real alpha = (depth - 20.0_rt) / (50.0_rt - 20.0_rt);
                temp = 14.0_rt + alpha * (10.0_rt - 14.0_rt);
            } else if (depth <= 300.0_rt) {
                const Real alpha = (depth - 50.0_rt) / (300.0_rt - 50.0_rt);
                temp = 10.0_rt + alpha * (6.5_rt - 10.0_rt);
            }

            Real salt = 33.5_rt;
            if (depth <= 50.0_rt) {
                salt = 33.5_rt;
            } else if (depth <= 150.0_rt) {
                const Real alpha = (depth - 50.0_rt) / (150.0_rt - 50.0_rt);
                salt = 33.5_rt + alpha * (34.3_rt - 33.5_rt);
            } else {
                const Real alpha = amrex::min((depth - 150.0_rt) / (300.0_rt - 150.0_rt), 1.0_rt);
                salt = 34.3_rt + alpha * (34.6_rt - 34.3_rt);
            }

            state(i, j, k, Temp_comp) = temp;
            if (l_use_salt) {
                state(i, j, k, Salt_comp) = salt;
            }

            state(i, j, k, Tracer_comp) = 0.0_rt;
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
    for (MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& Akv = mf_Akv.array(mfi);
        Array4<Real> const& Akt = mf_Akt.array(mfi);
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
    for (MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& visc2_p = mf_visc2_p.array(mfi);
        Array4<Real> const& visc2_r = mf_visc2_r.array(mfi);
        Array4<Real> const& diff2   = mf_diff2.array(mfi);
        Box bx = mfi.tilebox();
        bx.grow(IntVect(NGROW,NGROW,0));
        Gpu::streamSynchronize();

        int ncomp = mf_diff2.nComp();

        amrex::ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            visc2_p(i,j,0) = 2.0_rt;
            visc2_r(i,j,0) = 2.0_rt;

            for (int n = 0; n < ncomp; n++) {
                diff2(i,j,0,n) = 0.0_rt;
            }
        });
    }
}

void Problem::init_analytic_smflux(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    for (MFIter mfi((mf_sustr), TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Array4<Real> const& sustr = (mf_sustr).array(mfi);
      Array4<Real> const& svstr = (mf_svstr).array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      amrex::ParallelFor(makeSlab(bx,2,0),
      [=] AMREX_GPU_DEVICE (int i, int j, int)
      {
        sustr(i,j,0) = 0.0_rt;
        svstr(i,j,0) = 0.0_rt;
      });
    }
}

#include "prob.H"
#include "REMORA_prob_common.H"

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
    // Parse params
    ParmParse pp("remora.prob");

    pp.query("u_0", parms.u_0);
    pp.query("v_0", parms.v_0);
}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int /*lev*/, const amrex::Geometry& geom,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_h)
{
    const auto & geomdata = geom.data();
    mf_h.setVal(geomdata.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Gpu::streamSynchronize();
      amrex::ParallelFor(Box(h),
      [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
          h(i,j,0,0) = -geomdata.ProbLo(2);
          if (k==0) {
              h(i,j,0,1) = h(i,j,0,0);
          }
      });
    } // mfi
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
        int /*lev*/,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel,
        amrex::MultiFab& mf_zvel)
{
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);
        Array4<      Real> const& z_vel = mf_zvel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();

            state(i, j, k, Temp_comp) = 1.;

            state(i,j,k,Temp_comp)=T0; //+8.0*std::exp(z/50.0_rt);

            // Set tracer = 0 everywhere
            const Real xcent = 0.5*(prob_lo[0] + prob_hi[0]);
            const Real ycent = 0.5*(prob_lo[1] + prob_hi[1]);

            const Real x  = prob_lo[0] + (i + 0.5) * dx[0] - xcent;
            const Real y  = prob_lo[1] + (j + 0.5) * dx[1] - ycent;
            const Real r2 = x*x + y*y;
            const Real rad = 0.1 * (prob_hi[0]-prob_lo[0]);
            const Real radsq = rad*rad;
            const Real rad_inner = 0.05 * (prob_hi[0]-prob_lo[0]);
            [[maybe_unused]] const Real rad_inner_sq = rad_inner*rad_inner;

            if (l_use_salt) {
                state(i,j,k,Salt_comp)= S0;
            }

            // Single circle of tracer (default)
            state(i, j, k, Tracer_comp) = std::exp(-r2/(2.*radsq));

            // Donut of tracer
            //state(i, j, k, Tracer_comp) = 1.25 * (std::exp(-r2/(2.*radsq)) - std::exp(-r2/(2*rad_inner_sq)));
        });

        // Construct a box that is on x-faces
        const Box& xbx = surroundingNodes(bx,0);
        // Set the x-velocity
        ParallelFor(xbx, [=, parms_gpu=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              x_vel(i, j, k) = parms_gpu.u_0;
        });

        // Construct a box that is on y-faces
        const Box& ybx = surroundingNodes(bx,1);

        // Set the y-velocity
        ParallelFor(ybx, [=, parms_gpu=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              y_vel(i, j, k) = parms_gpu.v_0;
        });

        // Construct a box that is on z-faces
        const Box& zbx = surroundingNodes(bx,2);

        // Set the z-velocity
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel(i, j, k) = 0.0;
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

      Array4<const Real> const& z_w = remora.vec_z_w[lev]->array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = 2.0e-03+8.0e-03*std::exp(z_w(i,j,k)/150.0);

        Akt(i,j,k,Temp_comp) = 1.0e-6;
        Akt(i,j,k,Salt_comp) = 1.0e-6;
        Akt(i,j,k,Tracer_comp) = 0.0;
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
        visc2_p(i,j,0) = 5.0;
        visc2_r(i,j,0) = 5.0;

        for (int n = 0; n < ncomp; n++) {
            diff2(i,j,0,n) = 0.0;
        }
      });
    }
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

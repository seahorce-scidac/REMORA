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
        int /*lev*/, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        amrex::MultiFab& mf_h)
{
    mf_h.setVal(geom.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = grow(bx,IntVect(NGROW+1,NGROW+1,0));

      const auto & geomdata = geom.data();

      amrex::Real Xsize = 320000.0_rt;
      amrex::Real Esize = 320000.0_rt;
      amrex::Real depth = 5000.0_rt;

      if(!m_solverChoice.flat_bathymetry)
      {
          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
          {
              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              const Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
              const Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];

              Real val1, val2;
              val1 = (x-0.5_rt*Xsize)/40000.0_rt;
              val2 = (y-0.5_rt*Esize)/40000.0_rt;
              h(i,j,0) = depth - 4500.0_rt * std::exp(-(val1*val1+val2*val2));
          });

      } else {

          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              h(i,j,0) = -geomdata.ProbLo(2);
              h(i,j,0,1) = h(i,j,0);
          });
      }
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

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            // const auto prob_lo         = geomdata.ProbLo();
            // const auto dx              = geomdata.CellSize();

            // const Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            // const Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const Real z = z_r(i,j,k);

            state(i, j, k, Temp_comp) = 1.;

            state(i,j,k,Temp_comp)=T0+7.5_rt*std::exp(z/1000.0_rt);
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=S0;
            }

            // Set tracer = 0 everywhere
            state(i, j, k, Tracer_comp) = 0.0;
        });

        // Construct a box that is on x-faces
        const Box& xbx = surroundingNodes(bx,0);
        // Set the x-velocity
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              // Set the x-velocity
              x_vel(i, j, k) = 0.0;
        });

        // Construct a box that is on y-faces
        const Box& ybx = surroundingNodes(bx,1);

        // Set the y-velocity
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              y_vel(i, j, k) = 0.0;
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
          visc2_p_arr(i,j,0) = 0.0;
          visc2_r_arr(i,j,0) = 0.0;

          for (int n = 0; n < ncomp; n++) {
              diff2_arr(i,j,0,n) = 0.0;
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

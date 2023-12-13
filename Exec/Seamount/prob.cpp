#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"
#include "DepthStretchTransform.H"

using namespace amrex;

ProbParm parms;

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  ParmParse pp("prob");
}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void
init_custom_bathymetry (const Geometry& geom,
                        MultiFab& mf_h,
                        MultiFab& mf_Zt_avg1,
                        const SolverChoice& m_solverChoice)
{
    mf_h.setVal(geom.ProbHi(2));

    //HACK HACK manually setting zeta to 0
    mf_Zt_avg1.setVal(0.0);

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

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

              const Real x = prob_lo[0] + (i + 0.5) * dx[0];
              const Real y = prob_lo[1] + (j + 0.5) * dx[1];

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

    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        // const auto prob_lo         = geomdata.ProbLo();
        // const auto dx              = geomdata.CellSize();

        // const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        // const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = z_r(i,j,k);

        state(i, j, k, Temp_comp) = 1.;

        state(i,j,k,Temp_comp)=m_solverChoice.T0+7.5_rt*std::exp(z/1000.0_rt);
        if (l_use_salt) {
            state(i,j,k,Salt_comp)=m_solverChoice.S0;
        }

        // Set scalar = 0 everywhere
        state(i, j, k, Scalar_comp) = 0.0;
    });

  // Construct a box that is on x-faces
  const Box& xbx = surroundingNodes(bx,0);
  // Set the x-velocity
  ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
        // Set the x-velocity
        x_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on y-faces
  const Box& ybx = surroundingNodes(bx,1);

  // Set the y-velocity
  ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
        y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const Box& zbx = surroundingNodes(bx,2);

  // Set the z-velocity
  ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });

  Gpu::streamSynchronize();
}

void
init_custom_vmix(const Geometry& /*geom*/, MultiFab& mf_Akv, MultiFab& mf_Akt,
                 MultiFab& /*mf_z_w*/, const SolverChoice& /*m_solverChoice*/)
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
        Akt(i,j,k) = 1.0e-6_rt;
      });
    }
}

void
init_custom_hmix(const Geometry& /*geom*/, MultiFab& mf_visc2_p, MultiFab& mf_visc2_r,
                 MultiFab& mf_diff2, const SolverChoice& /*m_solverChoice*/)
{
    for ( MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& visc2_p_arr = mf_visc2_p.array(mfi);
      Array4<Real> const& visc2_r_arr = mf_visc2_r.array(mfi);
      Array4<Real> const& diff2_arr   = mf_diff2.array(mfi);

      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      int ncomp = 2; // temperature and salt
      Gpu::streamSynchronize();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
          visc2_p_arr(i,j,k) = 0.0;
          visc2_r_arr(i,j,k) = 0.0;

          for (int n = 0; n < ncomp; n++) {
              diff2_arr(i,j,k,n) = 0.0;
          }
      });
    } // mfi
}

void
init_custom_smflux(const Geometry& /*geom*/, const Real /*time*/,
                   MultiFab& mf_sustr, MultiFab& mf_svstr, const SolverChoice& /*m_solverChoice*/)
{
    mf_sustr.setVal(0.0);
    mf_svstr.setVal(0.0);
}

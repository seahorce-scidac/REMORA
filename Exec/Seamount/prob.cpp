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
  pp.query("R0", parms.R0);
  pp.query("S0", parms.S0);
  pp.query("T0", parms.T0);

  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
  pp.query("A_0", parms.A_0);
  pp.query("B_0", parms.B_0);
  pp.query("u_0", parms.u_0);
  pp.query("v_0", parms.v_0);
  pp.query("rad_0", parms.rad_0);
  pp.query("z0", parms.z0);
  pp.query("zRef", parms.zRef);
  pp.query("uRef", parms.uRef);

  pp.query("xc_frac", parms.xc_frac);
  pp.query("yc_frac", parms.yc_frac);
  pp.query("zc_frac", parms.zc_frac);

  pp.query("prob_type", parms.prob_type);
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
    //std::unique_ptr<MultiFab>& mf_z_w = vec_z_w[lev];
    //std::unique_ptr<MultiFab>& mf_h  = vec_hOfTheConfusingName[lev];
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    mf_h.setVal(geom.ProbHi(2));
    Real depth = geom.ProbHi(2);
    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    //HACK HACK manually setting zeta to 0
    mf_Zt_avg1.setVal(0.0);

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

      const auto & geomdata = geom.data();

      int ncomp = 1;
      //NOTE: this will need to be updated for multilevel
      Box subdomain = geom.Domain();

      int nx = subdomain.length(0);
      int ny = subdomain.length(1);
      int nz = subdomain.length(2);

      auto N = nz; // Number of vertical "levels" aka, NZ
      bool NSPeriodic = geomdata.isPeriodic(1);
      bool EWPeriodic = geomdata.isPeriodic(0);

      amrex::Real Xsize = 320000.0_rt;
      amrex::Real Esize = 320000.0_rt;
      amrex::Real depth = 5000.0_rt;
      amrex::Real f0 = 1e-4;
      amrex::Real beta = 1.0_rt;

      if(!m_solverChoice.flat_bathymetry) {
      Gpu::streamSynchronize();
      amrex::ParallelFor(gbx2, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
      {
          const auto prob_lo         = geomdata.ProbLo();
          const auto dx              = geomdata.CellSize();

          const Real x = prob_lo[0] + (i + 0.5) * dx[0];
          const Real y = prob_lo[1] + (j + 0.5) * dx[1];

          Real val1, val2;
          int iFort = i+1;
          int jFort = j+1;
          val1 = (x-0.5_rt*Xsize)/40000.0_rt;
          val2 = (y-0.5_rt*Esize)/40000.0_rt;
          h(i,j,0) = depth - 4500.0_rt * std::exp(-(val1*val1+val2*val2));
      });
      } else {
      Gpu::streamSynchronize();
      amrex::ParallelFor(gbx2, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
      {
          h(i,j,0) = -geomdata.ProbLo(2);
          if (k==0) {
              h(i,j,0,1) = h(i,j,0);
          }
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
        Array4<Real const> const& Hz,
        Array4<Real const> const& h,
        Array4<Real const> const& Zt_avg1,
        GeometryData const& geomdata,
        const SolverChoice& m_solverChoice)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  const Real& rho_sfc   = p_0 / (R_d*parms.T0);
  const Real& thetabar  = parms.T0;
  const Real& dz        = geomdata.CellSize()[2];
  const Real& el        = geomdata.ProbHi()[1];
  const Real& prob_lo_z = geomdata.ProbLo()[2];

    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = z_r(i,j,k);

        state(i, j, k, Temp_comp) = 1.;
        state(i, j, k, Rho_comp) = 1.;

        state(i,j,k,Temp_comp)=parms.T0+7.5_rt*std::exp(z/1000.0_rt);
#ifdef ROMSX_USE_SALINITY
        state(i,j,k,Salt_comp)=parms.S0;
#endif

        // Set scalar = 0 everywhere
        state(i, j, k, RhoScalar_comp) = parms.rho0;
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
init_custom_vmix(const Geometry& geom, MultiFab& mf_Akv, MultiFab& mf_Akt, MultiFab& mf_z_w, const SolverChoice& m_solverChoice)
{
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Array4<Real> const& z_w = (mf_z_w).array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      const auto & geomdata = geom.data();
      int ncomp = 1; Gpu::streamSynchronize();
      amrex::ParallelFor(bx, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
      {
        Akv(i,j,k) = 1.0e-5_rt;
        Akt(i,j,k) = 1.0e-6_rt;
      });
    }
}

void
init_custom_hmix(const Geometry& geom, MultiFab& mf_visc2_p, MultiFab& mf_visc2_r,
    MultiFab& mf_diff2_salt, MultiFab& mf_diff2_temp, const SolverChoice& m_solverChoice)
{
    for ( MFIter mfi((mf_visc2_p), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& visc2_p = (mf_visc2_p).array(mfi);
      Array4<Real> const& visc2_r = (mf_visc2_r).array(mfi);
      Array4<Real> const& diff2_salt = (mf_diff2_salt).array(mfi);
      Array4<Real> const& diff2_temp = (mf_diff2_temp).array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));
      const auto & geomdata = geom.data();
      int ncomp = 1;
      Gpu::streamSynchronize();
      amrex::ParallelFor(bx, ncomp,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
      {
          visc2_p(i,j,k) = 0.0;
          visc2_r(i,j,k) = 0.0;

          diff2_salt(i,j,k) = 0.0;
          diff2_temp(i,j,k) = 0.0;
      });
    }
}

void
init_custom_smflux(const Geometry& geom, const Real time, MultiFab& mf_sustr, MultiFab& mf_svstr,
    const SolverChoice& m_solverChoice)
{
    mf_sustr.setVal(0.0);
    mf_svstr.setVal(0.0);
}

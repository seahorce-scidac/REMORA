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

void
init_custom_prob(
        const Box& bx,
        Array4<Real      > const& state,
        Array4<Real      > const& x_vel,
        Array4<Real      > const& y_vel,
        Array4<Real      > const& z_vel,
        Array4<Real      > const& r_hse,
        Array4<Real      > const& p_hse,
        Array4<Real const> const& z_nd,
        Array4<Real const> const& z_cc,
        Array4<Real const> const& z_w,
        Array4<Real const> const& z_r,
        Array4<Real const> const& Hz,
        Array4<Real const> const& h,
        Array4<Real const> const& Zt_avg1,
        GeometryData const& geomdata)
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

        state(i,j,k,Temp_comp)=parms.T0+8.0*std::exp(z/50.0_rt);
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
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = -z_r(i,j,k);

	// Set the x-velocity
	x_vel(i, j, k) = parms.u_0 + parms.uRef *
	                 std::log((z + parms.z0)/parms.z0)/
                         std::log((parms.zRef +parms.z0)/parms.z0);
	//x_vel(i, j, k) = 1.0 / x;
  });

  // Construct a box that is on y-faces
  const Box& ybx = surroundingNodes(bx,1);

  // Set the y-velocity
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        y_vel(i, j, k) = parms.v_0;
        //y_vel(i, j, k) = 1.0 / y;
  });

  // Construct a box that is on z-faces
  const Box& zbx = surroundingNodes(bx,2);

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });

  Gpu::streamSynchronize();
}

void
init_custom_terrain(const Geometry& /*geom*/, MultiFab& z_phys_nd,
                    const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

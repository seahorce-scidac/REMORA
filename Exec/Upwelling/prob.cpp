#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"

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
  /*
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
  */
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
        GeometryData const& geomdata)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  const Real& rho_sfc   = p_0 / (R_d*parms.T0);
  const Real& thetabar  = parms.T0;
  const Real& dz        = geomdata.CellSize()[2];
  const Real& el        = geomdata.ProbHi()[1];
  const Real& prob_lo_z = geomdata.ProbLo()[2];

    amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        state(i, j, k, RhoTheta_comp) = 1.;
        state(i, j, k, Rho_comp) = 1.;

#if 0
	const Real val1=(44.69_rt/39.382_rt)*(44.69_rt/39.382_rt);
	const Real val2=val1*(parms.rho0*100.0_rt/parms.g)*(5.0E-5_rt/((42.689_rt/44.69_rt)*(42.689_rt/44.69_rt)));
        Real val3=parms.T0+val2*std::exp(z_r(i,j,k)/100.0_rt)*
	  (10.0_rt-0.4_rt*tanh(z_r(i,j,k)/100.0_rt));
	Real val4=yr(i,j)/el;
        state(i,j,k,RhoTheta_comp)=val3-3.0_rt*val4; // This may be missing rho effects
#ifdef ROMSX_USE_SALINITY
        state(i,j,k,Salt_comp)=34.5_rt-0.001_rt*z_r(i,j,k)-val4;
#endif
#else
	const Real val1=(44.69_rt/39.382_rt)*(44.69_rt/39.382_rt);
	const Real val2=val1*(parms.rho0*100.0_rt/parms.g)*(5.0E-5_rt/((42.689_rt/44.69_rt)*(42.689_rt/44.69_rt)));
        Real val3=parms.T0+val2*std::exp(z/100.0_rt)*
	  (10.0_rt-0.4_rt*tanh(z/100.0_rt));
	Real val4=y/el;
        state(i,j,k,RhoTheta_comp)=val3-3.0_rt*val4; // This may be missing rho effects
#ifdef ROMSX_USE_SALINITY
        state(i,j,k,Salt_comp)=34.5_rt-0.001_rt*z-val4;
#endif
#endif

        // Set scalar = 0 everywhere
        state(i, j, k, RhoScalar_comp) = parms.rho0;
    });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
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
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

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

  pp.query("u_0", parms.u_0);
  pp.query("v_0", parms.v_0);
  pp.query("z0", parms.z0);
  pp.query("zRef", parms.zRef);
  pp.query("velRef", parms.velRef);

}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void
init_custom_bathymetry (const Geometry& geom,
                        MultiFab& mf_h,
                        const SolverChoice& m_solverChoice)
{
    auto geomdata = geom.data();
    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW,NGROW,0));

      Box gbx2D = gbx2;
      gbx2D.makeSlab(2,0);

      Gpu::streamSynchronize();

      if (NSPeriodic) {

          ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              int iFort = i+1; // (+1 is to match the Fortran indexing in ROMS)

              Real val1 = (iFort <= Lm/2.0) ? iFort : Lm+1-iFort;

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5+66.526*std::tanh((val1-10.0)/7.0)));
          });

      } else if (EWPeriodic) {

          ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              int jFort = j+1; // (+1 is to match the Fortran indexing in ROMS)

              Real val1 = (jFort<=Mm/2.0) ? jFort : Mm+1-jFort;

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5+66.526*std::tanh((val1-10.0)/7.0)));
          });
      }
    } // mfi
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

    bool EWPeriodic = geomdata.isPeriodic(0);
    bool NSPeriodic = geomdata.isPeriodic(1);

    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo         = geomdata.ProbLo();
        const auto prob_hi         = geomdata.ProbHi();
        const auto dx              = geomdata.CellSize();

        const Real z = z_r(i,j,k);

        state(i, j, k, Temp_comp) = 1.;

        state(i,j,k,Temp_comp)=m_solverChoice.T0+8.0*std::exp(z/50.0_rt);
        if (l_use_salt) {
            state(i,j,k,Salt_comp)=m_solverChoice.S0;
        }

        // Set scalar = 0 everywhere
        const Real xcent = 0.5*(prob_lo[0] + prob_hi[0]);
        const Real ycent = 0.5*(prob_lo[1] + prob_hi[1]);

        const Real x  = prob_lo[0] + (i + 0.5) * dx[0] - xcent;
        const Real y  = prob_lo[1] + (j + 0.5) * dx[1] - ycent;
        const Real r2 = x*x + y*y;
        const Real rad = 0.1 * (prob_hi[0]-prob_lo[0]);
        const Real radsq = rad*rad;

        state(i, j, k, Scalar_comp) = 0.0;
    });

  const Box& xbx = surroundingNodes(bx,0);
  const Box& ybx = surroundingNodes(bx,1);

   // Below we set the velocity to be non-zero only in the periodic direction, and
   //      zero it in the non-periodic direction

  if (NSPeriodic) {
      ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          const Real z = -z_r(i,j,k);

          x_vel(i, j, k) = parms.u_0 + parms.velRef *
                           std::log((z + parms.z0)/parms.z0)/
                           std::log((parms.zRef +parms.z0)/parms.z0);
      });
  } else {
      ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          x_vel(i, j, k) = 0.0;
      });
  }


  if (EWPeriodic) {
      ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          const Real z = -z_r(i,j,k);

          x_vel(i, j, k) = parms.v_0 + parms.velRef *
                           std::log((z + parms.z0)/parms.z0)/
                           std::log((parms.zRef +parms.z0)/parms.z0);
      });
  } else {
      ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          y_vel(i, j, k) = 0.0;
      });
  }

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
        Akv(i,j,k) = 2.0e-03+8.0e-03*std::exp(z_w(i,j,k)/150.0);

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
        visc2_p(i,j,k) = 5.0;
        visc2_r(i,j,k) = 5.0;

        for (int n = 0; n < ncomp; n++) {
            diff2(i,j,k,n) = 0.0;
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
    Real pi = 3.14159265359;
    Real tdays=time/(24.0*60.0*60.0);
    Real dstart=0.0;
    Real windamp;
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

    // Flow in x-direction (EW):
    if (NSPeriodic) {
        mf_sustr.setVal(0.0);
    }
    else if (EWPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/m_solverChoice.rho0;
        else
            windamp=-0.1/m_solverChoice.rho0;
        mf_sustr.setVal(windamp);
    }

    // Flow in y-direction (NS):
    if (NSPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/m_solverChoice.rho0;
        else
            windamp=-0.1/m_solverChoice.rho0;
        mf_svstr.setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr.setVal(0.0);
    }
}

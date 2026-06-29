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
        int lev, const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_h)
{
    mf_h.setVal(geom.ProbHi(2));
    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

      const auto & geomdata = geom.data();

      bool NSPeriodic = geomdata.isPeriodic(1);
      bool EWPeriodic = geomdata.isPeriodic(0);

      Gpu::streamSynchronize();

      if (!m_solverChoice.flat_bathymetry) {
          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              Real val1;
              int iFort = i+1;
              int jFort = j+1;
              if(NSPeriodic) {
                  if (iFort<=Lm/2.0)
                      val1=iFort;
                  else
                      val1=Lm+1-iFort;
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1-10.0_rt)/7.0_rt)));

              } else if(EWPeriodic) {

                  if (jFort<=Mm/2.0)
                      val1=jFort;
                  else
                      val1=Mm+1-jFort;
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1-10.0_rt)/7.0_rt)));
              }
          });
      } else {

          ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
              h(i,j,0) = -geomdata.ProbLo(2);
              if (k==0) {
                  h(i,j,0,1) = h(i,j,0);
              }
          });
      }
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
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel)
{
    ParmParse pp("remora.prob");
    Real u_0  = 0.0; pp.query("u_0", u_0);
    Real v_0  = 0.0; pp.query("v_0", v_0);
    Real z0   = 0.1;  pp.query("z0", z0);     // Surface Roughness
    Real zRef = 80.0; pp.query("zRef", zRef); // Reference Height
    Real uRef = 0.0;  pp.query("uRef", uRef); // Reference Wind Speed

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    bool l_use_salt = m_solverChoice.use_salt;

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

        Real S0 = m_solverChoice.S0;
        Real T0 = m_solverChoice.T0;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Geometry (note we must include these here to get the data on device)
            // const auto prob_lo         = geomdata.ProbLo();
            // const auto dx              = geomdata.CellSize();

            // const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            // const Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const Real z = z_r(i,j,k);

            state(i, j, k, Temp_comp) = 1.;

            state(i,j,k,Temp_comp)=T0+8.0_rt*std::exp(z/50.0_rt);
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=S0;
            }

            // Set tracer = 0 everywhere
            state(i, j, k, Tracer_comp) = 0.0_rt;
        });

        // Construct a box that is on x-faces
        const Box& xbx = surroundingNodes(bx,0);
        // Set the x-velocity
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real z = -z_r(i,j,k);

            // Set the x-velocity
            x_vel(i, j, k) = u_0 + uRef * std::log((z + z0)/z0) / std::log((zRef +z0)/z0);
        });

        // Construct a box that is on y-faces
        const Box& ybx = surroundingNodes(bx,1);

        // Set the y-velocity
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
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Array4<const Real> const& z_w = remora.vec_z_w[lev]->const_array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      Gpu::streamSynchronize();

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
        visc2_p(i,j,0) = 5.0_rt;
        visc2_r(i,j,0) = 5.0_rt;

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
    auto geomdata = geom.data();
    bool NSPeriodic = geomdata.isPeriodic(1);
    bool EWPeriodic = geomdata.isPeriodic(0);
    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359;
    Real tdays=remora.get_t_old(lev)/(24.0*60.0*60.0);
    //this is a hack because time is off by dt. this needs to be fixed for non-fixed dt
    Real rho0=m_solverChoice.rho0;
    Real windamp;
    Real dstart = 0.0;
    //It's possible these should be set to be nonzero only at the boundaries they affect
    if(NSPeriodic) {
        mf_sustr.setVal(0.0);
    }
    else if(EWPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/rho0;
        else
            windamp=-0.1/rho0;
        mf_sustr.setVal(windamp);
    }
    if(NSPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/rho0;
        else
            windamp=-0.1/rho0;
        mf_svstr.setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr.setVal(0.0);
    }
}

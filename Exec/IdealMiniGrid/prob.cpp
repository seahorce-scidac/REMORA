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
init_custom_bathymetry (const Geometry& /*geom*/,
                        MultiFab& /*mf_h*/,
                        const SolverChoice& /*m_solverChoice*/)
{
    Abort("Shouldn't be in init_custom_bathymetry!");
}

/**
 * \brief Initializes coriolis forcing
 */
void
init_custom_coriolis (const Geometry& /*geom*/,
                      MultiFab& /*mf_fcor*/,
                      const SolverChoice& /*m_solverChoice*/) {}

/**
 * \brief Initializes custom sea surface height
 */
void
init_custom_zeta (const Geometry& geom,
                      MultiFab& mf_zeta,
                      const SolverChoice& m_solverChoice) {}

void
init_custom_prob(
        const Box& /*bx*/,
        Array4<Real      > const& /*state*/,
        Array4<Real      > const& /*x_vel*/,
        Array4<Real      > const& /*y_vel*/,
        Array4<Real      > const& /*z_vel*/,
        Array4<Real const> const& /*z_w*/,
        Array4<Real const> const& /*z_r*/,
        Array4<Real const> const& /*Hz*/,
        Array4<Real const> const& /*h*/,
        Array4<Real const> const& /*Zt_avg1*/,
        GeometryData const& /*geomdata*/,
        const SolverChoice& /*m_solverChoice*/)
{
    Abort("Shouldn't be in init_custom_prob!");
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
      Array4<Real> const& visc2_p_arr = mf_visc2_p.array(mfi);
      Array4<Real> const& visc2_r_arr = mf_visc2_r.array(mfi);
      Array4<Real> const& diff2_arr   = mf_diff2.array(mfi);

      Box bx = mfi.tilebox();
      bx.grow(IntVect(NGROW,NGROW,0));

      int ncomp = mf_diff2.nComp(); // temperature and salt and scalar
      Gpu::streamSynchronize();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
          visc2_p_arr(i,j,k) = 5.0e0;
          visc2_r_arr(i,j,k) = 5.0e0;

          for (int n = 0; n < ncomp; n++) {
              diff2_arr(i,j,k,n) = 5.0e3;
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

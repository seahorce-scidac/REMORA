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
    ParmParse pp("remora.prob");

    pp.query("traditional", parms.traditional);
}

/**
 * \brief Initializes bathymetry h and surface height Zeta
 */
void Problem::init_analytic_bathymetry (
        int lev, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        amrex::MultiFab& mf_h)
{
    if (parms.traditional) {
        mf_h.setVal(10.0_rt);
    } else {
        mf_h.setVal(0.0_rt);
        for (MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

            Array4<      Real> const& h = mf_h.array(mfi);
            Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);
            Array4<const Real> const& y_r  = remora.vec_yr[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
            {
                h(i,j,0) = 10.0_rt + (x_r(i,j,0) - 4200._rt) * (x_r(i,j,0) - 4200._rt) * 5.0_rt / (4200.0_rt * 4200.0_rt) -
                            (y_r(i,j,0) - 375._rt) * (y_r(i,j,0) - 375._rt) * 5.0_rt / (375._rt * 375._rt);
            });
        }
    }
}

void Problem::init_analytic_zeta (
        int lev, const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_zeta)
{
    mf_zeta.setVal(0.0_rt);

    for (MFIter mfi(mf_zeta, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.growntilebox(IntVect(1,1,0));

        Array4<      Real> const& zeta = mf_zeta.array(mfi);
        Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);

        ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE(int i, int j, int , int n) noexcept
        {
            if (x_r(i,j,0) < 1100.0_rt) {
                zeta(i,j,0,n) = -0.00125_rt * (x_r(i,j,0) - 100.0_rt) + 1.25_rt;
            }
        });
    }

}

void Problem::init_analytic_prob(
        int lev,
        const amrex::Geometry& geom,
        SolverChoice const& m_solverChoice,
        REMORA const& remora,
        amrex::MultiFab& mf_cons,
        amrex::MultiFab& mf_xvel,
        amrex::MultiFab& mf_yvel,
        amrex::MultiFab& mf_zvel)
{
    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];
    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);
        Array4<      Real> const& z_vel = mf_zvel.array(mfi);

        Array4<const Real> const& x_r = remora.vec_xr[lev]->const_array(mfi);
        if (parms.traditional) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (l_use_salt) {
                    state(i, j, k, Salt_comp) = 35.0_rt;
                }

                state(i,j,k,Temp_comp)= (x_r(i,j,0) < 1000._rt) ? 5.0_rt : 10.0_rt;

                // Set tracer = 0 everywhere
                state(i, j, k, Tracer_comp) = 0.0_rt;
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                [[maybe_unused]] const Real x = x_r(i,j,0);

                if (l_use_salt) {
                    state(i, j, k, Salt_comp) = 10.0_rt - 10.0_rt * (x_r(i,j,0) - 8400.0_rt) / 8400.0_rt;
                }

                state(i,j,k,Temp_comp)= 10.0_rt;

                // Set tracer = 0 everywhere
                state(i, j, k, Tracer_comp) = 0.0_rt;
            });
        }

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

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel(i, j, k) = 0.0_rt;
        });
    }
    Gpu::streamSynchronize();
}

void Problem::init_analytic_vmix(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& m_solverChoice,
        REMORA const& /*remora*/,
        MultiFab& mf_Akv, MultiFab& mf_Akt)
{
    mf_Akv.setVal(m_solverChoice.Akv_bak);
    mf_Akt.setVal(m_solverChoice.Akt_bak);
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
    mf_visc2_p.setVal(0.0_rt);
    mf_visc2_r.setVal(0.0_rt);
    mf_diff2.setVal(0.0_rt);
}

void Problem::init_analytic_smflux(
        int /*lev*/,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& /*remora*/,
        MultiFab& mf_sustr, MultiFab& mf_svstr)
{
    mf_svstr.setVal(0.0_rt);
    mf_sustr.setVal(0.0_rt);
}

void Problem::init_analytic_masks(
        int lev,
        const amrex::Geometry& /*geom*/,
        SolverChoice const& /*m_solverChoice*/,
        REMORA const& remora,
        MultiFab& mf_mskr)
{
    for (MFIter mfi(mf_mskr, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        Array4<      Real> const& mskr = mf_mskr.array(mfi);
        Array4<const Real> const& x_r  = remora.vec_xr[lev]->const_array(mfi);
        Array4<const Real> const& y_r  = remora.vec_yr[lev]->const_array(mfi);

        // mskr is 1 everywhere already
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            if ((x_r(i,j,0) > 2800.0_rt && x_r(i,j,0) < 5600.0_rt) &&
                (y_r(i,j,0) < 250.0_rt || y_r(i,j,0) > 500.0_rt)) {
                mskr(i,j,0) = 0.0_rt;
            }
        });
    }
}


#ifdef REMORA_USE_PARTICLES

#include <AMReX_ParticleInterpolators.H>
#include <REMORA_Constants.H>
#include <REMORA_PC.H>

using namespace amrex;

/**
 * @param a_mf     Multifab to calculate over
 * @param a_lev    level to calculate on
 * @param a_comp   component to calculate
 */
void REMORAPC::massDensity ( MultiFab&  a_mf,
                          const int& a_lev,
                          const int& a_comp ) const
{
    BL_PROFILE("REMORAPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(zero);

    ParticleToMesh( *this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (  const REMORAPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const REMORAPC::ParticleType&, int)
                {
                    auto mass = ptd.m_rdata[REMORAParticlesRealIdxSoA::mass][i];
                    return static_cast<Real>(mass*inv_cell_volume);
                });
        });

    return;
}

/**
 * Write the particles to a plotfile with x and y in physical coordinates. The particles
 * are advected in the uniform coordinates of the AMReX geometry; for output, x and y are
 * replaced by the grid's x_rho and y_rho interpolated bilinearly to each particle, as ROMS
 * does for floats. z is already the physical height. Only a copy of the particles is
 * converted, so writing a plotfile leaves the particles, and so the run, unchanged.
 *
 * The interpolation stencil is kept inside the domain, extrapolating linearly over the
 * outer half cells, so it never reads x_rho or y_rho outside the domain, where a periodic
 * fill would hold values from the opposite side.
 *
 * @param a_fname   plotfile directory
 * @param a_name    particle species name
 * @param a_xr      x at rho points, on each level
 * @param a_yr      y at rho points, on each level
 */
void REMORAPC::WritePlotFile ( const std::string& a_fname, const std::string& a_name,
                               const Vector<std::unique_ptr<MultiFab>>& a_xr,
                               const Vector<std::unique_ptr<MultiFab>>& a_yr ) const
{
    BL_PROFILE("REMORAPC::WritePlotFile()");

    auto pc = make_alike();
    pc.copyParticles(*this, true);

    using PCType = decltype(pc);
    for (int lev = 0; lev <= pc.finestLevel(); ++lev)
    {
        const Geometry& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const Box& dom = geom.Domain();
        const int ilo = dom.smallEnd(0);
        const int ihi = dom.bigEnd(0);
        const int jlo = dom.smallEnd(1);
        const int jhi = dom.bigEnd(1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (typename PCType::ParIterType pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& aos = pti.GetArrayOfStructs();
            const int np = aos.numParticles();
            auto* p_pbox = aos().data();
            const auto xr = (*a_xr[lev])[pti.index()].const_array();
            const auto yr = (*a_yr[lev])[pti.index()].const_array();

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int n)
            {
                auto& p = p_pbox[n];
                if (p.id() <= 0) { return; }

                // Cell-centered index coordinates, with the lower stencil index clamped so
                // that both stencil cells are inside the domain
                const Real lx = static_cast<Real>((p.pos(0)-plo[0])*dxi[0]) - myhalf;
                const Real ly = static_cast<Real>((p.pos(1)-plo[1])*dxi[1]) - myhalf;
                const int i0 = amrex::max(ilo, amrex::min(static_cast<int>(amrex::Math::floor(lx)), ihi-1));
                const int j0 = amrex::max(jlo, amrex::min(static_cast<int>(amrex::Math::floor(ly)), jhi-1));
                const int i1 = amrex::min(i0+1, ihi);
                const int j1 = amrex::min(j0+1, jhi);
                const Real sx = (ihi > ilo) ? lx - static_cast<Real>(i0) : zero;
                const Real sy = (jhi > jlo) ? ly - static_cast<Real>(j0) : zero;

                const Real x = xr(i0,j0,0) * (one-sx) * (one-sy) + xr(i1,j0,0) * sx * (one-sy)
                             + xr(i0,j1,0) * (one-sx) *      sy  + xr(i1,j1,0) * sx *      sy;
                const Real y = yr(i0,j0,0) * (one-sx) * (one-sy) + yr(i1,j0,0) * sx * (one-sy)
                             + yr(i0,j1,0) * (one-sx) *      sy  + yr(i1,j1,0) * sx *      sy;
                p.pos(0) = static_cast<ParticleReal>(x);
                p.pos(1) = static_cast<ParticleReal>(y);
            });
        }
    }

    pc.Checkpoint(a_fname, a_name, true, varNames());
}

#endif

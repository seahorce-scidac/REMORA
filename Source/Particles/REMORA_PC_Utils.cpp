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
 * Write the particles to a plotfile with x and y in physical coordinates. The particles move
 * in the uniform coordinates of the AMReX geometry; for output, x and y are replaced by the
 * physical coordinates of the cell corners (psi points) interpolated bilinearly to each
 * particle. These are the corner positions the plotfile's nodal displacement gives a viewer,
 * so each particle is drawn inside the cell that holds it. z is already the physical height.
 * Only a copy of the particles is converted, so writing a plotfile does not change the run.
 *
 * @param a_fname   plotfile directory
 * @param a_name    particle species name
 * @param a_xy      per level, physical x and y of the psi points (REMORA::plot_node_coords)
 */
void REMORAPC::WritePlotFile ( const std::string& a_fname, const std::string& a_name,
                               const Vector<MultiFab>& a_xy ) const
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (typename PCType::ParIterType pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& aos = pti.GetArrayOfStructs();
            const int np = aos.numParticles();
            auto* p_pbox = aos().data();
            const auto xy = a_xy[lev].const_array(pti.index());
            const Box& vbx = pti.validbox();
            const int ilo = vbx.smallEnd(0);
            const int ihi = vbx.bigEnd(0);
            const int jlo = vbx.smallEnd(1);
            const int jhi = vbx.bigEnd(1);

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int n)
            {
                auto& p = p_pbox[n];
                if (p.id() <= 0) { return; }

                // Nodal index coordinates, kept in ParticleReal: rounding them to a
                // single-precision Real first can carry a particle at the high edge of its grid
                // into the next cell. The cell is clamped to the particle's grid, so both corner
                // nodes are always inside the grid's nodal box (a_xy has no ghost cells).
                const ParticleReal lx = (p.pos(0)-plo[0])*dxi[0];
                const ParticleReal ly = (p.pos(1)-plo[1])*dxi[1];
                const int i = amrex::max(ilo, amrex::min(static_cast<int>(amrex::Math::floor(lx)), ihi));
                const int j = amrex::max(jlo, amrex::min(static_cast<int>(amrex::Math::floor(ly)), jhi));
                const auto sx = static_cast<Real>(lx - static_cast<ParticleReal>(i));
                const auto sy = static_cast<Real>(ly - static_cast<ParticleReal>(j));

                for (int d = 0; d < 2; ++d) {
                    p.pos(d) = static_cast<ParticleReal>(
                        xy(i  ,j  ,0,d) * (one-sx) * (one-sy) + xy(i+1,j  ,0,d) * sx * (one-sy) +
                        xy(i  ,j+1,0,d) * (one-sx) *      sy  + xy(i+1,j+1,0,d) * sx *      sy);
                }
            });
        }
    }

    pc.Checkpoint(a_fname, a_name, true, varNames());
}

#endif

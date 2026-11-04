#ifdef REMORA_USE_PARTICLES

#include <REMORA_PC.H>
#include <REMORA_IndexDefines.H>
#include <REMORA_Constants.H>

using namespace amrex;

/**
 * @param a_lev         level at which to evolve particles
 * @param a_dt_lev      time step
 * @param a_flow_vel    flow velocities
 * @param a_z_phys_nd   z nodal positions of grid
 * @param a_pm          1/dx, the inverse physical cell width in x
 * @param a_pn          1/dy, the inverse physical cell width in y
 */
void REMORAPC::EvolveParticles ( int                                        a_lev,
                                 Real                                       a_dt_lev,
                                 Vector<MultiFab const*>&                   a_flow_vel,
                                 const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd,
                                 MultiFab const&                            a_pm,
                                 MultiFab const&                            a_pn )
{
    BL_PROFILE("REMORAPCPC::EvolveParticles()");

    if (m_advect_w_flow) {
        AdvectWithFlow( a_lev, a_dt_lev, a_flow_vel, a_z_phys_nd[a_lev], a_pm, a_pn );
    }

    Redistribute();
    return;
}

/**
 * Midpoint-method advection in terrain-following coordinates. Horizontally the particle
 * moves with u and v at a fixed fractional w-level index (see REMORA_PC.H), so it follows
 * the s-levels over sloping bathymetry; vertically it moves across the levels at
 * omega / Hz. The index is reflected off the bed and the free surface and then mapped
 * back to a height, so a particle can never leave the water column.
 *
 * Particle positions are in the uniform coordinates of the AMReX geometry, whose cell size
 * is (prob_hi - prob_lo) / n_cell, while u and v are physical velocities along the grid
 * directions. They are converted with the local metric, dx/dt = u * pm * dx_geom and
 * dy/dt = v * pn * dy_geom, so particles cover the right physical distance on stretched or
 * curvilinear grids; on a uniform grid the factors are 1.
 *
 * The index is recovered from the stored height on the current grid, so a particle stays
 * at a fixed height, not a fixed index, while the free surface moves between steps. Both
 * stages of the midpoint method run in one kernel; nothing between them needs the particles
 * redistributed, since the midpoint stays within the ghost cells of the particle's grid.
 *
 * @param a_lev         level at which to evolve particles
 * @param a_dt          time step
 * @param a_umac        u, v and omega, the velocity across the w-levels (z-nodal, m/s)
 * @param a_z_height    z nodal positions of grid
 * @param a_pm          1/dx, the inverse physical cell width in x
 * @param a_pn          1/dy, the inverse physical cell width in y
 */
void REMORAPC::AdvectWithFlow ( int                                 a_lev,
                                Real                                a_dt,
                                Vector<MultiFab const*>&            a_umac,
                                const std::unique_ptr<MultiFab>&    a_z_height,
                                MultiFab const&                     a_pm,
                                MultiFab const&                     a_pn )
{
    BL_PROFILE("REMORAPCPC::AdvectWithUmac()");
    AMREX_ASSERT(OK(a_lev, a_lev, a_umac[0]->nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

    AMREX_ALWAYS_ASSERT(a_z_height != nullptr);
    AMREX_D_TERM(AMREX_ASSERT(a_umac[0]->nGrow() >= 1);,
                 AMREX_ASSERT(a_umac[1]->nGrow() >= 1);,
                 AMREX_ASSERT(a_umac[2]->nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!a_umac[0]->contains_nan());,
                 AMREX_ASSERT(!a_umac[1]->contains_nan());,
                 AMREX_ASSERT(!a_umac[2]->contains_nan()););

    const auto      strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const int nz = geom.Domain().length(2);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        int grid    = pti.index();
        auto& ptile = ParticlesAt(a_lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        auto& soa  = ptile.GetStructOfArrays();
        const int n = aos.numParticles();
        auto *p_pbox = aos().data();

        // Set at the end of the step: the physical velocities u and v at the midpoint (m/s,
        // not the rate of change of the geometry coordinates) and vz, the rate of change of
        // the particle's height over the step
        auto* vx_ptr = soa.GetRealData(REMORAParticlesRealIdxSoA::vx).data();
        auto* vy_ptr = soa.GetRealData(REMORAParticlesRealIdxSoA::vy).data();
        auto* vz_ptr = soa.GetRealData(REMORAParticlesRealIdxSoA::vz).data();

        const auto uarr = (*a_umac[0])[grid].const_array();
        const auto varr = (*a_umac[1])[grid].const_array();
        const auto warr = (*a_umac[2])[grid].const_array();
        const auto zheight = (*a_z_height)[grid].const_array();
        const auto pmarr = a_pm[grid].const_array();
        const auto pnarr = a_pn[grid].const_array();

        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            // u, v and the rates of change of x, y and s at (a_x, a_y) and level index a_loc
            const auto rates = [&] (ParticleReal a_x, ParticleReal a_y, LevelIndex const& a_loc,
                                    Real& a_u, Real& a_v, Real& a_xdot, Real& a_ydot, Real& a_sdot)
            {
                a_u = interpolate_to_level_index(a_x, a_y, a_loc.s, nz, IntVect(1,0,0),
                                                 plo, dxi, uarr);
                a_v = interpolate_to_level_index(a_x, a_y, a_loc.s, nz, IntVect(0,1,0),
                                                 plo, dxi, varr);
                const Real w = interpolate_to_level_index(a_x, a_y, a_loc.s, nz, IntVect(0,0,1),
                                                          plo, dxi, warr);
                // A cell of zero thickness (a dry or degenerate column) does not move s
                a_sdot = (a_loc.dz > zero) ? w / a_loc.dz : zero;
                // Physical velocity to the rate of change of the geometry coordinate:
                // pm * dx_geom is the number of geometry metres per physical metre
                const HorizontalStencil cc = horizontal_stencil(a_x, a_y, false, false, plo, dxi);
                a_xdot = a_u * interpolate_horizontal(cc, 0, pmarr) / dxi[0];
                a_ydot = a_v * interpolate_horizontal(cc, 0, pnarr) / dxi[1];
            };

            // omega vanishes at the bed and the surface, so only the time
            // discretization can carry a particle past them: reflect it back
            const auto reflect = [=] (Real a_s)
            {
                if (a_s < zero) { a_s = -a_s; }
                if (a_s > static_cast<Real>(nz)) { a_s = Real(2*nz) - a_s; }
                return amrex::max(Real(0.0), amrex::min(static_cast<Real>(nz), a_s));
            };

            Real u, v, xdot, ydot, sdot;

            // Half step to the midpoint, from the start of the step
            const ParticleReal x0 = p.pos(0);
            const ParticleReal y0 = p.pos(1);
            const LevelIndex loc0 = height_to_level_index(x0, y0, p.pos(2),
                                                          p.idata(REMORAParticlesIntIdxAoS::k),
                                                          nz, plo, dxi, zheight);
            rates(x0, y0, loc0, u, v, xdot, ydot, sdot);

            const ParticleReal xh = x0 + static_cast<ParticleReal>(myhalf*a_dt*xdot);
            const ParticleReal yh = y0 + static_cast<ParticleReal>(myhalf*a_dt*ydot);
            const Real sh = reflect(loc0.s + myhalf*a_dt*sdot);

            // Full step with the rates at the midpoint, whose level index is recovered from
            // its height as at the start of a step
            const auto zh = static_cast<ParticleReal>(level_index_to_height(xh, yh, sh, nz,
                                                                            plo, dxi, zheight));
            const LevelIndex loch = height_to_level_index(xh, yh, zh, level_index_cell(sh, nz),
                                                          nz, plo, dxi, zheight);
            rates(xh, yh, loch, u, v, xdot, ydot, sdot);

            p.pos(0) = x0 + static_cast<ParticleReal>(a_dt*xdot);
            p.pos(1) = y0 + static_cast<ParticleReal>(a_dt*ydot);
            const Real s1 = reflect(loc0.s + a_dt*sdot);
            p.pos(2) = static_cast<ParticleReal>(level_index_to_height(p.pos(0), p.pos(1), s1, nz,
                                                                       plo, dxi, zheight));
            p.idata(REMORAParticlesIntIdxAoS::k) = level_index_cell(s1, nz);

            vx_ptr[i] = static_cast<ParticleReal>(u);
            vy_ptr[i] = static_cast<ParticleReal>(v);
            vz_ptr[i] = (p.pos(2) - static_cast<ParticleReal>(loc0.z)) / a_dt;
        });
    }

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

#ifdef AMREX_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                    ParallelContext::CommunicatorSub());

                Print() << "REMORAPC::AdvectWithFlow() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

#endif

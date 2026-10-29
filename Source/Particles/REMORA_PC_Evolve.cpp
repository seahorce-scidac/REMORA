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
 * at a fixed height, not a fixed index, while the free surface moves between steps.
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

    for (int ipass = 0; ipass < 2; ipass++)
    {
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

            // Between the passes these hold the starting x, y and level index; afterwards
            // the physical velocities u and v (m/s, not the rate of change of the geometry
            // coordinates) and vz, the rate of change of the particle's height
            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(REMORAParticlesRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(REMORAParticlesRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(REMORAParticlesRealIdxSoA::vz).data();

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

                const ParticleReal x = p.pos(0);
                const ParticleReal y = p.pos(1);
                int kfound;
                Real dz_found;
                const Real s = height_to_level_index(x, y, p.pos(2),
                                                     p.idata(REMORAParticlesIntIdxAoS::k),
                                                     nz, plo, dxi, zheight, &kfound, &dz_found);

                const Real u = interpolate_to_level_index(x, y, s, nz, IntVect(1,0,0),
                                                          plo, dxi, uarr);
                const Real v = interpolate_to_level_index(x, y, s, nz, IntVect(0,1,0),
                                                          plo, dxi, varr);
                // The search's cell is the one holding s except when k + f rounds, or f = 1
                // puts s on the next level; only then are the level heights recomputed
                const int ks = level_index_cell(s, nz);
                const Real dz = (ks == kfound) ? dz_found
                              : level_height(x, y, ks+1, plo, dxi, zheight)
                              - level_height(x, y, ks  , plo, dxi, zheight);
                const Real sdot = interpolate_to_level_index(x, y, s, nz, IntVect(0,0,1),
                                                             plo, dxi, warr) / dz;
                // Physical velocity to the rate of change of the geometry coordinate:
                // pm * dx_geom is the number of geometry metres per physical metre
                const Real xdot = u * interpolate_horizontal(x, y, 0, false, plo, dxi, pmarr) / dxi[0];
                const Real ydot = v * interpolate_horizontal(x, y, 0, false, plo, dxi, pnarr) / dxi[1];

                Real s_new;
                ParticleReal z_old = zero;
                if (ipass == 0) {
                    v_ptr[0][i] = x;
                    v_ptr[1][i] = y;
                    v_ptr[2][i] = static_cast<ParticleReal>(s);
                    p.pos(0) += static_cast<ParticleReal>(myhalf*a_dt*xdot);
                    p.pos(1) += static_cast<ParticleReal>(myhalf*a_dt*ydot);
                    s_new = s + myhalf*a_dt*sdot;
                } else {
                    const auto s_old = static_cast<Real>(v_ptr[2][i]);
                    z_old = level_index_to_height(v_ptr[0][i], v_ptr[1][i], s_old, nz,
                                                  plo, dxi, zheight);
                    p.pos(0) = v_ptr[0][i] + static_cast<ParticleReal>(a_dt*xdot);
                    p.pos(1) = v_ptr[1][i] + static_cast<ParticleReal>(a_dt*ydot);
                    s_new = s_old + a_dt*sdot;
                }

                // omega vanishes at the bed and the surface, so only the time
                // discretization can carry a particle past them: reflect it back
                if (s_new < zero) { s_new = -s_new; }
                if (s_new > static_cast<Real>(nz)) { s_new = Real(2*nz) - s_new; }
                s_new = amrex::max(zero, amrex::min(static_cast<Real>(nz), s_new));

                p.pos(2) = level_index_to_height(p.pos(0), p.pos(1), s_new, nz,
                                                 plo, dxi, zheight);
                p.idata(REMORAParticlesIntIdxAoS::k) = level_index_cell(s_new, nz);

                if (ipass == 1) {
                    v_ptr[0][i] = u;
                    v_ptr[1][i] = v;
                    v_ptr[2][i] = (p.pos(2) - z_old) / a_dt;
                }
            });
        }
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

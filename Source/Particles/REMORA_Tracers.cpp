#include <string>
#include <REMORA.H>
#include <REMORA_PC.H>

#ifdef REMORA_USE_PARTICLES

using namespace amrex;

/*! Read tracer and hydro particles parameters */
void REMORA::readTracersParams ()
{
    ParmParse pp(pp_prefix);

    m_use_tracer_particles = 0;
    m_use_hydro_particles = 0;

    pp.queryAdd(std::string("use_"+REMORAParticleNames::tracers).c_str(), m_use_tracer_particles);
    pp.queryAdd(std::string("use_"+REMORAParticleNames::hydro).c_str(), m_use_hydro_particles);

    if (m_use_tracer_particles) {
        particleData.addName(REMORAParticleNames::tracers);
    }

    if (m_use_hydro_particles) {
        particleData.addName(REMORAParticleNames::hydro);
    }
    return;
}

/*! Initialize tracer and hydro particles */
void REMORA::initializeTracers ( ParGDBBase* a_gdb,
                              const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    auto& namelist_unalloc( particleData.getNamesUnalloc() );

    for (auto it = namelist_unalloc.begin(); it != namelist_unalloc.end(); ++it) {

        std::string species_name( *it );

        if (species_name == REMORAParticleNames::tracers) {

            AMREX_ASSERT(m_use_tracer_particles);
            REMORAPC* pc = new REMORAPC(a_gdb, REMORAParticleNames::tracers);
            pc->InitializeParticles(a_z_phys_nd[0]);
            amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " tracer particles.\n";
            particleData.pushBack(REMORAParticleNames::tracers, pc);

        } else if (species_name == REMORAParticleNames::hydro) {

            AMREX_ASSERT(m_use_hydro_particles);
            REMORAPC* pc = new REMORAPC(a_gdb, REMORAParticleNames::hydro);
            pc->InitializeParticles(a_z_phys_nd[0]);
            amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " hydro particles.\n";
            particleData.pushBack(REMORAParticleNames::hydro, pc);

        }
    }

    if (m_use_tracer_particles) namelist_unalloc.remove( REMORAParticleNames::tracers );
    if (m_use_hydro_particles)  namelist_unalloc.remove( REMORAParticleNames::hydro );

    return;
}

/*! The particle species, in a fixed order */
static Vector<REMORAPC*> particle_species (ParticleData& a_data, bool a_tracers, bool a_hydro)
{
    Vector<REMORAPC*> species;
    if (a_tracers) { species.push_back(a_data[REMORAParticleNames::tracers]); }
    if (a_hydro)   { species.push_back(a_data[REMORAParticleNames::hydro]); }
    return species;
}

/**
 * Whether each species moves with the flow on level a_lev, and whether any does. One with
 * no particles there needs no flow fields: particles are seeded on level 0 and usually stay
 * there, so this skips them on the finer levels. The count is a global reduction, so every
 * rank agrees.
 */
static bool advecting_species (const Vector<REMORAPC*>& a_species, int a_lev,
                               Vector<int>& a_advect)
{
    a_advect.resize(a_species.size());
    bool any_advect = false;
    for (int n = 0; n < a_species.size(); ++n) {
        a_advect[n] = a_species[n]->advectsWithFlow()
                   && a_species[n]->NumberOfParticlesAtLevel(a_lev) > 0;
        any_advect = any_advect || a_advect[n];
    }
    return any_advect;
}

/**
 * Particles are advected in terrain-following coordinates, so they need omega, the velocity
 * across the w-levels, rather than a Cartesian vertical velocity. As in ROMS floats, it comes
 * from the time-averaged fluxes update_massflux_3d leaves, the ones the tracers move with,
 * rather than the diagnostic omega from the current velocities. advance_3d integrates those
 * fluxes into W, the same volume flux compute_omega would build, its first ghost column
 * included, so this keeps W scaled by pm*pn instead of integrating the columns again. The
 * scaling and fill are those of compute_omega, so the result is the same to the bit.
 *
 * @param[in] lev   level
 * @param[in] mf_W  advance_3d's W, on the valid region and the first ghost column
 */
void REMORA::store_particle_omega (int lev, const MultiFab& mf_W)
{
    if (!m_use_tracer_particles && !m_use_hydro_particles) { return; }

    if (int(particle_omega_stored.size()) <= lev) {
        particle_omega_stored.resize(lev+1, 0);
        particle_advect.resize(lev+1);
    }
    const Vector<REMORAPC*> species = particle_species(particleData, m_use_tracer_particles,
                                                       m_use_hydro_particles);
    particle_omega_stored[lev] = 1;
    if (!advecting_species(species, lev, particle_advect[lev])) { return; }

    if (int(particle_omega.size()) <= lev) { particle_omega.resize(lev+1); }
    const BoxArray ba = convert(grids[lev],IntVect(0,0,1));
    const IntVect ng(NGROW,NGROW,0);
    auto& mf = particle_omega[lev];
    if (!mf || !mf->boxArray().CellEqual(ba) || mf->boxArray().ixType() != ba.ixType()
            || mf->DistributionMap() != dmap[lev] || mf->nGrowVect() != ng) {
        mf = std::make_unique<MultiFab>(ba, dmap[lev], 1, ng);
    }
    MultiFab& mf_omega = *mf;
    const int N = Geom(lev).Domain().size()[2]-1;

    mf_omega.setVal(zero);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_omega, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& omega = mf_omega.array(mfi);
        Array4<Real const> const& W     = mf_W.const_array(mfi);
        Array4<Real const> const& pm    = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn    = vec_pn[lev]->const_array(mfi);

        // The tile, with the first ghost column where it touches the edge of the grid
        const Box bx2d = makeSlab(mfi.growntilebox(IntVect(1,1,0)), 2, 0);
        ParallelFor(bx2d, N, [=] AMREX_GPU_DEVICE (int i, int j, int, int kk)
        {
            omega(i,j,kk+1) = W(i,j,kk+1) * pm(i,j,0) * pn(i,j,0);
        });
    }

    mf_omega.FillBoundary(geom[lev].periodicity());
}

/*! Evolve tracers and hydro particles for one time step*/
void REMORA::evolveTracers (int a_lev, Real a_dt_lev)
{
    if (!m_use_tracer_particles && !m_use_hydro_particles) { return; }

    const Vector<REMORAPC*> species = particle_species(particleData, m_use_tracer_particles,
                                                       m_use_hydro_particles);

    // Normally advance_3d has just stored omega and the species flags (store_particle_omega).
    // Otherwise compute them here.
    const bool stored = a_lev < int(particle_omega_stored.size()) && particle_omega_stored[a_lev];
    Vector<int> advect_here;
    bool any_advect;
    if (stored) {
        advect_here = particle_advect[a_lev];
        any_advect = false;
        for (int a : advect_here) { any_advect = any_advect || a; }
        particle_omega_stored[a_lev] = 0;
    } else {
        any_advect = advecting_species(species, a_lev, advect_here);
    }

    if (any_advect) {
        // Kept between steps rather than allocated on each one, and redefined when the
        // level's grids change. store_particle_omega (or compute_omega) and the fills below
        // overwrite them entirely.
        if (int(particle_omega.size()) <= a_lev) { particle_omega.resize(a_lev+1); }
        if (int(particle_u.size()) <= a_lev) {
            particle_u.resize(a_lev+1);
            particle_v.resize(a_lev+1);
        }
        auto fit = [&] (std::unique_ptr<MultiFab>& mf, const BoxArray& ba, const IntVect& ng)
        {
            if (!mf || !mf->boxArray().CellEqual(ba) || mf->boxArray().ixType() != ba.ixType()
                    || mf->DistributionMap() != dmap[a_lev] || mf->nGrowVect() != ng) {
                mf = std::make_unique<MultiFab>(ba, dmap[a_lev], 1, ng);
            }
        };
        if (!stored) {
            fit(particle_omega[a_lev], convert(grids[a_lev],IntVect(0,0,1)), IntVect(NGROW,NGROW,0));
            compute_omega(a_lev, *particle_omega[a_lev], *vec_Huon[a_lev], *vec_Hvom[a_lev]);
        }
        fit(particle_u[a_lev], xvel_new[a_lev]->boxArray(), xvel_new[a_lev]->nGrowVect());
        fit(particle_v[a_lev], yvel_new[a_lev]->boxArray(), yvel_new[a_lev]->nGrowVect());
        MultiFab& mf_u     = *particle_u[a_lev];
        MultiFab& mf_v     = *particle_v[a_lev];

        // A particle within half a cell of a wall or a coarse/fine boundary interpolates u and
        // v from ghost cells, which so far this step have only been filled across grids and
        // periodic boundaries. Fill copies with the physical and coarse/fine boundary
        // conditions at the end of the step. The model's own velocity ghosts are left alone:
        // refilling them with the physical boundary conditions here changes the answer on
        // open boundaries (see timeStepML).
        //
        // As at initialization and regrid, each copy is its own "calc" field, with no time
        // step: the NetCDF radiation conditions read the calc field unguarded. With calc
        // equal to the field, a radiation (orlanski_rad*) boundary leaves the ghost as it
        // was, so there the particles still see the start-of-step ghost value.
        MultiFab::Copy(mf_u, *xvel_new[a_lev], 0, 0, 1, mf_u.nGrowVect());
        MultiFab::Copy(mf_v, *yvel_new[a_lev], 0, 0, 1, mf_v.nGrowVect());
        FillPatch(a_lev, t_new[a_lev], mf_u, xvel_new, xvel_bc(), BdyVars::u, 0, true, false,
                  0, 0, zero, mf_u, xvel_old, xvel_new);
        FillPatch(a_lev, t_new[a_lev], mf_v, yvel_new, yvel_bc(), BdyVars::v, 0, true, false,
                  0, 0, zero, mf_v, yvel_old, yvel_new);
    }

    for (int n = 0; n < species.size(); ++n) {
        if (advect_here[n]) {
            Vector<MultiFab const*> flowvel = {particle_u[a_lev].get(), particle_v[a_lev].get(),
                                               particle_omega[a_lev].get()};
            species[n]->EvolveParticles(a_lev, a_dt_lev, flowvel, vec_z_phys_nd,
                                        *vec_pm[a_lev], *vec_pn[a_lev]);
        } else {
            // Nothing to advect on this level. A species that does not advect still keeps
            // its level index, so its heights follow the levels; then redistribute, as
            // EvolveParticles does.
            if (!species[n]->advectsWithFlow()) {
                species[n]->SetHeightsFromLevelIndex(a_lev, a_dt_lev, vec_z_phys_nd[a_lev]);
            }
            species[n]->Redistribute();
        }
    }
    return;
}

#endif

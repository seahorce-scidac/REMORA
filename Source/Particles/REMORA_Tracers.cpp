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

/*! Evolve tracers and hydro particles for one time step*/
void REMORA::evolveTracers (int a_lev, Real a_dt_lev)
{
    if (!m_use_tracer_particles && !m_use_hydro_particles) { return; }

    Vector<REMORAPC*> species;
    if (m_use_tracer_particles) { species.push_back(particleData[REMORAParticleNames::tracers]); }
    if (m_use_hydro_particles)  { species.push_back(particleData[REMORAParticleNames::hydro]); }

    // Whether each species moves with the flow on this level. One with no particles here needs
    // no flow fields: particles are seeded on level 0 and usually stay there, so this skips
    // them on the finer levels. The count is a global reduction, so every rank agrees.
    Vector<int> advect_here(species.size());
    bool any_advect = false;
    for (int n = 0; n < species.size(); ++n) {
        advect_here[n] = species[n]->advectsWithFlow()
                      && species[n]->NumberOfParticlesAtLevel(a_lev) > 0;
        any_advect = any_advect || advect_here[n];
    }

    // Particles are advected in terrain-following coordinates, so they need omega, the
    // velocity across the w-levels, rather than a Cartesian vertical velocity. As in ROMS
    // floats, it comes from the time-averaged fluxes update_massflux_3d leaves, the ones the
    // tracers move with, rather than the diagnostic omega from the current velocities.
    if (any_advect) {
        // Kept between steps rather than allocated on each one, and redefined when the
        // level's grids change. compute_omega and the fills below overwrite them entirely.
        if (int(particle_omega.size()) <= a_lev) {
            particle_omega.resize(a_lev+1);
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
        fit(particle_omega[a_lev], convert(grids[a_lev],IntVect(0,0,1)), IntVect(NGROW,NGROW,0));
        fit(particle_u[a_lev], xvel_new[a_lev]->boxArray(), xvel_new[a_lev]->nGrowVect());
        fit(particle_v[a_lev], yvel_new[a_lev]->boxArray(), yvel_new[a_lev]->nGrowVect());
        MultiFab& mf_omega = *particle_omega[a_lev];
        MultiFab& mf_u     = *particle_u[a_lev];
        MultiFab& mf_v     = *particle_v[a_lev];

        compute_omega(a_lev, mf_omega, *vec_Huon[a_lev], *vec_Hvom[a_lev]);

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

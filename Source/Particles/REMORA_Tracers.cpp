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

    pp.query(std::string("use_"+REMORAParticleNames::tracers).c_str(), m_use_tracer_particles);
    pp.query(std::string("use_"+REMORAParticleNames::hydro).c_str(), m_use_hydro_particles);

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
void REMORA::evolveTracers ( int                                        a_lev,
                             Real                                       a_dt_lev,
                             Vector<MultiFab const*>&                   a_flowvel,
                             const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    if (m_use_tracer_particles) {
      particleData[REMORAParticleNames::tracers]->EvolveParticles( a_lev,
                                                                   a_dt_lev,
                                                                   a_flowvel,
                                                                   a_z_phys_nd );
    }
    if (m_use_hydro_particles) {
      particleData[REMORAParticleNames::hydro]->EvolveParticles( a_lev,
                                                                 a_dt_lev,
                                                                 a_flowvel,
                                                                 a_z_phys_nd );
    }
    return;
}

#endif

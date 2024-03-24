#ifdef REMORA_USE_PARTICLES

#include <REMORA_PC.H>
#include <IndexDefines.H>
#include <REMORA_Constants.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

/*! Evolve particles for one time step */
void REMORAPC::EvolveParticles ( int                                        a_lev,
                                 Real                                       a_dt_lev,
                                 Vector<MultiFab const*>&                   a_flow_vel,
                                 const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("REMORAPCPC::EvolveParticles()");

    if (m_advect_w_flow) {
        AdvectWithFlow( a_lev, a_dt_lev, a_flow_vel, a_z_phys_nd[a_lev] );
    }

    Redistribute();
    return;
}

/*! Uses midpoint method to advance particles using flow velocity. */
void REMORAPC::AdvectWithFlow ( int                                 a_lev,
                                Real                                a_dt,
                                Vector<MultiFab const*>&            a_umac,
                                const std::unique_ptr<MultiFab>&    a_z_height )
{
    BL_PROFILE("REMORAPCPC::AdvectWithUmac()");
    AMREX_ASSERT(OK(a_lev, a_lev, a_umac[0]->nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());

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

            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(REMORAParticlesRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(REMORAParticlesRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(REMORAParticlesRealIdxSoA::vz).data();

            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*a_umac[0])[grid]),
                                                                  &((*a_umac[1])[grid]),
                                                                  &((*a_umac[2])[grid])) };

            // Array of these pointers to pass to the GPU
            GpuArray<Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {fab[0]->array(), fab[1]->array(), fab[2]->array()};

            bool use_terrain = (a_z_height != nullptr);
            auto zheight = use_terrain ? (*a_z_height)[grid].array() : Array4<Real>{};

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                if (use_terrain) {
                    mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
                } else {
                    mac_interpolate(p, plo, dxi, umacarr, v);
                }

                if (ipass == 0) {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        v_ptr[dim][i] = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);
                } else {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = v_ptr[dim][i] + static_cast<ParticleReal>(a_dt*v[dim]);
                        v_ptr[dim][i] = v[dim];
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);
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

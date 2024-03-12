#include "TracerPC.H"

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

void
TracerPC::
InitParticles (const MultiFab& a_z_height)
{
    BL_PROFILE("TracerPC::InitParticles");

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const auto& height = a_z_height[mfi];
        const FArrayBox* height_ptr = nullptr;
#ifdef AMREX_USE_GPU
        std::unique_ptr<FArrayBox> hostfab;
        if (height.arena()->isManaged() || height.arena()->isDevice()) {
            hostfab = std::make_unique<FArrayBox>(height.box(), height.nComp(),
                                                  The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(hostfab->dataPtr(), height.dataPtr(),
                                   height.size()*sizeof(Real));
            Gpu::streamSynchronize();
            height_ptr = hostfab.get();
        }
#else
        height_ptr = &height;
#endif
        Gpu::HostVector<ParticleType> host_particles;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            if (iv[0] == 3) {
                // Real r[3] = {0.5_rt, 0.5_rt, 0.4_rt};  // this means place just below the cell center
                // Real r[3] = {0.5_rt, 0.5_rt, 0.6_rt};  // this means place just above the cell center
                Real r[3] = {0.5_rt, 0.5_rt, 0.5_rt};  // this means place just at the cell center

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = (*height_ptr)(iv) + r[2]*((*height_ptr)(iv + IntVect(AMREX_D_DECL(0, 0, 1))) - (*height_ptr)(iv));

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;

                p.rdata(TracerRealIdx::old_x) = p.pos(0);
                p.rdata(TracerRealIdx::old_y) = p.pos(1);
                p.rdata(TracerRealIdx::old_z) = p.pos(2);

                p.idata(TracerIntIdx::k) = iv[2];  // particles carry their z-index

                host_particles.push_back(p);
           }
        }

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
}

/*
  /brief Uses midpoint method to advance particles using umac.
*/
void
TracerPC::AdvectWithUmac (Array<MultiFab const*, AMREX_SPACEDIM> umac,
                          int lev, Real dt, MultiFab& a_z_height)
{
    BL_PROFILE("TracerPC::AdvectWithUmac()");
    AMREX_ASSERT(OK(lev, lev, umac[0]->nGrow()-1));
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(umac[0]->nGrow() >= 1);,
                 AMREX_ASSERT(umac[1]->nGrow() >= 1);,
                 AMREX_ASSERT(umac[2]->nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!umac[0]->contains_nan());,
                 AMREX_ASSERT(!umac[1]->contains_nan());,
                 AMREX_ASSERT(!umac[2]->contains_nan()););

    const auto      strttime = amrex::second();
    const Geometry& geom = m_gdb->Geom(lev);
    const Box& domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto dx  = geom.CellSizeArray();

    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab const*> umac_pointer(AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(OnSameGrids(lev, *umac[0]));

     for (int i = 0; i < AMREX_SPACEDIM; i++) {
         umac_pointer[i] = umac[i];
     }

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };

            //array of these pointers to pass to the GPU
            GpuArray<Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                             (*fab[1]).array(),
                                             (*fab[2]).array() )}};

            auto const& zheight = a_z_height[grid].array();

            // We set this flag because the velocity arrays do not extend outside
            //    the domain in the z-direction
            int protect_against_out_of_bounds_in_vert = 1;

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                ParticleReal v[AMREX_SPACEDIM];

                // Interpolate the velocity from faces to the particle location
                mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v,
                                         protect_against_out_of_bounds_in_vert);

                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5_rt)*dt*v[dim]);
                    }

                    // Update z-coordinate stored in p.idata(0)
                    update_location_idata(p,plo,dxi,zheight);
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }

                    // Update z-coordinate stored in p.idata(0)
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

                amrex::Print() << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
#ifdef AMREX_LAZY
        });
#endif
    }
}

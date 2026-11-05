#ifdef REMORA_USE_PARTICLES

#include <AMReX_ParticleHeader.H>
#include <AMReX_ParticleMesh.H>
#include <REMORA_Constants.H>
#include <REMORA_PC.H>

using namespace amrex;

/**
 * Mass per unit physical volume, by cloud-in-cell deposit in index space. Horizontally the
 * weights come from the particle's position in the uniform coordinates it moves in, and
 * vertically from its level index s (cell k spans s in [k, k+1]), so the deposit lands in
 * the cells around the particle whatever the stretching of the levels. Weight that would
 * fall below the bed, above the free surface or across a non-periodic side of the domain
 * goes to the edge cell instead, so no mass is lost. The deposited mass is then divided by
 * each cell's volume, Hz/(pm*pn).
 *
 * @param a_mf     Multifab to calculate over
 * @param a_lev    level to calculate on
 * @param a_Hz     cell thickness on this level
 * @param a_pm     inverse grid spacing in x on this level
 * @param a_pn     inverse grid spacing in y on this level
 * @param a_comp   component to calculate
 */
void REMORAPC::massDensity ( MultiFab&  a_mf,
                             const int& a_lev,
                             MultiFab const& a_Hz,
                             MultiFab const& a_pm,
                             MultiFab const& a_pn,
                             const int& a_comp ) const
{
    BL_PROFILE("REMORAPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto dlo = lbound(geom.Domain());
    const auto dhi = ubound(geom.Domain());
    const bool clamp_x = !geom.isPeriodic(0);
    const bool clamp_y = !geom.isPeriodic(1);

    a_mf.setVal(zero);

    ParticleToMesh( *this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (  const REMORAPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            const auto& p = ptd.m_aos[i];
            const auto mass = static_cast<Real>(ptd.m_rdata[REMORAParticlesRealIdxSoA::mass][i]);

            // Lower cell of the stencil and the weight of the upper one, from the
            // fractional cell-centred index a_xc
            auto cic = [] (Real a_xc, int a_lo, int a_hi, bool a_clamp, int& a_i0, Real& a_w1)
            {
                a_i0 = static_cast<int>(amrex::Math::floor(a_xc));
                a_w1 = a_xc - static_cast<Real>(a_i0);
                if (a_clamp) {
                    if (a_i0 < a_lo) {
                        a_i0 = a_lo; a_w1 = zero;
                    } else if (a_i0 >= a_hi) {
                        a_i0 = a_hi; a_w1 = zero;
                    }
                }
            };

            int i0, j0, k0;
            Real wx, wy, wz;
            cic((static_cast<Real>(p.pos(0)) - plo[0])*dxi[0] - myhalf, dlo.x, dhi.x, clamp_x, i0, wx);
            cic((static_cast<Real>(p.pos(1)) - plo[1])*dxi[1] - myhalf, dlo.y, dhi.y, clamp_y, j0, wy);
            cic(static_cast<Real>(ptd.m_rdata[REMORAParticlesRealIdxSoA::s][i]) - myhalf,
                dlo.z, dhi.z, true, k0, wz);

            for (int kk = 0; kk <= 1; ++kk) {
                const Real fz = kk ? wz : one - wz;
                if (fz == zero) { continue; }
                for (int jj = 0; jj <= 1; ++jj) {
                    const Real fy = jj ? wy : one - wy;
                    if (fy == zero) { continue; }
                    for (int ii = 0; ii <= 1; ++ii) {
                        const Real fx = ii ? wx : one - wx;
                        if (fx == zero) { continue; }
                        Gpu::Atomic::AddNoRet(&rho(i0+ii, j0+jj, k0+kk, a_comp), fx*fy*fz*mass);
                    }
                }
            }
        });

    // Mass per cell to mass per unit volume
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real> rho = a_mf.array(mfi);
        const Array4<Real const> Hz = a_Hz.const_array(mfi);
        const Array4<Real const> pm = a_pm.const_array(mfi);
        const Array4<Real const> pn = a_pn.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho(i,j,k,a_comp) = (Hz(i,j,k) > zero)
                              ? rho(i,j,k,a_comp) * pm(i,j,0) * pn(i,j,0) / Hz(i,j,k) : zero;
        });
    }

    return;
}

/**
 * Checkpoints written before the level index s was stored have one real component fewer.
 * Such a checkpoint is read into a container with that layout, the particles are copied into
 * this one, and s is computed from each particle's height on a_z_phys_nd, the grid read from
 * the same checkpoint. That is the index the particle's next step would have recovered.
 *
 * @param a_dir         checkpoint directory
 * @param a_name        particle species name
 * @param a_z_phys_nd   z nodal positions of the grid, per level
 */
void REMORAPC::Restart ( const std::string& a_dir, const std::string& a_name,
                         const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    BL_PROFILE("REMORAPC::Restart()");

    using BasePC   = ParticleContainer<REMORAParticlesRealIdxAoS::ncomps, REMORAParticlesIntIdxAoS::ncomps,
                                       REMORAParticlesRealIdxSoA::ncomps, REMORAParticlesIntIdxSoA::ncomps,
                                       DefaultAllocator, REMORAParticlesAssignor>;
    // The layout before s was added: the same components, with s missing at the end
    constexpr int nreal_old = REMORAParticlesRealIdxSoA::s;
    using LegacyPC = ParticleContainer<REMORAParticlesRealIdxAoS::ncomps, REMORAParticlesIntIdxAoS::ncomps,
                                       nreal_old, REMORAParticlesIntIdxSoA::ncomps,
                                       DefaultAllocator, REMORAParticlesAssignor>;

    Vector<char> header_chars;
    ParallelDescriptor::ReadAndBcastFile(a_dir + "/" + a_name + "/Header", header_chars);
    std::istringstream header_stream(std::string(header_chars.dataPtr()), std::istringstream::in);
    ParticleHeader header;
    header.parse(header_stream);

    if (header.num_real == REMORAParticlesRealIdxSoA::ncomps) {
        BasePC::Restart(a_dir, a_name);
        return;
    }
    if (header.num_real != nreal_old) {
        Abort("REMORAPC::Restart(): " + a_name + " in " + a_dir + " has "
              + std::to_string(header.num_real) + " real components; expected "
              + std::to_string(REMORAParticlesRealIdxSoA::ncomps) + ", or "
              + std::to_string(nreal_old) + " for a checkpoint written before level_index");
    }
    Print() << "Restarting " << a_name << " from a checkpoint without level_index;"
            << " computing it from the particle heights\n";

    LegacyPC old_pc(m_gdb);
    old_pc.Restart(a_dir, a_name);

    resizeData();
    for (int lev = 0; lev <= old_pc.finestLevel(); ++lev)
    {
        AMREX_ALWAYS_ASSERT(a_z_phys_nd[lev] != nullptr);
        const Geometry& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const int nz = geom.Domain().length(2);

        for (const auto& kv : old_pc.GetParticles(lev))
        {
            const auto& src = kv.second;
            const int np = src.numParticles();
            auto& dst = DefineAndReturnParticleTile(lev, kv.first.first, kv.first.second);
            dst.resize(np);
            if (np == 0) { continue; }

            const auto& src_aos = src.GetArrayOfStructs();
            Gpu::copyAsync(Gpu::deviceToDevice, src_aos.begin(), src_aos.end(),
                           dst.GetArrayOfStructs().begin());
            for (int comp = 0; comp < nreal_old; ++comp) {
                const auto& src_comp = src.GetStructOfArrays().GetRealData(comp);
                Gpu::copyAsync(Gpu::deviceToDevice, src_comp.begin(), src_comp.end(),
                               dst.GetStructOfArrays().GetRealData(comp).begin());
            }

            auto* p_pbox = dst.GetArrayOfStructs()().data();
            auto* s_ptr = dst.GetStructOfArrays().GetRealData(REMORAParticlesRealIdxSoA::s).data();
            const auto zheight = (*a_z_phys_nd[lev])[kv.first.first].const_array();
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                const ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { s_ptr[i] = zero; return; }
                s_ptr[i] = static_cast<ParticleReal>(
                    height_to_level_index(p.pos(0), p.pos(1), p.pos(2),
                                          p.idata(REMORAParticlesIntIdxAoS::k),
                                          nz, plo, dxi, zheight).s);
            });
        }
    }
    Gpu::streamSynchronize();
}

/**
 * Write the particles of one species, for a checkpoint or a plotfile. With a_xy null they are
 * written as they are, in the uniform coordinates of the AMReX geometry they move in. With a_xy
 * given (plotfiles), x and y are replaced by the physical coordinates of the cell corners (psi
 * points) interpolated bilinearly to each particle. These are the corner positions the
 * plotfile's nodal displacement gives a viewer, so each particle is drawn inside the cell that
 * holds it. z is already the physical height. The uniform x and y are saved before the
 * conversion and restored bit for bit after the write, so writing a plotfile does not change
 * the run.
 *
 * Each particle stays filed under the grid that physically holds it, though its physical
 * position need not map to that grid through the uniform geometry, so pc.OK() does not hold
 * during the write. AMReX's particle writers do not require it (AMReX-Codes/amrex#6175).
 *
 * @param a_fname   checkpoint or plotfile directory
 * @param a_name    particle species name
 * @param a_xy      null, or per level the physical x and y of the psi points
 *                  (REMORA::plot_node_coords)
 */
void REMORAPC::WriteOutput ( const std::string& a_fname, const std::string& a_name,
                             const Vector<MultiFab>* a_xy )
{
    BL_PROFILE("REMORAPC::WriteOutput()");

    if (a_xy == nullptr) {
        Checkpoint(a_fname, a_name, true, varNames());
        return;
    }

    // The uniform x and y of each tile, keyed like the particle tiles by (grid, tile). Every
    // entry is made here, so the threaded loops below only look entries up.
    Vector<std::map<std::pair<int,int>, Gpu::DeviceVector<ParticleReal>>> saved_xy(finestLevel()+1);
    for (int lev = 0; lev <= finestLevel(); ++lev) {
        for (const auto& kv : GetParticles(lev)) {
            saved_xy[lev][kv.first].resize(2 * kv.second.numParticles());
        }
    }

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        const Geometry& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& aos = pti.GetArrayOfStructs();
            const int np = aos.numParticles();
            auto* p_pbox = aos().data();
            auto* p_save = saved_xy[lev].at(pti.GetPairIndex()).data();
            const auto xy = (*a_xy)[lev].const_array(pti.index());
            const Box& vbx = pti.validbox();
            const int ilo = vbx.smallEnd(0);
            const int ihi = vbx.bigEnd(0);
            const int jlo = vbx.smallEnd(1);
            const int jhi = vbx.bigEnd(1);

            ParallelFor(np, [=] AMREX_GPU_DEVICE (int n)
            {
                auto& p = p_pbox[n];
                p_save[2*n  ] = p.pos(0);
                p_save[2*n+1] = p.pos(1);
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

    Checkpoint(a_fname, a_name, true, varNames());

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& aos = pti.GetArrayOfStructs();
            const int np = aos.numParticles();
            auto* p_pbox = aos().data();
            const auto* p_save = saved_xy[lev].at(pti.GetPairIndex()).data();
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int n)
            {
                p_pbox[n].pos(0) = p_save[2*n  ];
                p_pbox[n].pos(1) = p_save[2*n+1];
            });
        }
    }
    // saved_xy is freed on return. The restore kernels may run on several streams, so wait for
    // the whole device rather than relying on ParIter's own synchronization.
    Gpu::synchronize();
}

#endif

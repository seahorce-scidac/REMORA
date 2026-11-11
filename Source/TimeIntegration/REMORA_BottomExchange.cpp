#include <REMORA.H>
#include <REMORA_prob_common.H>
#include <REMORA_BottomExchange.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Reduce.H>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

/** Return the level-0 rho layout accepted by the coupler.
 * @param[out  ] ba horizontal cell layout, with one cell in z
 * @param[out  ] dm distribution of those cells across ranks
 */
void
REMORA::GetBottomExchangeLayout (BoxArray& ba, DistributionMapping& dm) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level >= 0, "Initialize REMORA before querying bottom layout");
    ba = vec_mskr[0]->boxArray();
    dm = vec_mskr[0]->DistributionMap();
}

/** Copy a coupler snapshot of water flux and incoming concentrations.
 * @param[in   ] q water flux per horizontal area (m/s), positive into the ocean
 * @param[in   ] concentrations incoming values in GetTracerNames() order
 * @param[in   ] begin beginning of validity in elapsed model seconds
 * @param[in   ] end end of validity in elapsed model seconds
 */
void
REMORA::SetBottomExchangeData (const MultiFab& q, const MultiFab& concentrations,
                             Real begin, Real end)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.bottom_exchange_type == BottomExchangeType::coupler,
        "SetBottomExchangeData requires remora.bottom_exchange_type=coupler");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(begin) && std::isfinite(end) && end >= begin,
        "Bottom exchange needs a finite validity interval");
    BoxArray ba;
    DistributionMapping dm;
    GetBottomExchangeLayout(ba, dm);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(q.boxArray() == ba && concentrations.boxArray() == ba &&
        q.nComp() == 1 && concentrations.nComp() == ncons,
        "Bottom input must cover the rho layout, with 1 water and ncons concentration components");
    m_bottom_coupled = std::make_unique<MultiFab>(ba, dm, ncons+1, IntVect(0));
    m_bottom_coupled->ParallelCopy(q, 0, 0, 1);
    m_bottom_coupled->ParallelCopy(concentrations, 0, 1, ncons);
    m_bottom_begin = begin;
    m_bottom_end = end;
    invalidate_bottom_exchange(0);
}

/** Export bottom-cell concentrations on the requested rho layout, with land zeroed.
 * @param[in   ] lev refinement level
 * @param[out  ] concentrations values in GetTracerNames() order
 */
void
REMORA::PackBottomState (int lev, MultiFab& concentrations) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev >= 0 && lev <= finest_level &&
        concentrations.boxArray() == vec_mskr[lev]->boxArray() && concentrations.nComp() == ncons,
        "PackBottomState needs the requested level's rho slab with ncons components");
    MultiFab bottom(vec_mskr[lev]->boxArray(), dmap[lev], ncons, 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(bottom, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto out = bottom.array(mfi);
        const auto state = cons_new[lev]->const_array(mfi);
        const auto mask = vec_mskr[lev]->const_array(mfi);
        ParallelFor(mfi.tilebox(), ncons, [=] AMREX_GPU_DEVICE (int i, int j, int, int n) {
            out(i,j,0,n) = mask(i,j,0) > 0 ? state(i,j,0,n) : Real(0);
        });
    }
    concentrations.ParallelCopy(bottom, 0, 0, ncons);
}

/** Return prepared integrated water and incoming tracer transports.
 * @param[in   ] lev refinement level
 * @return water (m^3/s), then tracer transports in registry order; null if unavailable
 */
const MultiFab*
REMORA::GetBottomExchangeTransport (int lev) const
{
    if (lev < 0 || lev >= int(m_bottom_transport.size()) || !m_bottom_transport[lev] ||
        m_bottom_transport[lev]->boxArray() != vec_mskr[lev]->boxArray() ||
        m_bottom_transport[lev]->DistributionMap() != dmap[lev]) { return nullptr; }
    return m_bottom_transport[lev].get();
}

/** Invalidate derived transports after new coupled data or a level-layout change.
 * @param[in   ] lev changed level; level 0 also invalidates the shared forcing surface
 */
void
REMORA::invalidate_bottom_exchange (int lev)
{
    if (lev == 0) {
        m_bottom_base_ready = false;
        m_bottom_transport.clear();
        m_bottom_level_times.clear();
    } else if (lev < int(m_bottom_transport.size())) {
        m_bottom_transport[lev].reset();
        if (lev < int(m_bottom_work.size())) { m_bottom_work[lev].reset(); }
    }
}

/** Sample the provider and conservatively distribute bed transports to a level.
 * @param[in   ] lev refinement level
 * @param[in   ] time beginning of the slow step in elapsed model seconds
 * @param[in   ] dt_lev slow-step duration (s); zero samples at time for diagnostics
 */
void
REMORA::prepare_bottom_exchange (int lev, Real time, Real dt_lev)
{
    if (solverChoice.bottom_exchange_type == BottomExchangeType::none) { return; }
    BL_PROFILE("REMORA::prepare_bottom_exchange()");
    if (solverChoice.bottom_exchange_type == BottomExchangeType::coupler) {
        const Real tol = Real(32)*std::numeric_limits<Real>::epsilon()*
                         amrex::max(Real(1), amrex::Math::abs(time+dt_lev));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_bottom_coupled && time >= m_bottom_begin-tol &&
            time+dt_lev <= m_bottom_end+tol,
            "Coupled bottom forcing is missing or expired; submit data covering the whole slow step"
            " (also after restart)");
    }
    const int nc = ncons + 1;
    const int salt_comp = component_index("salt") + 1;
    // compute_W and the free-surface update read the first horizontal ghost ring.
    const IntVect ng(1, 1, 0);
    const auto& ba = vec_mskr[0]->boxArray();
    const auto& dm = vec_mskr[0]->DistributionMap();
    if (!m_bottom_raw || m_bottom_raw->boxArray() != ba || m_bottom_raw->DistributionMap() != dm) {
        invalidate_bottom_exchange(0);
        m_bottom_work.clear();
        m_bottom_work_generation.clear();
        m_bottom_raw = std::make_unique<MultiFab>(ba, dm, nc, ng);
        m_bottom_base_transport = std::make_unique<MultiFab>(ba, dm, nc, ng);
#ifdef REMORA_USE_NETCDF
        m_bottom_series.clear();
        m_bottom_file_fields.clear();
#endif
    }
    const Real sample_time = time + Real(0.5)*dt_lev;
    const bool changing = solverChoice.bottom_exchange_type == BottomExchangeType::analytic ||
                          solverChoice.bottom_exchange_type == BottomExchangeType::netcdf;
    if (GetBottomExchangeTransport(lev) &&
        (!changing || m_bottom_level_times[lev] == sample_time)) { return; }
    Vector<BCRec> bcs(nc);
    for (auto& bc : bcs) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bc.setLo(d, geom[0].isPeriodic(d) ? BCType::int_dir : BCType::foextrap);
            bc.setHi(d, geom[0].isPeriodic(d) ? BCType::int_dir : BCType::foextrap);
        }
    }
    if (!m_bottom_base_ready || (changing && m_bottom_base_time != sample_time)) {
        auto& raw = *m_bottom_raw;
        raw.setVal(solverChoice.bottom_q, 0, 1, ng);
        for (int n = 0; n < ncons; ++n) {
            raw.setVal(solverChoice.bottom_cons[n], n + 1, 1, ng);
        }
        if (solverChoice.bottom_exchange_type == BottomExchangeType::analytic) {
            // One physical forcing surface, evaluated at each level's own slow-step time.
            prob->bottom_exchange(0, sample_time, geom[0], *this, raw);
        } else if (solverChoice.bottom_exchange_type == BottomExchangeType::coupler) {
            raw.ParallelCopy(*m_bottom_coupled, 0, 0, nc);
        }
#ifdef REMORA_USE_NETCDF
        else if (solverChoice.bottom_exchange_type == BottomExchangeType::netcdf) {
            if (m_bottom_series.empty()) {
                m_bottom_series.resize(nc);
                m_bottom_file_fields.resize(nc);
                for (int n = 0; n < nc; ++n) {
                    const std::string name = n == 0 ? "bottom_q" : "bottom_" + cons_names[n-1];
                    if (n > 0 && !NCTimeSeries::has_field(solverChoice.bottom_exchange_files, name)) {
                        continue;
                    }
                    m_bottom_file_fields[n] = std::make_unique<MultiFab>(ba, dm, 1, ng);
                    m_bottom_series[n] = std::make_unique<NCTimeSeries>(
                        solverChoice.bottom_exchange_files, name, solverChoice.bottom_time_name,
                        geom[0].Domain(), m_bottom_file_fields[n].get(), true, false,
                        solverChoice.time_ref);
                    m_bottom_series[n]->Initialize();
                }
            }
            for (int n = 0; n < nc; ++n) {
                if (!m_bottom_series[n]) { continue; }
                m_bottom_series[n]->update_interpolated_to_time(model_time(sample_time), 0,
                    m_bottom_file_fields[n].get(), geom, ref_ratio);
                MultiFab::Copy(raw, *m_bottom_file_fields[n], 0, n, 1, ng);
            }
        }
#endif
        Vector<MultiFab*> raw_levels{&raw};
        const int prescribed_ring = solverChoice.bottom_exchange_type == BottomExchangeType::netcdf
                                  ? 1 : 0;
        for (int n = 0; n < nc; ++n) {
            FillPatch(0, sample_time, raw, raw_levels, foextrap_periodic_bc(),
                      BdyVars::null, n, false, false, prescribed_ring);
        }
        // Validate only wet cells and incoming concentrations actually used. Never multiply
        // missing concentrations by zero: an outflow-only file needs no tracer variables.
        ReduceOps<ReduceOpMax> op;
        ReduceData<int> data(op);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(raw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = raw.const_array(mfi);
            const auto mask = vec_mskr[0]->const_array(mfi);
            op.eval(mfi.tilebox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> GpuTuple<int> {
                if (mask(i,j,0) <= 0) { return {0}; }
                if (!std::isfinite(a(i,j,0,0))) { return {1}; }
                if (a(i,j,0,0) > 0) {
                    for (int n = 1; n < nc; ++n) {
                        if (!std::isfinite(a(i,j,0,n))) { return {1}; }
                    }
                    if (a(i,j,0,salt_comp) < 0) { return {1}; }
                }
                return {0};
            });
        }
        int bad = amrex::get<0>(data.value(op));
        ParallelDescriptor::ReduceIntMax(bad);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!bad,
            "Bottom exchange needs finite q and explicit finite concentrations for every tracer"
            " on inflow; salinity must be nonnegative");
        auto& base = *m_bottom_base_transport;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(base, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = raw.const_array(mfi);
            const auto out = base.array(mfi);
            const auto mask = vec_mskr[0]->const_array(mfi);
            const auto pm = vec_pm[0]->const_array(mfi);
            const auto pn = vec_pn[0]->const_array(mfi);
            ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int) {
                const Real q = mask(i,j,0) > 0 ? a(i,j,0,0)/(pm(i,j,0)*pn(i,j,0)) : Real(0);
                out(i,j,0,0) = q;
                for (int n = 1; n < nc; ++n) { out(i,j,0,n) = q > 0 ? q*a(i,j,0,n) : Real(0); }
            });
        }
        m_bottom_base_ready = true;
        m_bottom_base_time = sample_time;
        ++m_bottom_generation;
    }
    // Fill one instantaneous hierarchy in level-0 transport units. The coarse
    // fields have this sample time, independently of the parent's held slow-step data.
    m_bottom_work.resize(finest_level + 1);
    m_bottom_work_generation.resize(finest_level + 1, -1);
    Vector<MultiFab*> fields(lev + 1, nullptr);
    fields[0] = m_bottom_base_transport.get();
    IntVect ratio(1);
    const Geometry coarse_geom(geom[0].Domain(), geom[0].ProbDomain(), geom[0].Coord(),
                               Array<int,AMREX_SPACEDIM>{0,0,0});
    for (int level = 1; level <= lev; ++level) {
        ratio *= ref_ratio[level - 1];
        auto& work = m_bottom_work[level];
        if (!work || work->boxArray() != vec_mskr[level]->boxArray() ||
            work->DistributionMap() != dmap[level]) {
            work = std::make_unique<MultiFab>(vec_mskr[level]->boxArray(), dmap[level], nc, ng);
            m_bottom_work_generation[level] = -1;
        }
        fields[level] = work.get();
        if (m_bottom_work_generation[level] != m_bottom_generation) {
            InterpFromCoarseLevel(*work, ng, ng, *fields[0], 0, 0, nc, coarse_geom,
                                 geom[level], ratio, &pc_interp, bcs, 0);
            for (int n = 0; n < nc; ++n) {
                FillPatch(level, sample_time, *work, fields, foextrap_periodic_bc(),
                          BdyVars::null, n, false, false, 1, 0, Real(0), MultiFab(),
                          {}, {}, -1, &pc_interp);
            }
            m_bottom_work_generation[level] = m_bottom_generation;
        }
    }

    m_bottom_transport.resize(finest_level + 1);
    m_bottom_level_times.resize(finest_level + 1);
    auto& target = m_bottom_transport[lev];
    if (!target || target->boxArray() != vec_mskr[lev]->boxArray() ||
        target->DistributionMap() != dmap[lev]) {
        target = std::make_unique<MultiFab>(vec_mskr[lev]->boxArray(), dmap[lev], nc, ng);
    }
    MultiFab::Copy(*target, *fields[lev], 0, 0, nc, ng);
    if (lev > 0) {
        // Divide integrated Q and Q*C among descendant cells after filling halos.
        const Real fraction = Real(1) / Real(ratio[0] * ratio[1]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*target, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = target->array(mfi);
            ParallelFor(mfi.growntilebox(), nc,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                a(i,j,k,n) *= fraction;
            });
        }
        // A finer land mask must not silently delete part of a water source.
        ReduceOps<ReduceOpMax> mask_op;
        ReduceData<int> mask_data(mask_op);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*target, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = target->const_array(mfi), mask = vec_mskr[lev]->const_array(mfi);
            mask_op.eval(mfi.tilebox(), mask_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int) -> GpuTuple<int> {
                    return {mask(i,j,0) <= 0 && a(i,j,0) != 0 ? 1 : 0};
                });
        }
        int inconsistent = get<0>(mask_data.value(mask_op));
        ParallelDescriptor::ReduceIntMax(inconsistent);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!inconsistent,
            "Bottom exchange requires consistent coarse/fine wet masks over source cells");
    }
    m_bottom_level_times[lev] = sample_time;
}

/** Check the source Courant limit and optionally report mechanical GLS measures.
 * @param[in   ] lev refinement level
 * @param[in   ] dt_lev slow-step duration (s)
 */
void
REMORA::check_bottom_exchange (int lev, Real dt_lev) const
{
    const MultiFab* mf = GetBottomExchangeTransport(lev);
    if (!mf) { return; }
    const bool report = solverChoice.bottom_exchange_diagnostics && verbose > 0;
    const bool gls = solverChoice.vert_mixing_type == VertMixingType::GLS;
    const int record = istep[lev] % 2;
    const Real cmu3 = std::pow(solverChoice.gls_cmu0, Real(3));
    ReduceOps<ReduceOpMax,ReduceOpMax,ReduceOpMax,ReduceOpMax> op;
    ReduceData<Real,Real,Real,Real> data(op);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto q = mf->const_array(mfi);
        const auto pm = vec_pm[lev]->const_array(mfi);
        const auto pn = vec_pn[lev]->const_array(mfi);
        const auto hz = vec_Hz[lev]->const_array(mfi);
        const auto mask = vec_mskr[lev]->const_array(mfi);
        const auto km = vec_Akv[lev]->const_array(mfi);
        const auto kk = gls ? vec_Akk[lev]->const_array(mfi) : Array4<Real const>();
        const auto kp = gls ? vec_Akp[lev]->const_array(mfi) : Array4<Real const>();
        const auto tke = gls ? vec_tke[lev]->const_array(mfi) : Array4<Real const>();
        const auto len = gls ? vec_Lscale[lev]->const_array(mfi) : Array4<Real const>();
        const auto u = xvel_old[lev]->const_array(mfi);
        const auto v = yvel_old[lev]->const_array(mfi);
        const auto tx = vec_bustr[lev]->const_array(mfi);
        const auto ty = vec_bvstr[lev]->const_array(mfi);
        op.eval(mfi.tilebox(), data,
        [=] AMREX_GPU_DEVICE (int i, int j, int) -> GpuTuple<Real,Real,Real,Real>
        {
            if (mask(i,j,0) <= 0) { return {0,0,0,0}; }
            const Real speed = std::abs(q(i,j,0))*pm(i,j,0)*pn(i,j,0);
            const Real courant = speed*dt_lev/hz(i,j,0);
            if (!report || !gls) { return {courant,0,0,0}; }
            const Real k = tke(i,j,1,record);
            const Real l = len(i,j,1);
            const Real eps = cmu3*k*std::sqrt(k)/amrex::max(l,Real(1.e-20));
            const Real ut = Real(0.5)*(u(i,j,0)+u(i+1,j,0));
            const Real vt = Real(0.5)*(v(i,j,0)+v(i,j+1,0));
            const Real stress = std::hypot(Real(0.5)*(tx(i,j,0)+tx(i+1,j,0)),
                                          Real(0.5)*(ty(i,j,0)+ty(i,j+1,0)));
            const auto ratios = BottomExchange::scales(speed, Real(0.5)*hz(i,j,0),
                km(i,j,1), kk(i,j,1), kp(i,j,1), k, eps, stress, std::hypot(ut,vt));
            return {courant,ratios.transport,ratios.adjustment,ratios.wall};
        });
    }
    const auto result = data.value(op);
    Real values[4] = {get<0>(result),get<1>(result),get<2>(result),get<3>(result)};
    ParallelDescriptor::ReduceRealMax(values,4);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(values[0]) &&
        values[0] <= solverChoice.bottom_exchange_max_courant,
        "Bottom exchange exceeds bottom_exchange_max_courant; reduce the slow timestep");
    if (report) {
        Print() << "Bottom exchange level " << lev << ": Courant=" << values[0];
        if (gls) {
            Print() << " transport=" << values[1] << " adjustment=" << values[2]
                    << " wall=" << values[3]
                    << " (wall undefined at zero stress; TKE floors are numerical)";
        }
        Print() << '\n';
    }
}

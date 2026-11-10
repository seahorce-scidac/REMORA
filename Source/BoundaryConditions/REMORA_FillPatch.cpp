#include <REMORA.H>
#include <REMORA_PhysBCFunct.H>
#include <REMORA_IndexDefines.H>
#include <REMORA_TimeInterpolatedData.H>

using namespace amrex;

PhysBCFunctNoOp null_bc;

namespace {
/** \brief dst(comp) = src(comp) wherever mask is set, ghost cells included */
void copy_on_mask (MultiFab& dst, const MultiFab& src, const iMultiFab& mask, int comp)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto d = dst.array(mfi, comp);
        auto a = src.const_array(mfi, comp);
        auto m = mask.const_array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (m(i,j,k)) { d(i,j,k) = a(i,j,k); }
        });
    }
}

/** \brief number of cells where mask is set and a(comp) and b(comp) differ in value
 *         (two NaNs count as equal), over all ranks */
Long count_mask_diff (const MultiFab& a, const MultiFab& b, const iMultiFab& mask, int comp)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(a); mfi.isValid(); ++mfi) {
        auto x = a.const_array(mfi, comp);
        auto y = b.const_array(mfi, comp);
        auto m = mask.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            // Outside the mask the cells may never have been written (a signalling NaN in a
            // debug build), so do not even compare them
            if (!m(i,j,k)) { return { Long(0) }; }
            const Real xv = x(i,j,k);
            const Real yv = y(i,j,k);
            const bool same = (xv == yv) || (xv != xv && yv != yv);
            return { same ? Long(0) : Long(1) };
        });
    }
    Long n = amrex::get<0>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceLongSum(n);
    return n;
}

/**
 * ba with each box extended outward by ng across the non-periodic domain faces it touches.
 * The result is kept and handed out again for the same ba, so a MultiFab built on it has the
 * same BoxArray from call to call and FillBoundary reuses its communication pattern rather
 * than rebuilding it on every fill. The few layouts a run uses at a time are kept; the cache
 * is emptied at finalize, before AMReX is torn down.
 */
BoxArray
extended_boxarray (const BoxArray& ba, const Geometry& geom, const IntVect& ng)
{
    struct Entry { BoxArray ba; Box dom; IntVect ng; BoxArray ext; };
    static Vector<Entry> cache;
    static bool registered = false;
    if (!registered) {
        amrex::ExecOnFinalize([] () { cache.clear(); });
        registered = true;
    }

    const Box dom = amrex::convert(geom.Domain(), ba.ixType());
    for (const auto& e : cache) {
        if (e.ng == ng && e.dom == dom && e.ba == ba) { return e.ext; }
    }

    BoxList bl(ba.ixType());
    for (int i = 0; i < ba.size(); ++i) {
        Box b = ba[i];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (geom.isPeriodic(d)) { continue; }
            if (b.smallEnd(d) <= dom.smallEnd(d)) { b.growLo(d, ng[d]); }
            if (b.bigEnd(d)   >= dom.bigEnd(d))   { b.growHi(d, ng[d]); }
        }
        bl.push_back(b);
    }
    BoxArray ext(std::move(bl));

    constexpr int max_entries = 16;
    if (cache.size() >= max_entries) { cache.erase(cache.begin()); }
    cache.push_back({ba, dom, ng, ext});
    return ext;
}

/**
 * A copy of mf whose ghost cells beyond a non-periodic domain face agree with their periodic
 * images, for a copy that reads ghost cells with periodicity.
 *
 * FillBoundary copies only from valid cells, so a ghost cell that lies outside the domain in
 * a non-periodic direction and beyond a periodic edge -- a corner where a periodic seam meets
 * a wall, or the vertical ghosts at a seam -- keeps whatever the boundary conditions put there,
 * which need not match the cell it images. Face data normal to a periodic direction has the
 * same problem one step earlier: below a wall, the faces at the two ends of the periodic width
 * are both owned by a box and need not agree either. InterpFromCoarseLevel copies coarse ghost
 * cells with periodicity, so such a cell then has two sources, the box's own ghost and its
 * periodic image. ParallelCopy does not say which one wins. On a CPU the order is fixed, but on
 * a GPU whichever thread gets there first wins, so the result varies from run to run (#768).
 *
 * Here each box is extended outward across its non-periodic domain faces, so those ghosts
 * become valid cells of a temporary. OverrideSync makes the shared periodic faces agree, and
 * FillBoundary then fills the images. The result goes into scratch rather than back into mf,
 * so the coarse level's own ghost cells are left as they were.
 *
 * Doing this with AMReX calls alone gets most of the way. EnforcePeriodicity takes its
 * sources from each box grown by its ghost cells -- Periodicity::Domain() is unbounded in the
 * non-periodic directions -- so on a copy of mf, FillBoundary followed by EnforcePeriodicity
 * already makes cell-centered ghosts beyond a wall match their images. On Channel_Test_ML_seam
 * that gives a plotfile bitwise identical to this function's. It does not cover face data
 * normal to a periodic direction: faces 0 and nx both lie inside the nodal periodic domain, so
 * EnforcePeriodicity never writes them, and OverrideSync reconciles only their valid copies.
 * Below a wall the two then still differ (up to 3e-15 for xvel in that test), leaving two
 * sources. The extended temporary makes those faces valid so OverrideSync reaches them.
 *
 * @param[in   ] mf       data to read, its ghost cells already FillBoundary'd
 * @param[in   ] geom     geometry of mf's level
 * @param[inout] scratch  storage for the copy, used only when one is needed
 * @return mf itself when no ghost cell can have two sources, otherwise scratch
 */
const MultiFab&
periodic_consistent_ghosts (const MultiFab& mf, const Geometry& geom, MultiFab& scratch)
{
    const IntVect ng = mf.nGrowVect();
    bool any_periodic = false;
    bool any_outside  = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) {
            any_periodic = any_periodic || ng[d] > 0;
        } else {
            any_outside  = any_outside  || ng[d] > 0;
        }
    }
    if (!any_periodic || !any_outside) { return mf; }

    const BoxArray& ba = mf.boxArray();
    IntVect per_ng(0);
    bool nodal_periodic = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) {
            per_ng[d] = ng[d];
            nodal_periodic = nodal_periodic || mf.ixType().nodeCentered(d);
        }
    }

    // The temporary's ghosts start from mf's too: FillBoundary leaves alone a ghost with no
    // valid cell to copy from -- one at a coarse-fine interface of a level that does not span
    // the domain -- and that must come out as it went in, not as uninitialized memory.
    MultiFab ext(extended_boxarray(ba, geom, ng), mf.DistributionMap(), mf.nComp(), per_ng);
    const int ncomp = mf.nComp();
    for (MFIter mfi(ext); mfi.isValid(); ++mfi) {
        const auto& e = ext.array(mfi);
        const auto& a = mf.const_array(mfi);
        ParallelFor(mfi.fabbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            e(i,j,k,n) = a(i,j,k,n);
        });
    }
    if (nodal_periodic) { ext.OverrideSync(geom.periodicity()); }
    ext.FillBoundary(geom.periodicity());

    // Ghosts the temporary does not cover -- those of a box away from the domain faces --
    // come from mf unchanged.
    scratch.define(ba, mf.DistributionMap(), ncomp, ng);
    MultiFab::Copy(scratch, mf, 0, 0, ncomp, ng);
    for (MFIter mfi(ext); mfi.isValid(); ++mfi) {
        const auto& e = ext.const_array(mfi);
        const auto& a = scratch.array(mfi);
        ParallelFor(mfi.fabbox(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            a(i,j,k,n) = e(i,j,k,n);
        });
    }
    return scratch;
}
} // namespace

/**
 * Coarse data ready for a coarse-to-fine interpolation (FillPatchTwoLevels,
 * InterpFromCoarseLevel), which copies the coarse cells, ghost cells included, with
 * periodicity. The ghost cells are exchanged across grids and periodic edges, and those that
 * also lie beyond a non-periodic domain face are made to agree with their periodic images
 * (periodic_consistent_ghosts), so that no cell of the copy has two different sources (#768,
 * #774, #776). Every coarse-to-fine fill goes through here.
 *
 * The ghost cells beyond a wall are otherwise taken as the coarse level's own fill left them;
 * they are not refilled from the boundary conditions here.
 *
 * While crse_lev is frozen (freeze_coarse_level) and reuse is set, the result for a given
 * crse is made once and handed out again: crse does not change in between, so its exchange
 * would change nothing and the copy would come out the same. The copy then lives in the
 * freeze's entry rather than in scratch. Only pass reuse for level data that lives through
 * the freeze, since entries are found by address.
 *
 * @param[in   ] crse_lev  level of crse
 * @param[inout] crse      coarse data; its ghost cells are exchanged in place
 * @param[inout] scratch   storage for a consistent copy, used only when one is needed
 * @param[in   ] reuse     whether the result may be kept and reused while crse_lev is frozen
 * @return crse itself, or the consistent copy
 */
MultiFab&
REMORA::coarse_for_interp (int crse_lev, MultiFab& crse, MultiFab& scratch, bool reuse)
{
    CrseInterpEntry* e = reuse ? crse_interp_entry(crse_lev, crse) : nullptr;
    if (e && e->ready) { return e->copy ? *e->copy : crse; }

    exchange_coarse(crse_lev, crse, reuse);
    MultiFab& dst = e ? *(e->copy = std::make_unique<MultiFab>()) : scratch;
    const MultiFab& src = periodic_consistent_ghosts(crse, geom[crse_lev], dst);
    if (e) {
        e->ready = true;
        if (&src == &crse) { e->copy.reset(); }
    }
    return (&src == &crse) ? crse : dst;
}

/**
 * @param[in   ] crse_lev  level of crse
 * @param[inout] crse      coarse data; its ghost cells are exchanged in place
 * @param[in   ] reuse     whether an exchange already made while crse_lev is frozen counts
 */
void
REMORA::exchange_coarse (int crse_lev, MultiFab& crse, bool reuse)
{
    CrseInterpEntry* e = reuse ? crse_interp_entry(crse_lev, crse) : nullptr;
    if (e && e->exchanged) { return; }
    crse.FillBoundary(geom[crse_lev].periodicity());
    if (e) { e->exchanged = true; }
}

/**
 * timeStep freezes a level while its child substeps: the parent has finished its step, and
 * nothing writes it again until the refluxing and averaging down that follow the substeps.
 * Ending the freeze drops what coarse_for_interp kept.
 *
 * @param[in   ] crse_lev  level to freeze or thaw
 * @param[in   ] on        true to freeze
 */
void
REMORA::freeze_coarse_level (int crse_lev, bool on)
{
    if (int(crse_frozen.size()) <= crse_lev) {
        crse_frozen.resize(crse_lev+1, 0);
        crse_interp_entries.resize(crse_lev+1);
    }
    crse_frozen[crse_lev] = on ? 1 : 0;
    crse_interp_entries[crse_lev].clear();
}

/**
 * @param[in   ] crse_lev  coarse level
 * @param[in   ] crse      coarse data
 * @return the entry for crse if crse_lev is frozen, else nullptr
 */
REMORA::CrseInterpEntry*
REMORA::crse_interp_entry (int crse_lev, const MultiFab& crse)
{
    if (crse_lev >= int(crse_frozen.size()) || !crse_frozen[crse_lev]) { return nullptr; }
    auto& entries = crse_interp_entries[crse_lev];
    for (auto& e : entries) {
        if (e.src == &crse && e.ba == crse.boxArray() && e.dm == crse.DistributionMap()) {
            return &e;
        }
    }
    entries.emplace_back();
    CrseInterpEntry& e = entries.back();
    e.src = &crse;
    e.ba  = crse.boxArray();
    e.dm  = crse.DistributionMap();
    return &e;
}

/**
 * Fill valid and ghost data in the MultiFab "mf"
 * This version fills the MultiFab mf_to_fill in valid regions with
 * the "state data" at the given time;
 * values in mf when it is passed in are *not* used.
 *
 * @param[in   ] lev            level to FillPatch on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[in   ] mfs            vector over levels of multifabs associated with mf_to_fill
 * @param[in   ] bccomp         index into domain_bcs_type
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          component to fill
 * @param[in   ] fill_all       whether to fill all components
 * @param[in   ] fill_set       whether to use FillPatchers
 * @param[in   ] icomp_calc     component to use in RHS boundary calculation
 * @param[in   ] dt_lev         time step at this level
 * @param[in   ] mf_calc        data to use in RHS boundary calculation
 */
void
REMORA::FillPatch (int lev, Real time, MultiFab& mf_to_fill, Vector<MultiFab*> const& mfs,
                  const int  bccomp,
                  const int  bdy_var_type,
                  const int  icomp,
                  const bool fill_all,
                  const bool fill_set,
                  const int  n_not_fill,
                  const int  icomp_calc,
                  const Real dt_lev,
                  const MultiFab& mf_calc,
                  Vector<MultiFab*> const& mfs_crse_old,
                  Vector<MultiFab*> const& mfs_crse_new,
                  const int cf_cache_slot)
{
    BL_PROFILE_VAR("REMORA::FillPatch()",REMORA_FillPatch);
    amrex::Interpolater* mapper = nullptr;

    Box mf_box(mf_to_fill.boxArray()[0]);
    bool is_2d = mf_box.length(2) == 1;

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written below because the FillPatch operators see these as
    // valid faces.
    // ***************************************************************************
    if (lev>0 && fill_set) {
        if (cf_set_width > 0 &&
            mf_box.ixType() == IndexType(IntVect(0,0,0))) {
            FPr_c[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
        } else if (cf_set_width >= 0) {
            if (!is_2d) {
                if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                    FPr_u[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                    FPr_v[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,0,1))) {
                    FPr_w[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                }
            } else {
                if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                    FPr_ubar[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                    FPr_vbar[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                }
            }
        }
    }

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill.nComp();
    } else {
        ncomp = 1;
    }

    if (mf_box.ixType() == IndexType(IntVect(0,0,0)))
    {
        mapper = &cell_cons_interp;
    }
    else {
        mapper = &face_cons_linear_interp;
    }

    if (lev == 0)
    {
        Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        amrex::FillPatchSingleLevel(mf_to_fill, time, fmf, ftime, icomp, icomp, ncomp,
                                    geom[lev], null_bc, bccomp);
    }
    else
    {
        // The barotropic fills (advance_2d) repeat within a parent step for a handful of
        // (component, time) pairs, and what the two-level fill leaves in the coarse-fine
        // ghost cells depends only on the parent's data (fixed until register_coarse_data
        // replaces it), the component and the time -- for faces also on the fine interface
        // faces, which FillSet has just set from the same parent data and time. So keep the
        // result per pair and replay it with the fine-fine exchange that ends the fill,
        // instead of exchanging and interpolating the parent again on every fast step.
        const bool cacheable = cf_2d_cacheable(lev, cf_cache_slot, mf_to_fill, ncomp, fill_set);
        CF2DCache* cache = nullptr;
        long gen = 0;
        bool hit = false;
        if (cacheable) {
            cf_2d_cache_entry(lev, cf_cache_slot, mf_to_fill, *mfs[lev-1], time, mapper, bccomp);
            cache = &cf_2d_cache[lev][cf_cache_slot];
            gen = cf_crse_gen[lev-1];
            hit = (cache->gen[icomp] == gen) && (cache->time[icomp] == time);
        }

        if (hit && cf_2d_fill_cache == 1) {
            copy_on_mask(mf_to_fill, *cache->vals, *cache->mask, icomp);
            mf_to_fill.FillBoundary(icomp, 1, mf_to_fill.nGrowVect(), geom[lev].periodicity());
        } else {
            Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            // A subcycled fine level asks for times inside the parent's step, so the coarse
            // contribution has to be interpolated. Passing one MultiFab twice makes
            // FillPatchTwoLevels return it whatever the time: fine in lockstep, wrong for
            // subcycling.
            MultiFab crse_scratch, crse_old_scratch, crse_new_scratch;
            Vector<MultiFab*> cmf;
            Vector<Real> ctime;
            if (do_substep && int(mfs_crse_old.size()) >= lev && int(mfs_crse_new.size()) >= lev) {
                cmf   = {&coarse_for_interp(lev-1, *mfs_crse_old[lev-1], crse_old_scratch, true),
                         &coarse_for_interp(lev-1, *mfs_crse_new[lev-1], crse_new_scratch, true)};
                ctime = {t_old[lev-1], t_new[lev-1]};
                exchange_coarse(lev-1, *mfs[lev-1], true);
            } else {
                MultiFab& crse = coarse_for_interp(lev-1, *mfs[lev-1], crse_scratch, true);
                cmf   = {&crse, &crse};
                ctime = {time, time};
            }

            amrex::FillPatchTwoLevels(mf_to_fill, mf_to_fill.nGrowVect(), IntVect(0,0,0),
                                      time, cmf, ctime, fmf, ftime,
                                      icomp, icomp, ncomp, geom[lev-1], geom[lev],
                                      //null_bc, bccomp, null_bc, bccomp,
                                      refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);

            if (cacheable && hit) {
                const Long ndiff = count_mask_diff(mf_to_fill, *cache->vals, *cache->mask, icomp);
                if (ndiff > 0) {
                    amrex::Abort("remora.cf_2d_fill_cache = 2: the cached coarse-fine fill differs "
                                 "from a fresh one in " + std::to_string(ndiff) + " cells (level "
                                 + std::to_string(lev) + ", slot " + std::to_string(cf_cache_slot)
                                 + ")");
                }
            } else if (cacheable) {
                copy_on_mask(*cache->vals, mf_to_fill, *cache->mask, icomp);
                cache->gen[icomp]  = gen;
                cache->time[icomp] = time;
            }
        }
    } // lev > 0

    FillPatchBCs(lev, time, mf_to_fill, bccomp, bdy_var_type, icomp, ncomp, n_not_fill,
                 icomp_calc, dt_lev, mf_calc);
}

/**
 * Whether FillPatch may keep and replay the coarse-fine fill of mf in cf_2d_cache. Faces
 * qualify only when FillSet has rewritten the interface faces the fill reads first.
 *
 * @param[in] lev       fine level
 * @param[in] slot      cache slot, or -1 for none
 * @param[in] mf        the MultiFab filled
 * @param[in] ncomp     number of components filled
 * @param[in] fill_set  whether the fill starts with FillSet
 */
bool
REMORA::cf_2d_cacheable (int lev, int slot, const MultiFab& mf, int ncomp, bool fill_set) const
{
    const bool is_face = mf.ixType() != IndexType(IntVect(0,0,0));
    return slot >= 0 && cf_2d_fill_cache > 0 && do_substep
        && ncomp == 1 && int(cf_crse_gen.size()) >= lev
        && (!is_face || (fill_set && cf_width >= 0 && cf_set_width >= 0));
}

/**
 * The barotropic fills of a fine level in advance_2d -- ubar and vbar with FillSet, zeta
 * without, component knew at time -- done together when all three are cache hits. That is
 * what FillPatch does for each in turn, with the three fine-fine exchanges posted together
 * in one FillBoundary, as at level 0. The reordering changes no value: FillSet writes only
 * its own field (the interface faces, in every component), the replay only its field's knew
 * ghost cells, and an exchange only its own field; the boundary conditions of each field read
 * only that field, zeta's know component (Flather and Chapman), h and the masks, and still
 * run in the same order, each after its own field's FillSet, replay and exchange.
 *
 * @param[in] lev        fine level
 * @param[in] time       fill time
 * @param[in] knew       component to fill
 * @param[in] know       component the boundary conditions read
 * @param[in] dt2d       barotropic step the boundary conditions use
 * @param[in] ubar_know  component know of ubar, for the boundary conditions
 * @param[in] vbar_know  component know of vbar
 * @param[in] zeta_know  component know of zeta
 * @return whether the fills were done; if false nothing was changed
 */
bool
REMORA::replay_cf_2d_fills (int lev, Real time, int knew, int know, Real dt2d,
                            const MultiFab& ubar_know, const MultiFab& vbar_know,
                            const MultiFab& zeta_know)
{
    if (cf_2d_fill_cache != 1) { return false; }

    MultiFab* mfs[3]     = {vec_ubar[lev].get(), vec_vbar[lev].get(), vec_zeta[lev].get()};
    MultiFab* crse[3]    = {vec_ubar[lev-1].get(), vec_vbar[lev-1].get(), vec_zeta[lev-1].get()};
    const int bccomp[3]  = {ubar_bc(), vbar_bc(), zeta_bc()};
    const int bdy[3]     = {bdy_ubar(), bdy_vbar(), bdy_zeta()};
    const bool fset[3]   = {true, true, false};
    const MultiFab* mf_calc[3] = {&ubar_know, &vbar_know, &zeta_know};

    CF2DCache* cache[3];
    for (int n = 0; n < 3; ++n) {
        if (!cf_2d_cacheable(lev, n, *mfs[n], 1, fset[n])) { return false; }
        Interpolater* mapper = (mfs[n]->ixType() == IndexType(IntVect(0,0,0)))
                             ? static_cast<Interpolater*>(&cell_cons_interp)
                             : static_cast<Interpolater*>(&face_cons_linear_interp);
        cf_2d_cache_entry(lev, n, *mfs[n], *crse[n], time, mapper, bccomp[n]);
        cache[n] = &cf_2d_cache[lev][n];
        if (cache[n]->gen[knew] != cf_crse_gen[lev-1] || cache[n]->time[knew] != time) {
            return false;
        }
    }

    BL_PROFILE("REMORA::replay_cf_2d_fills()");
    if (cf_set_width >= 0) {
        FPr_ubar[lev-1].FillSet(*mfs[0], time, null_bc, domain_bcs_type);
        FPr_vbar[lev-1].FillSet(*mfs[1], time, null_bc, domain_bcs_type);
    }
    for (int n = 0; n < 3; ++n) {
        copy_on_mask(*mfs[n], *cache[n]->vals, *cache[n]->mask, knew);
    }
    const Periodicity& per = geom[lev].periodicity();
    amrex::FillBoundary(Vector<MultiFab*>{mfs[0], mfs[1], mfs[2]},
                        Vector<int>{knew, knew, knew}, Vector<int>{1, 1, 1},
                        Vector<IntVect>{mfs[0]->nGrowVect(), mfs[1]->nGrowVect(),
                                        mfs[2]->nGrowVect()},
                        Vector<Periodicity>{per, per, per});
    for (int n = 0; n < 3; ++n) {
        FillPatchBCs(lev, time, *mfs[n], bccomp[n], bdy[n], knew, 1, 0, know, dt2d, *mf_calc[n]);
    }
    return true;
}

/**
 * Fit cf_2d_cache[lev][slot] to mf's layout. When the layout changes, the cached values are
 * dropped and the mask is rebuilt: the ghost cells a two-level fill of mf writes, which is
 * the coarse-fine ghost region together with whatever the interpolater's patches overlap.
 * Rather than restate AMReX's region logic, run the same fill on a probe whose ghost cells
 * hold a sentinel no fill can produce (valid cells 0, parent 1) and mark what changed.
 *
 * @param[in] lev     fine level
 * @param[in] slot    cache slot (0 ubar, 1 vbar, 2 zeta)
 * @param[in] mf      the MultiFab the fills write
 * @param[in] crse    the parent level's MultiFab of the same field, for its layout
 * @param[in] time    a fill time (only the layout matters)
 * @param[in] mapper  the interpolater the fills use
 * @param[in] bccomp  index into domain_bcs_type the fills use
 */
void
REMORA::cf_2d_cache_entry (int lev, int slot, const MultiFab& mf, const MultiFab& crse,
                           Real time, Interpolater* mapper, int bccomp)
{
    if (int(cf_2d_cache.size()) <= lev) { cf_2d_cache.resize(lev+1); }
    CF2DCache& c = cf_2d_cache[lev][slot];
    if (c.vals && c.vals->boxArray() == mf.boxArray()
               && c.vals->DistributionMap() == mf.DistributionMap()
               && c.vals->nComp() == mf.nComp() && c.vals->nGrowVect() == mf.nGrowVect()) {
        return;
    }

    c.vals = std::make_unique<MultiFab>(mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                                        mf.nGrowVect());
    c.mask = std::make_unique<iMultiFab>(mf.boxArray(), mf.DistributionMap(), 1,
                                         mf.nGrowVect());
    c.time.assign(mf.nComp(), zero);
    c.gen.assign(mf.nComp(), -1);

    constexpr Real sentinel = Real(-1.0e30);
    MultiFab probe(mf.boxArray(), mf.DistributionMap(), 1, mf.nGrowVect());
    probe.setVal(sentinel);
    probe.setVal(zero, 0, 1, 0);
    MultiFab cprobe(crse.boxArray(), crse.DistributionMap(), 1, crse.nGrowVect());
    cprobe.setVal(one);

    Vector<MultiFab*> fmf = {&probe, &probe};
    Vector<Real> ftime    = {time, time};
    Vector<MultiFab*> cmf = {&cprobe, &cprobe};
    Vector<Real> ctime    = {time, time};
    amrex::FillPatchTwoLevels(probe, probe.nGrowVect(), IntVect(0,0,0),
                              time, cmf, ctime, fmf, ftime,
                              0, 0, 1, geom[lev-1], geom[lev],
                              refRatio(lev-1),
                              mapper, domain_bcs_type, bccomp);

    // Ghost cells another fine box covers are left out: the exchange that ends every fill
    // overwrites them with that box's current values, whatever the interpolation put there.
    iMultiFab covered(mf.boxArray(), mf.DistributionMap(), 1, mf.nGrowVect());
    covered.setVal(0);
    covered.setVal(1, 0, 1, 0);
    covered.FillBoundary(geom[lev].periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*c.mask, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto m = c.mask->array(mfi);
        auto cv = covered.const_array(mfi);
        auto pr = probe.const_array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            m(i,j,k) = (cv(i,j,k) == 0 && pr(i,j,k) != sentinel) ? 1 : 0;
        });
    }
}

/**
 * The boundary-condition half of FillPatch, for a MultiFab whose fine-fine (and, above
 * level 0, coarse-fine) ghost cells are already filled: physical boundary conditions,
 * boundary data from NetCDF files, periodic corners, and free-slip at the top for the
 * face-centered velocities. Split out so that several fields can share one exchange and
 * then take their conditions in turn (advance_2d does this for ubar, vbar and zeta).
 *
 * @param[in   ] lev            level to fill on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[in   ] bccomp         index into domain_bcs_type
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          first component to fill
 * @param[in   ] ncomp          number of components to fill
 * @param[in   ] n_not_fill     number of ghost cells not to fill
 * @param[in   ] icomp_calc     component to use in RHS boundary calculation
 * @param[in   ] dt_lev         time step at this level
 * @param[in   ] mf_calc        data to use in RHS boundary calculation
 */
void
REMORA::FillPatchBCs (int lev, Real time, MultiFab& mf_to_fill,
                      const int  bccomp,
                      const int  bdy_var_type,
                      const int  icomp,
                      const int  ncomp,
                      const int  n_not_fill,
                      const int  icomp_calc,
                      const Real dt_lev,
                      const MultiFab& mf_calc)
{
    Box mf_box(mf_to_fill.boxArray()[0]);

    MultiFab* mask;
    if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
        mask = vec_msku[lev].get();
    } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
        mask = vec_mskv[lev].get();
    } else {
        mask = vec_mskr[lev].get();
    }

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************

    // Enforce physical boundary conditions
    (*physbcs[lev])(mf_to_fill,*mask,icomp,ncomp,mf_to_fill.nGrowVect(),time,bccomp,n_not_fill,mf_calc,
            *vec_msku[lev].get(), *vec_mskv[lev].get());

#ifdef REMORA_USE_NETCDF
    // Fill the data which is stored in the boundary data read from netcdf files
    if ( (solverChoice.boundary_from_netcdf) && (bdy_var_type != BdyVars::null) )
    {
        fill_from_bdyfiles(lev, mf_to_fill,*mask,time,bccomp,bdy_var_type, icomp,icomp_calc,mf_calc,dt_lev);
    }
#else
    amrex::ignore_unused(bdy_var_type, icomp_calc, dt_lev);
#endif
    // Fill corners of the domain with periodic data
    if  ( mf_box.ixType() == IndexType(IntVect(0,0,0)) ) {
        mf_to_fill.EnforcePeriodicity(geom[lev].periodicity());
    }

    // Also enforce free-slip at top boundary (on xvel or yvel). Without ghost cells in z
    // (ubar and vbar) there is nothing above the top to set, so skip the loop and the
    // device synchronization that comes with it.
    if ( ((mf_box.ixType() == IndexType(IntVect(1,0,0))) ||
          (mf_box.ixType() == IndexType(IntVect(0,1,0)))) && mf_to_fill.nGrowVect()[2] > 0 )
    {
        int khi = geom[lev].Domain().bigEnd(2);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf_to_fill); mfi.isValid(); ++mfi)
        {
            Box gbx  = mfi.growntilebox(); // Note this is face-centered since vel is
            gbx.setSmall(2,khi+1);
            if (gbx.ok()) {
                Array4<Real> vel_arr = mf_to_fill.array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    vel_arr(i,j,k) = vel_arr(i,j,khi);
                });
            }
        }
    }
}

/**
 * Fill valid and ghost data in the MultiFab "mf"
 * This version fills the MultiFab mf_to_fill in valid regions with
 * the "state data" at the given time;
 * values in mf when it is passed in are *not* used.
 * Unlike FillPatch, FillPatchNoBC does not apply boundary conditions
 *
 * @param[in   ] lev            level to FillPatch on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[in   ] mfs            vector over levels of multifabs associated with mf_to_fill
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          component to fill
 * @param[in   ] fill_all       whether to fill all components
 * @param[in   ] fill_set       whether to use FillPatchers
 */
void
REMORA::FillPatchNoBC (int lev, Real time, MultiFab& mf_to_fill, Vector<MultiFab*> const& mfs,
#ifdef REMORA_USE_NETCDF
                  const int /*bdy_var_type*/,
#else
                  const int /*bdy_var_type*/,
#endif
                  const int  icomp,
                  const bool fill_all,
                  const bool fill_set,
                  Vector<MultiFab*> const& mfs_crse_old,
                  Vector<MultiFab*> const& mfs_crse_new)
{
    // HACK: Note that this is hacky; should be able to have a single call to FillPatch with a
    // flag for bcs, but for some reason it was acting weird, so we're splitting this out into
    // two functions with repeated code for the time being
    BL_PROFILE_VAR("REMORA::FillPatchNoBC()",REMORA_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    Box mf_box(mf_to_fill.boxArray()[0]);
    bool is_2d = mf_box.length(2) == 1;

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written below because the FillPatch operators see these as
    // valid faces.
    // ***************************************************************************
    if (lev>0 && fill_set) {
        if (cf_set_width > 0 &&
            mf_box.ixType() == IndexType(IntVect(0,0,0))) {
            FPr_c[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
        } else if (cf_set_width >= 0) {
            if (!is_2d) {
                if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                    FPr_u[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                    FPr_v[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,0,1))) {
                    FPr_w[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                }
            } else {
                if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                    FPr_ubar[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                    FPr_vbar[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
                }
            }
        }
    }

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill.nComp();
    } else {
        ncomp = 1;
    }

    if (mf_box.ixType() == IndexType(IntVect(0,0,0)))
    {
        bccomp = 0;
        mapper = &cell_cons_interp;
    }
    else if (mf_box.ixType() == IndexType(IntVect(1,0,0)))
    {
        bccomp = xvel_bc();
        mapper = &face_cons_linear_interp;
    }
    else if (mf_box.ixType() == IndexType(IntVect(0,1,0)))
    {
        bccomp = yvel_bc();
        mapper = &face_cons_linear_interp;
    }
    else {
        bccomp = zvel_bc();
        mapper = &face_cons_linear_interp;
    }

    if (lev == 0)
    {
        Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        amrex::FillPatchSingleLevel(mf_to_fill, time, fmf, ftime, icomp, icomp, ncomp,
                                    geom[lev], null_bc, bccomp);
    }
    else
    {
        Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        // A subcycled fine level asks for times inside the parent's step, so the coarse
        // contribution has to be interpolated. Passing one MultiFab twice makes
        // FillPatchTwoLevels return it whatever the time: fine in lockstep, wrong for
        // subcycling.
        MultiFab crse_scratch, crse_old_scratch, crse_new_scratch;
        Vector<MultiFab*> cmf;
        Vector<Real> ctime;
        if (do_substep && int(mfs_crse_old.size()) >= lev && int(mfs_crse_new.size()) >= lev) {
            cmf   = {&coarse_for_interp(lev-1, *mfs_crse_old[lev-1], crse_old_scratch, true),
                     &coarse_for_interp(lev-1, *mfs_crse_new[lev-1], crse_new_scratch, true)};
            ctime = {t_old[lev-1], t_new[lev-1]};
            exchange_coarse(lev-1, *mfs[lev-1], true);
        } else {
            MultiFab& crse = coarse_for_interp(lev-1, *mfs[lev-1], crse_scratch, true);
            cmf   = {&crse, &crse};
            ctime = {time, time};
        }

        amrex::FillPatchTwoLevels(mf_to_fill, mf_to_fill.nGrowVect(), IntVect(0,0,0),
                                  time, cmf, ctime, fmf, ftime,
                                  icomp, icomp, ncomp, geom[lev-1], geom[lev],
                                  //null_bc, bccomp, null_bc, bccomp,
                                  refRatio(lev-1),
                                  mapper, domain_bcs_type, bccomp);
    } // lev > 0
}


// utility to copy in data from old/new data into a struct that holds data for FillPatching
TimeInterpolatedData
REMORA::GetDataAtTime (int /*lev*/, Real /*time*/)
{
    BL_PROFILE_VAR("GetDataAtTime()",GetDataAtTime);
    TimeInterpolatedData data;

// HACK HACK HACK
#if 0
    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3_rt;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_new[lev][i], data.non_owning);
        }
        data.set_time(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_old[lev][i], data.non_owning);
        }
        data.set_time(t_old[lev]);
    }
    else if (time > t_old[lev] && time < t_new[lev])
    {
        // do first order interpolation in time between [t_old[lev], t_new[lev]]
        // time interpolation includes the ghost cells
        for (int i = 0; i < Vars::NumTypes; ++i) {
            MultiFab* mf_tmp = new MultiFab(vars_new[lev][i].boxArray(),
                                            vars_new[lev][i].DistributionMap(),
                                            vars_new[lev][i].nComp(), vars_new[lev][i].nGrowVect());
            mf_tmp->setVal(0.0_rt);

            const Real dt_fraction = (time - t_old[lev]) / (t_new[lev] - t_old[lev]);
            MultiFab::Saxpy(*mf_tmp, 1.0_rt - dt_fraction, vars_old[lev][i], 0, 0, mf_tmp->nComp(), mf_tmp->nGrowVect());
            MultiFab::Saxpy(*mf_tmp,          dt_fraction, vars_new[lev][i], 0, 0, mf_tmp->nComp(), mf_tmp->nGrowVect());

            data.add_var(mf_tmp, data.owning);
        }
        data.set_time(time);
    }
    else
    {
        amrex::Error("Requested data at a time outside the interval [t_old, t_new]");
    }

    // We need to make sure to fill these before we compute the viscosity
    for (int i = 0; i < Vars::NumTypes; ++i) {
        data.get_var(Vars::xvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::yvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::zvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::cons).FillBoundary(geom[lev].periodicity());
    }
#endif
    return data;
}


/**
 * Like FillCoarsePatch but uses piecewise constant interpolater
 *
 * @param[in   ] lev            level to FillPatch on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[inout] mf_crse        MultiFab of coarse data
 * @param[in   ] bccomp         index into domain_bcs_type
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          component to fill
 * @param[in   ] fill_all       whether to fill all components
 * @param[in   ] icomp_calc     component to use in RHS boundary calculation
 * @param[in   ] dt_lev         time step at this level
 * @param[in   ] mf_calc        data to use in RHS boundary calculation
 */
void
REMORA::FillCoarsePatchPC (int lev, Real time, MultiFab* mf_to_fill,
                           MultiFab* mf_crse,
                           const int bccomp,
                           const int bdy_var_type,
                           const int  icomp,
                           const bool fill_all,
                           const int  n_not_fill,
                           const int  icomp_calc,
                           const Real dt_lev,
                           const MultiFab& mf_calc) {
    amrex::Interpolater* mapper = &pc_interp;
    FillCoarsePatchMap(lev, time, mf_to_fill, mf_crse, bccomp, bdy_var_type,icomp,
                    fill_all, n_not_fill, icomp_calc, dt_lev, mf_calc, mapper);
}

/**
 * Fill an entire multifab by interpolating from the coarser level --
 * this is used only when a new level of refinement is being created
 * during a run (i.e not at initialization)
 * This will never be used with static refinement.
 *
 * @param[in   ] lev            level to FillPatch on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[inout] mf_crse        MultiFab of coarse data
 * @param[in   ] bccomp         index into domain_bcs_type
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          component to fill
 * @param[in   ] fill_all       whether to fill all components
 * @param[in   ] icomp_calc     component to use in RHS boundary calculation
 * @param[in   ] dt_lev         time step at this level
 * @param[in   ] mf_calc        data to use in RHS boundary calculation
 */
void
REMORA::FillCoarsePatch (int lev, Real time, MultiFab* mf_to_fill,
                         MultiFab* mf_crse,
                         const int bccomp,
                         const int bdy_var_type,
                         const int  icomp,
                         const bool fill_all,
                         const int  n_not_fill,
                         const int  icomp_calc,
                         const Real dt_lev,
                         const MultiFab& mf_calc) {
    amrex::Interpolater* mapper = nullptr;
    FillCoarsePatchMap(lev, time, mf_to_fill, mf_crse, bccomp, bdy_var_type,icomp,
                    fill_all, n_not_fill, icomp_calc, dt_lev, mf_calc, mapper);
}

/**
 * Fill an entire multifab by interpolating from the coarser level --
 * this is used only when a new level of refinement is being created
 * during a run (i.e not at initialization)
 * This will never be used with static refinement.
 *
 * @param[in   ] lev            level to FillPatch on
 * @param[in   ] time           current time
 * @param[inout] mf_to_fill     MultiFab to fill
 * @param[inout] mf_crse        MultiFab of coarse data
 * @param[in   ] bccomp         index into domain_bcs_type
 * @param[in   ] bdy_var_type   when filling from NetCDF, which boundary data
 * @param[in   ] icomp          component to fill
 * @param[in   ] fill_all       whether to fill all components
 * @param[in   ] icomp_calc     component to use in RHS boundary calculation
 * @param[in   ] dt_lev         time step at this level
 * @param[in   ] mf_calc        data to use in RHS boundary calculation
 */
void
REMORA::FillCoarsePatchMap (int lev, Real time, MultiFab* mf_to_fill, MultiFab* mf_crse,
                         const int bccomp,
#ifdef REMORA_USE_NETCDF
                         const int bdy_var_type,
#else
                         const int /*bdy_var_type*/,
#endif
                         const int  icomp,
                         const bool fill_all,
                         const int  n_not_fill,
                         const int  icomp_calc,
                         const Real dt_lev,
                         const MultiFab& mf_calc,
                               amrex::Interpolater* mapper)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    MultiFab* mask = nullptr;

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill->nComp();
    } else {
        ncomp = 1;
    }

    Box box_mf(mf_to_fill->boxArray()[0]);

    if (box_mf.ixType() == IndexType(IntVect(0,0,0)))
    {
        mask = vec_mskr[lev].get();
    }
    else if (box_mf.ixType() == IndexType(IntVect(1,0,0)))
    {
        mask = vec_msku[lev].get();
    }
    else if (box_mf.ixType() == IndexType(IntVect(0,1,0)))
    {
        mask = vec_mskv[lev].get();
    }
    else if (box_mf.ixType() == IndexType(IntVect(0,0,1)))
    {
        mask = vec_mskr[lev].get();
    } else {
          amrex::Abort("Dont recognize this box type in REMORA_FillPatch");
    }

    if (mapper == nullptr) {
        if (box_mf.ixType() == IndexType(IntVect(0,0,0)))
        {
            mapper = &cell_cons_interp;
        }
        else if (box_mf.ixType() == IndexType(IntVect(1,0,0)))
        {
            mapper = &face_cons_linear_interp;
        }
        else if (box_mf.ixType() == IndexType(IntVect(0,1,0)))
        {
            mapper = &face_cons_linear_interp;
        }
        else if (box_mf.ixType() == IndexType(IntVect(0,0,1)))
        {
            mapper = &face_cons_linear_interp;
        } else {
              amrex::Abort("Dont recognize this box type in REMORA_FillPatch");
        }
    }

    // InterpFromCoarseLevel copies the coarse ghost cells with periodicity. When the levels
    // are first made the parent has had no FillPatch yet, so its ghosts beyond a wall and
    // beyond a periodic edge still hold initial values rather than their periodic images.
    // Where the new patch touches the seam the two then compete for one cell of the copy,
    // and on a GPU the winner varies from run to run (Channel_Test_ML_seam, issue #768).
    MultiFab crse_scratch;
    const MultiFab& crse_src = coarse_for_interp(lev-1, *mf_crse, crse_scratch);
    amrex::InterpFromCoarseLevel(*mf_to_fill, mf_to_fill->nGrowVect(), IntVect(0,0,0),
            crse_src, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                 refRatio(lev-1),
                                 mapper, domain_bcs_type, bccomp);

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************

    // The radiation conditions compare the field with its previous state (mf_calc),
    // which a level being created does not have. Give them the interpolated field
    // itself: the ghost cells then keep the coarse values and the nudging term,
    // with no time step, is zero. Only the radiation conditions read it, so the
    // copy is made only when one of them is set on this field, and only of the
    // components being filled. It cannot be an alias of the field: the radiation
    // kernels read the along-boundary neighbours of the column they write, so on a
    // GPU a thread could read a cell another thread has already overwritten.
    MultiFab calc_copy;
    const MultiFab* calc = &mf_calc;
    int icomp_calc_used = icomp_calc;
    if (!mf_calc.isDefined()) {
        bool radiation = false;
        for (int n = 0; n < ncomp; ++n) {
            const BCRec& bcr = domain_bcs_type[bccomp+n];
            for (int dir = 0; dir < 2; ++dir) {
                for (int bc : {bcr.lo(dir), bcr.hi(dir)}) {
                    radiation = radiation || bc == REMORABCType::orlanski_rad
                                          || bc == REMORABCType::orlanski_rad_nudge;
                }
            }
        }
        if (radiation) {
            calc_copy.define(mf_to_fill->boxArray(), mf_to_fill->DistributionMap(),
                             ncomp, mf_to_fill->nGrowVect());
            MultiFab::Copy(calc_copy, *mf_to_fill, icomp, 0, ncomp,
                           mf_to_fill->nGrowVect());
            calc = &calc_copy;
            icomp_calc_used = 0;
        }
    }

    // Enforce physical boundary conditions
    (*physbcs[lev])(*mf_to_fill,*mask,icomp,ncomp,mf_to_fill->nGrowVect(),
            time,bccomp,n_not_fill,*calc,
            *vec_msku[lev].get(), *vec_mskv[lev].get());

#ifdef REMORA_USE_NETCDF
    // Fill the data which is stored in the boundary data read from netcdf files
    if ( (solverChoice.boundary_from_netcdf) && (bdy_var_type != BdyVars::null) )
    {
        fill_from_bdyfiles(lev, *mf_to_fill,*mask,time,bccomp,bdy_var_type, icomp,icomp_calc_used,*calc,dt_lev);
    }
#else
    amrex::ignore_unused(icomp_calc_used, dt_lev);
#endif
    // Fill corners of the domain with periodic data
    if  ( box_mf.ixType() == IndexType(IntVect(0,0,0)) ) {
        mf_to_fill->EnforcePeriodicity(geom[lev].periodicity());
    } else {
        // The face-centred analogue. InterpFromCoarseLevel fills a periodic ghost from the
        // coarse level, which for a face tangential to the periodic direction does not
        // reproduce the fine value it wraps to: on a new level of Channel_Test, vbar's x
        // ghosts at the seam differed from their images by 2e-4 to 6e-4 where zeta's,
        // corrected above, were exact. The first barotropic substep then sees different vbar
        // at the two copies of the periodic u-face, and that one-shot asymmetry amplifies into
        // a blow-up. Copying the ghosts from this level's own valid data, as every later
        // FillPatch does, removes it.
        mf_to_fill->FillBoundary(geom[lev].periodicity());
    }

    // Enforce free-slip at top boundary (on xvel or yvel)
    if ( (box_mf.ixType() == IndexType(IntVect(1,0,0))) ||
         (box_mf.ixType() == IndexType(IntVect(0,1,0))) )
    {
        int khi = geom[lev].Domain().bigEnd(2);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*mf_to_fill); mfi.isValid(); ++mfi)
        {
            Box gbx  = mfi.growntilebox(); // Note this is face-centered since vel is
            gbx.setSmall(2,khi+1);
            if (gbx.ok()) {
                Array4<Real> vel_arr = mf_to_fill->array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    vel_arr(i,j,k) = vel_arr(i,j,khi);
                });
            }
        }
    }

}

namespace {
/**
 * Sign for reflecting one cell-centered velocity component into a ghost cell.
 * The normal velocity vanishes at a wall or symmetry plane; only a no-slip wall
 * zeroes the tangential ones too. Anything else copies the adjacent valid cell.
 *
 * @param[in   ] bc        physical BC type for that component on that face
 * @param[in   ] is_normal whether the component is normal to the face
 */
Real
cc_vel_bdy_mult (REMORA_BC bc, bool is_normal)
{
    if (is_normal) {
        return ( (bc == REMORA_BC::slip_wall   ) ||
                 (bc == REMORA_BC::no_slip_wall) ||
                 (bc == REMORA_BC::symmetry    ) ) ? -one : one;
    }
    return (bc == REMORA_BC::no_slip_wall) ? -one : one;
}
} // namespace

/**
 * Fill boundary ghost cells for cell-centered velocities. Each component takes
 * the sign implied by its own BC, which with remora.boundary_per_variable can
 * differ between u, v and w on one face.
 *
 * @param[in   ] lev        level to fill at
 * @param[inout] mf_cc_vel  MultiFab of cell-centered velocities
 */
void
REMORA::FillBdyCCVels (int lev, MultiFab& mf_cc_vel)
{
    // Impose bc's at domain boundaries
    Box domain(Geom(lev).Domain());

    int ihi = domain.bigEnd(0);
    int jhi = domain.bigEnd(1);
    int khi = domain.bigEnd(2);

    // Impose periodicity first
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    // Sign multipliers per face, per velocity component
    GpuArray<GpuArray<Real,AMREX_SPACEDIM>,AMREX_SPACEDIM*2> mult;
    for (OrientationIter oit; oit; ++oit) {
        Orientation ori = oit();
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            mult[ori][n] = cc_vel_bdy_mult(phys_bc_type[xvel_bc()+n][ori],
                                           n == ori.coordDir());
        }
    }

    // Callers pass different ghost widths (the plotfile ones carry none in z)
    const IntVect ng = mf_cc_vel.nGrowVect();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_cc_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Note that we don't fill corners here -- only the cells that share a face
        //      with interior cells -- this is all that is needed to calculate vorticity
        const Box& bx = mfi.tilebox();
        const Array4<Real>& vel_arr = mf_cc_vel.array(mfi);

        if (!Geom(lev).isPeriodic(0) && ng[0] > 0) {
            // Low-x side
            if (bx.smallEnd(0) <= domain.smallEnd(0)) {
                auto const m = mult[Orientation(Direction::x,Orientation::low)];
                ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                {
                    vel_arr(-1,j,k,0) = m[0]*vel_arr(0,j,k,0); // u
                    vel_arr(-1,j,k,1) = m[1]*vel_arr(0,j,k,1); // v
                    vel_arr(-1,j,k,2) = m[2]*vel_arr(0,j,k,2); // w
                });
            }

            // High-x side
            if (bx.bigEnd(0) >= domain.bigEnd(0)) {
                auto const m = mult[Orientation(Direction::x,Orientation::high)];
                ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                {
                    vel_arr(ihi+1,j,k,0) = m[0]*vel_arr(ihi,j,k,0); // u
                    vel_arr(ihi+1,j,k,1) = m[1]*vel_arr(ihi,j,k,1); // v
                    vel_arr(ihi+1,j,k,2) = m[2]*vel_arr(ihi,j,k,2); // w
                });
            }
        } // !periodic

        if (!Geom(lev).isPeriodic(1) && ng[1] > 0) {
            // Low-y side
            if (bx.smallEnd(1) <= domain.smallEnd(1)) {
                auto const m = mult[Orientation(Direction::y,Orientation::low)];
                ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int  , int k) noexcept
                {
                    vel_arr(i,-1,k,0) = m[0]*vel_arr(i,0,k,0); // u
                    vel_arr(i,-1,k,1) = m[1]*vel_arr(i,0,k,1); // v
                    vel_arr(i,-1,k,2) = m[2]*vel_arr(i,0,k,2); // w
                });
            }

            // High-y side
            if (bx.bigEnd(1) >= domain.bigEnd(1)) {
                auto const m = mult[Orientation(Direction::y,Orientation::high)];
                ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int , int k) noexcept
                {
                    vel_arr(i,jhi+1,k,0) = m[0]*vel_arr(i,jhi,k,0); // u
                    vel_arr(i,jhi+1,k,1) = m[1]*vel_arr(i,jhi,k,1); // v
                    vel_arr(i,jhi+1,k,2) = m[2]*vel_arr(i,jhi,k,2); // w
                });
            }
        } // !periodic

        if (!Geom(lev).isPeriodic(2) && ng[2] > 0) {
            // Low-z side
            if (bx.smallEnd(2) <= domain.smallEnd(2)) {
                auto const m = mult[Orientation(Direction::z,Orientation::low)];
                ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    vel_arr(i,j,-1,0) = m[0]*vel_arr(i,j,0,0); // u
                    vel_arr(i,j,-1,1) = m[1]*vel_arr(i,j,0,1); // v
                    vel_arr(i,j,-1,2) = m[2]*vel_arr(i,j,0,2); // w
                });
            }

            // High-z side
            if (bx.bigEnd(2) >= domain.bigEnd(2)) {
                auto const m = mult[Orientation(Direction::z,Orientation::high)];
                ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    vel_arr(i,j,khi+1,0) = m[0]*vel_arr(i,j,khi,0); // u
                    vel_arr(i,j,khi+1,1) = m[1]*vel_arr(i,j,khi,1); // v
                    vel_arr(i,j,khi+1,2) = m[2]*vel_arr(i,j,khi,2); // w
                });
            }
        } // !periodic
    } // MFIter
}

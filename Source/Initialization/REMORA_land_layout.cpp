/**
 * \file REMORA_land_layout.cpp
 *
 * How a level's boxes are laid out over land: the wet-cell weights the distribution mapping
 * can balance (remora.load_balance) and the all-land boxes level 0 can leave out
 * (remora.drop_dry_boxes). Both need the land mask of a layout before the level on it exists,
 * which probe_rho_mask provides.
 */

#include <REMORA.H>
#include <REMORA_prob_common.H>

#include <AMReX_ParallelReduce.H>

using namespace amrex;

namespace {
BoxArray slab_boxes (const BoxArray& ba)
{
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    return BoxArray(std::move(bl2d));
}
}

/**
 * Fill mskr, a rho-point mask on any layout of level lev, with the mask set_masks would give
 * a level on that layout, before such a level exists. On level 0 every ghost cell the level's
 * own mask carries is set, as set_masks sets it; on a finer level only the valid cells are,
 * from the parent's mask (piecewise constant, as set_masks injects it), which is all a weight
 * needs.
 *
 * @param[in   ] lev   level the layout belongs to
 * @param[inout] mskr  mask to fill, on 2D slab boxes, initially anything
 * @return false when the mask cannot be known before the level is made: level 0 with
 *         remora.hires_grid_level >= 0 at initialization, whose mask is coarsened from
 *         full-domain data that init_only builds
 */
bool
REMORA::probe_rho_mask (int lev, MultiFab& mskr)
{
    mskr.setVal(one);
    if (solverChoice.mask_type == MaskType::none) {
        return true;
    }

    if (lev <= hires_grid_level) {
        if (lev >= vec_mskr_full_domain.size() || !vec_mskr_full_domain[lev]) {
            return false;
        }
        // As set_masks_averaged_down
        mskr.ParallelCopy(*vec_mskr_full_domain[lev], 0, 0, 1,
                          vec_mskr_full_domain[lev]->nGrowVect(), mskr.nGrowVect());
        mskr.FillBoundary(geom[lev].periodicity());
        return true;
    }

    if (lev == 0) {
        if (solverChoice.mask_type == MaskType::analytic) {
            prob->init_analytic_masks(lev, geom[lev], solverChoice, *this, mskr);
            mskr.FillBoundary(geom[lev].periodicity());
            return true;
        }
#ifdef REMORA_USE_NETCDF
        fill_rho_mask_from_netcdf(lev, mskr);
        return true;
#else
        return false;
#endif
    }

    // A finer level's valid cells take the parent's mask, piecewise constant
    AMREX_ALWAYS_ASSERT(lev-1 < vec_mskr.size() && vec_mskr[lev-1]);
    const IntVect rr = refRatio(lev-1);
    MultiFab ctmp(amrex::coarsen(mskr.boxArray(), rr), mskr.DistributionMap(), 1, 0);
    ctmp.setVal(one);
    ctmp.ParallelCopy(*vec_mskr[lev-1], 0, 0, 1);
    const int rx = rr[0], ry = rr[1];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mskr, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<      Real> const& f = mskr.array(mfi);
        Array4<const Real> const& c = ctmp.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            f(i,j,k) = c(amrex::coarsen(i, rx), amrex::coarsen(j, ry), 0);
        });
    }
    Gpu::streamSynchronize();
    return true;
}

/**
 * Wet rho cells in each box of a prospective layout of level lev: the box's wet columns times
 * its vertical extent. The same on every rank.
 *
 * @param[in   ] lev  level the layout belongs to
 * @param[in   ] ba   the layout's (3D) boxes
 * @param[  out] wet  wet cells per box of ba
 * @return false when the mask cannot be known yet (see probe_rho_mask)
 */
bool
REMORA::wet_cells_per_box (int lev, const BoxArray& ba, Vector<Long>& wet)
{
    const BoxArray ba2d = slab_boxes(ba);
    MultiFab mskr(ba2d, DistributionMapping(ba2d), 1, IntVect(NGROW+1,NGROW+1,0));
    if (!probe_rho_mask(lev, mskr)) { return false; }

    wet.assign(ba.size(), 0);
    for (MFIter mfi(mskr); mfi.isValid(); ++mfi) {
        // The mask is exactly 0 or 1, so the sum counts wet columns exactly
        const Real ncol = mskr[mfi].sum<RunOn::Device>(mfi.validbox(), 0);
        wet[mfi.index()] = static_cast<Long>(ncol) * ba[mfi.index()].length(2);
    }
    ParallelAllReduce::Sum(wet.data(), static_cast<int>(wet.size()),
                           ParallelContext::CommunicatorSub());
    return true;
}

/**
 * The distribution mapping of a new BoxArray on level lev (overrides AmrMesh). Called by AMReX
 * for every new layout -- initialization, regrid, a level made from coarse -- and by the
 * checkpoint reader. remora.load_balance picks the strategy; each box weighs its wet cells
 * plus remora.load_balance_land_weight times its land cells. Where every box weighs its cell
 * count (no land, or a land weight of one) the SFC strategy is AMReX's own mapping, exactly as
 * remora.load_balance = amrex gives it. Ownership never changes the answer.
 *
 * @param[in   ] lev  level of the BoxArray
 * @param[in   ] ba   the BoxArray to distribute
 */
DistributionMapping
REMORA::MakeDistributionMap (int lev, BoxArray const& ba)
{
    BL_PROFILE("REMORA::MakeDistributionMap()");
    if (load_balance == LoadBalance::amrex) {
        return AmrCore::MakeDistributionMap(lev, ba);
    }

    // A land cell weighing as much as a wet one makes every weight the cell count, so the mask
    // is not needed (nor read from a grid file again).
    Vector<Long> wet;
    const bool known = (load_balance_land_weight != one) && wet_cells_per_box(lev, ba, wet);

    // Weights in hundredths of a wet cell, so a fractional land weight keeps its value
    const Long land_w = std::lround(Real(100.0) * load_balance_land_weight);
    std::vector<Long> wgts(ba.size());
    Long total_cells = 0, total_wet = 0;
    bool cell_count_weights = true;
    for (int b = 0; b < ba.size(); ++b) {
        const Long npts = ba[b].numPts();
        const Long nwet = known ? wet[b] : npts;
        wgts[b] = 100 * nwet + land_w * (npts - nwet);
        cell_count_weights = cell_count_weights && (wgts[b] == 100 * npts);
        total_cells += npts;
        total_wet += nwet;
    }

    if (!known && load_balance_land_weight != one) {
        amrex::Print() << "Level " << lev << ": the land mask is not known before the level is"
                       << " made (remora.hires_grid_level >= " << lev << "), so its boxes are"
                       << " weighted by cell count for remora.load_balance.\n";
    }

    if (cell_count_weights && load_balance == LoadBalance::sfc) {
        return AmrCore::MakeDistributionMap(lev, ba);
    }

    const int nprocs = ParallelContext::NProcsSub();
    DistributionMapping dm;
    if (load_balance == LoadBalance::sfc) {
        dm.SFCProcessorMap(ba, wgts, nprocs);
    } else if (load_balance == LoadBalance::knapsack) {
        dm.KnapSackProcessorMap(wgts, nprocs);
    } else {
        dm.RoundRobinProcessorMap(wgts, nprocs);
    }

    // How even the result is, by the weights balanced and by cell count
    std::vector<Long> rank_w(nprocs, 0), rank_c(nprocs, 0);
    for (int b = 0; b < ba.size(); ++b) {
        const int r = ParallelContext::global_to_local_rank(dm[b]);
        rank_w[r] += wgts[b];
        rank_c[r] += ba[b].numPts();
    }
    auto max_over_mean = [nprocs] (const std::vector<Long>& v) {
        Long tot = 0, mx = 0;
        for (auto x : v) { tot += x; mx = std::max(mx, x); }
        return (tot > 0) ? Real(mx) * nprocs / Real(tot) : one;
    };
    amrex::Print() << "Level " << lev << " load balance (" << load_balance_name()
                   << ", land weight " << load_balance_land_weight << "): " << ba.size()
                   << " boxes, " << total_wet << " of " << total_cells << " cells wet;"
                   << " max/mean rank load " << max_over_mean(rank_w)
                   << " by weight, " << max_over_mean(rank_c) << " by cell count\n";
    return dm;
}

/**
 * @return the remora.load_balance spelling of the current strategy
 */
std::string
REMORA::load_balance_name () const
{
    switch (load_balance) {
        case LoadBalance::amrex:      return "amrex";
        case LoadBalance::sfc:        return "sfc";
        case LoadBalance::knapsack:   return "knapsack";
        case LoadBalance::roundrobin: return "roundrobin";
    }
    return "unknown";
}

/**
 * The level-0 cells whose boxes a finer level needs kept: every refinement indicator's region,
 * in level-0 index space, as AMReX's grid generation sees it. AMReX buffers the tagged cells by
 * amr.n_error_buf, rounds them to blocks of amr.blocking_factor (of level 1) over the refinement
 * ratio, and clears every tag within amr.n_proper such blocks of a cell no level-0 box covers.
 * So a dropped box that close to an indicator would clip the fine grids the full-domain run
 * makes. The region is therefore the indicator's cells plus the buffer, rounded out to whole
 * blocks, grown by n_proper + 1 blocks, and then by the coarse stencil of the coarse-fine fills
 * and the flux register. A box touching it is never dropped, so the fine grids, and the coarse
 * data under and around them, are those of the full-domain run.
 *
 * @param[  out] unbounded  true when an indicator has no region, so a finer level may go anywhere
 */
BoxList
REMORA::level0_refine_footprint (bool& unbounded) const
{
    BoxList bl;
    unbounded = false;
    if (max_level < 1 || ref_tags.empty()) { return bl; }
    unbounded = refine_region_unbounded;
    if (unbounded) { return bl; }

    const IntVect rr = ref_ratio[0];
    IntVect bf(1), stencil(0);
    for (int d = 0; d < 2; ++d) {
        bf[d] = std::max(1, blocking_factor[1][d] / rr[d]);
        stencil[d] = (NGROW + 2 + rr[d] - 1) / rr[d] + 1;
    }
    const Box& dom = geom[0].Domain();
    const Real* plo = geom[0].ProbLo();
    const Real* dx  = geom[0].CellSize();
    for (const auto& rb : refine_regions) {
        IntVect lo = dom.smallEnd(), hi = dom.bigEnd();
        for (int d = 0; d < 2; ++d) {
            lo[d] = static_cast<int>(std::floor((rb.lo(d) - plo[d]) / dx[d])) - n_error_buf[0][d];
            hi[d] = static_cast<int>(std::floor((rb.hi(d) - plo[d]) / dx[d])) + n_error_buf[0][d];
        }
        Box b(lo, hi);
        // Whole blocks, then n_proper + 1 more of them, then the coarse stencil
        b.coarsen(bf);
        b.grow(IntVect(n_proper + 1, n_proper + 1, 0));
        b.refine(bf);
        b.grow(stencil);
        b &= dom;
        if (b.ok()) { bl.push_back(b); }
    }
    return bl;
}

/**
 * The boxes of a level-0 BoxArray that remora.drop_dry_boxes keeps: every box with a wet rho
 * cell in its footprint -- its valid cells and every ghost cell the mask carries, up to one
 * cell beyond a non-periodic domain edge -- and every
 * box a finer level may need (level0_refine_footprint). The rest are land over all of that, so
 * no wet cell's stencil reaches them; partly wet boxes are kept whole.
 *
 * @param[in   ] ba  candidate level-0 BoxArray
 * @return the kept boxes, in the order of ba; ba itself when nothing is dropped
 */
BoxArray
REMORA::drop_dry_level0_boxes (const BoxArray& ba)
{
    BL_PROFILE("REMORA::drop_dry_level0_boxes()");
    if (solverChoice.mask_type == MaskType::none) {
        amrex::Print() << "remora.drop_dry_boxes: no land mask (remora.mask_type = none), so no"
                       << " box is dropped.\n";
        return ba;
    }

    bool unbounded = false;
    const BoxList protect = level0_refine_footprint(unbounded);
    if (unbounded) {
        amrex::Print() << "remora.drop_dry_boxes: a refinement indicator has no box, so a finer"
                       << " level may cover any land; no level-0 box is dropped.\n";
        return ba;
    }

    const BoxArray ba2d = slab_boxes(ba);
    MultiFab mskr(ba2d, DistributionMapping(ba2d), 1, IntVect(NGROW+1,NGROW+1,0));
    if (!probe_rho_mask(0, mskr)) {
        amrex::Print() << "remora.drop_dry_boxes: the level-0 mask comes from remora.hires_grid_level"
                       << " and is not known before the level is made; no box is dropped.\n";
        return ba;
    }

    // The footprint stops one cell beyond a non-periodic domain edge: a grid file holds the
    // mask there and no further, and the rings beyond it are physical-boundary ghost cells,
    // filled from the inside rather than exchanged, so they keep the all-water placeholder.
    Box reach = makeSlab(geom[0].Domain(), 2, 0);
    for (int d = 0; d < 2; ++d) {
        reach.grow(d, geom[0].isPeriodic(d) ? NGROW+1 : 1);
    }
    Vector<int> keep(ba.size(), 0);
    for (MFIter mfi(mskr); mfi.isValid(); ++mfi) {
        const Real wet = mskr[mfi].max<RunOn::Device>(mfi.fabbox() & reach, 0);
        keep[mfi.index()] = (wet > zero) ? 1 : 0;
    }
    ParallelAllReduce::Max(keep.data(), static_cast<int>(keep.size()),
                           ParallelContext::CommunicatorSub());

    BoxList kept;
    Long ndrop = 0, cells_drop = 0;
    for (int b = 0; b < ba.size(); ++b) {
        bool k = keep[b] != 0;
        if (!k) {
            for (Box p : protect) {
                p.setRange(2,0);
                if (p.intersects(ba2d[b])) { k = true; break; }
            }
        }
        if (k) {
            kept.push_back(ba[b]);
        } else {
            ++ndrop;
            cells_drop += ba[b].numPts();
        }
    }
    if (kept.isEmpty()) {
        amrex::Abort("remora.drop_dry_boxes: every level-0 box is land; there is no water to model");
    }

    amrex::Print() << "remora.drop_dry_boxes: level 0 keeps " << kept.size() << " of " << ba.size()
                   << " boxes; dropped " << ndrop << " all-land boxes, " << cells_drop << " of "
                   << ba.numPts() << " cells (" << Real(100.0)*Real(cells_drop)/Real(ba.numPts())
                   << "%)\n";
    if (ndrop == 0) { return ba; }
    return BoxArray(std::move(kept));
}

/**
 * Whether a level's BoxArray leaves part of its domain out, as remora.drop_dry_boxes does to
 * level 0. Only meaningful for level 0, whose boxes are disjoint and should cover the domain.
 */
bool
REMORA::level0_has_holes () const
{
    return grids[0].numPts() != geom[0].Domain().numPts();
}

/**
 * Abort if level 1, or the coarse cells its fills and averaging use, reaches a dropped level-0
 * box. level0_refine_footprint keeps those boxes when the run starts; this catches a grid that
 * outgrew that footprint at a regrid.
 */
void
REMORA::check_dry_box_nesting () const
{
    if (finest_level < 1 || !level0_has_holes()) { return; }
    const IntVect rr = ref_ratio[0];
    IntVect ng(0);
    for (int d = 0; d < 2; ++d) { ng[d] = (NGROW + 2 + rr[d] - 1) / rr[d] + 1; }
    BoxArray need = amrex::coarsen(grids[1], rr);
    need.grow(ng);
    BoxList bl;
    for (int b = 0; b < need.size(); ++b) {
        Box bx = need[b] & geom[0].Domain();
        if (bx.ok()) { bl.push_back(bx); }
    }
    if (!grids[0].contains(BoxArray(std::move(bl)))) {
        amrex::Abort("remora.drop_dry_boxes: level 1, or the level-0 cells it needs, reaches a"
                     " dropped level-0 box. Set remora.drop_dry_boxes = false for this run.");
    }
}

/**
 * Regrid (overrides AmrCore), then make sure level 1 has not grown into a dropped level-0 box.
 */
void
REMORA::regrid (int lbase, Real time, bool initial)
{
    AmrCore::regrid(lbase, time, initial);
    check_dry_box_nesting();
}

/**
 * Level 0's mask, as set_masks built it on the kept boxes, has to be the mask the boxes were
 * dropped on; otherwise the decision rested on a different coastline. Compared everywhere the
 * mask is stored, ghost cells included.
 */
void
REMORA::check_dropped_level0_mask ()
{
    const MultiFab& mskr = *vec_mskr[0];
    MultiFab probe(mskr.boxArray(), mskr.DistributionMap(), 1, mskr.nGrowVect());
    AMREX_ALWAYS_ASSERT(probe_rho_mask(0, probe));
    MultiFab::Subtract(probe, mskr, 0, 0, 1, mskr.nGrowVect());
    if (probe.norminf(0, 1, mskr.nGrowVect()) != zero) {
        amrex::Abort("remora.drop_dry_boxes: level 0's mask differs from the one the boxes were"
                     " dropped on");
    }
}

namespace {
void add_if (Vector<MultiFab*>& v, MultiFab* mf) { if (mf) { v.push_back(mf); } }
void add_if (Vector<MultiFab*>& v, const std::unique_ptr<MultiFab>& mf) { if (mf) { v.push_back(mf.get()); } }
}

/**
 * remora.dry_ghost_test (debug). After initialization, overwrite every ghost cell of a kept
 * level-0 box that lies inside the domain but in no kept box -- the cells no exchange ever
 * fills once boxes are dropped -- in every per-level field: 1 writes a signalling NaN, which
 * with amrex.fpe_trap_invalid = 1 stops the run at the first arithmetic read of one; 2 writes
 * 1e30, so a wet-cell answer that still matches the full-domain run cannot depend on them.
 */
void
REMORA::poison_dry_ghost_cells ()
{
    if (dry_ghost_test == 0 || !level0_has_holes()) { return; }
    const int lev = 0;
    // The fields the time step evolves or rebuilds from them. remora.dry_ghost_fields = all
    // adds the static ones -- grid, metrics, bathymetry, masks, forcing, coefficients -- which
    // initialization computes over every ghost cell itself, so poisoning them tests nothing
    // about the exchanges.
    std::string which = "dynamic";
    ParmParse("remora").query("dry_ghost_fields", which);
    Vector<MultiFab*> f;
    for (auto* m : {cons_new[lev], cons_old[lev], xvel_new[lev], xvel_old[lev],
                    yvel_new[lev], yvel_old[lev], zvel_new[lev], zvel_old[lev]}) { add_if(f, m); }
    for (auto* v : {&vec_ubar, &vec_vbar, &vec_zeta, &vec_Zt_avg1, &vec_Huon, &vec_Hvom,
                    &vec_ru, &vec_rv, &vec_ru2d, &vec_rv2d, &vec_rufrc, &vec_rvfrc, &vec_Akv,
                    &vec_Akt, &vec_sustr, &vec_svstr, &vec_bustr, &vec_bvstr, &vec_stflx,
                    &vec_btflx, &vec_DUon, &vec_DVom, &vec_scratch_2d, &vec_DU_avg1,
                    &vec_DU_avg2, &vec_DV_avg1, &vec_DV_avg2, &vec_Dubar_old, &vec_Dubar_new,
                    &vec_Dvbar_old, &vec_Dvbar_new, &vec_rubar, &vec_rvbar, &vec_rzeta,
                    &vec_sstore, &vec_rhoS, &vec_rhoA, &vec_bvf, &vec_tke, &vec_gls,
                    &vec_Lscale, &vec_Akk, &vec_Akp}) {
        if (lev < v->size()) { add_if(f, (*v)[lev]); }
    }
    if (which == "all") {
        for (auto* v : {&vec_h, &vec_Hz, &vec_z_w, &vec_z_r, &vec_z_phys_nd, &vec_visc2_p,
                        &vec_visc2_r, &vec_diff2, &vec_stflux, &vec_btflux, &vec_pm, &vec_pn,
                        &vec_metric_r, &vec_metric_u, &vec_metric_v, &vec_metric_p, &vec_fcor,
                        &vec_xr, &vec_yr, &vec_xu, &vec_yu, &vec_xv, &vec_yv, &vec_xp, &vec_yp,
                        &vec_dndx, &vec_dmde, &vec_rdrag, &vec_rdrag2, &vec_ZoBot, &vec_uwind,
                        &vec_vwind, &vec_Tair, &vec_qair, &vec_Pair, &vec_srflx,
                        &vec_longwave_down, &vec_cloud, &vec_EminusP, &vec_alpha, &vec_beta,
                        &vec_lrflx, &vec_lhflx, &vec_shflx, &vec_rain, &vec_evap, &vec_mskr,
                        &vec_msku, &vec_mskv, &vec_mskp, &vec_mskr3d}) {
            if (lev < v->size()) { add_if(f, (*v)[lev]); }
        }
        for (auto& nv : vec_nudg_coeff) {
            if (lev < nv.size()) { add_if(f, nv[lev]); }
        }
    }

    Real big = Real(1.e30);
    ParmParse("remora").query("dry_ghost_value", big);
    const Real val = (dry_ghost_test == 1) ? std::numeric_limits<Real>::signaling_NaN() : big;
    Long npoison = 0;
    for (MultiFab* mf : f) {
        const IndexType it = mf->ixType();
        Box dom = amrex::convert(geom[lev].Domain(), it);
        BoxArray cov = amrex::convert(grids[lev], it);
        if (mf->boxArray().minimalBox().length(2) == 1) {
            // A 2D field lives on slabs at k = 0
            dom.setRange(2,0);
            BoxList bl = cov.boxList();
            for (auto& b : bl) { b.setRange(2,0); }
            cov = BoxArray(std::move(bl));
        }
        for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {
            const Box fab = mfi.fabbox() & dom;
            if (!fab.ok()) { continue; }
            const BoxList holes = cov.complementIn(fab);
            Array4<Real> const& a = mf->array(mfi);
            const int nc = mf->nComp();
            for (const Box& h : holes) {
                npoison += h.numPts() * nc;
                ParallelFor(h, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    a(i,j,k,n) = val;
                });
            }
        }
    }
    Gpu::streamSynchronize();
    ParallelDescriptor::ReduceLongSum(npoison);
    amrex::Print() << "remora.dry_ghost_test = " << dry_ghost_test << ": overwrote " << npoison
                   << " uncovered ghost values in " << f.size() << " level-0 fields\n";
}

/**
 * With level-0 boxes dropped, the ghost cells of a kept box that lie in a dropped box are never
 * filled by an exchange; they hold whatever the box's own kernels leave there, or what the
 * memory held when the field was allocated. The land stencils of the kept box read them, but
 * what they compute stays on land -- a wet cell is at least NGROW+2 land cells away -- so the
 * answer does not depend on the values (the regression tests overwrite them and compare). It
 * has to be a number, though: a NaN or an infinity there would survive the masks (0 * NaN is
 * NaN). So every FArrayBox is zeroed when it is allocated, as AMReX debug builds already do,
 * unless the input set fab.do_initval itself.
 */
void
REMORA::zero_new_fabs_if_holes ()
{
    if (!level0_has_holes()) { return; }
    if (amrex::InitSNaN()) {
        // AMReX then skips fab.do_initval, so those cells start as signalling NaNs
        amrex::Print() << "WARNING: amrex.init_snan with remora.drop_dry_boxes: the ghost cells"
                       << " of kept boxes inside dropped ones start as signalling NaNs, and"
                       << " the land stencils of the kept boxes read them, so"
                       << " amrex.fpe_trap_invalid will stop the run.\n";
        return;
    }
    ParmParse pp("fab");
    if (pp.contains("do_initval")) { return; }
    if (!FArrayBox::get_do_initval() || FArrayBox::get_initval() != zero) {
        FArrayBox::set_do_initval(true);
        FArrayBox::set_initval(zero);
        amrex::Print() << "remora.drop_dry_boxes: new FArrayBoxes are zeroed (fab.do_initval = 1,"
                       << " fab.initval = 0), so the ghost cells no exchange reaches hold a"
                       << " number\n";
    }
}

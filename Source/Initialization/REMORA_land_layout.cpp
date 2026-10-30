/**
 * \file REMORA_land_layout.cpp
 *
 * How a level's boxes are laid out over land: the wet-cell weights the distribution mapping
 * can balance (remora.load_balance). They need the land mask of a layout before the level on
 * it exists, which probe_rho_mask provides.
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

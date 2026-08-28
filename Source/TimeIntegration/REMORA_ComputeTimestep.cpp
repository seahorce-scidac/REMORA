#include <REMORA.H>

#include <cmath>

using namespace amrex;

void
REMORA::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = estTimeStep(lev);
    }

    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = amrex::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = amrex::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // dt_0 is identical on every rank here, so this aborts collectively. The
    // bogus_large_value test is the load-bearing one: the sentinel used to survive to
    // this point intact whenever stop_time was left at its default of Real::max(), and to
    // be laundered into stop_time - t_new[0] -- a single step covering the entire run --
    // whenever stop_time was set. The change_max cap above cannot catch either case,
    // since dt[] is itself seeded to bogus_large_value.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt_0 > zero && std::isfinite(dt_0) &&
                                     dt_0 < bogus_large_value,
                                     "REMORA::ComputeDt: computed a non-positive, "
                                     "non-finite, or unusably large level-0 dt");

    // Limit dt's by the value of stop_time.
    const Real eps = Real(1.e-3)*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

/**
 * Estimate the largest stable slow (baroclinic) timestep on a level.
 *
 * The estimate is the smaller of an advective limit and an external gravity wave limit,
 * each a maximum over the wet cells of the level:
 *
 *     dt_adv  = cfl / max( |u|*pm, |v|*pn, |w|/Hz )
 *     dt_grav = cfl * ndtfast / max( sqrt(g*|h|) * sqrt(pm^2 + pn^2) )
 *
 * The second is the Courant quantity ROMS forms as Cg_max in metrics.F. It is what keeps
 * the estimate finite for a run started from rest, where every velocity is zero and the
 * advective limit alone is undefined.
 *
 * pm and pn are the per-cell 1/dx and 1/dy metric terms rather than the geometry's
 * uniform cell size, so stretched and curvilinear grids give the right answer.
 *
 * @param[in] level    level of refinement
 */
Real
REMORA::estTimeStep(int level) const
{
    BL_PROFILE("REMORA::estTimeStep()");

    // The barotropic mode is substepped ndtfast times per baroclinic step
    // (REMORA_Advance.cpp), so the slow step only has to resolve the external gravity
    // wave to within that ratio. ReadParameters guarantees ndtfast is positive.

    // g is a file-scope constexpr; hoist it into a local for the device lambda.
    const Real grav = g;

    MultiFab ccvel(grids[level],dmap[level],3,0);

    average_face_to_cellcenter(ccvel,0,
        Array<const MultiFab*,3>{xvel_new[level], yvel_new[level], zvel_new[level]});

    // Two maxima over the same cells. amrex::ReduceMax's FabArray overloads take at most
    // three FabArrays and this needs six, so drive ReduceOps directly. Unlike
    // amrex::ReduceMax, ReduceOps::eval has no host fallback in a GPU build, so no
    // Gpu::LaunchSafeGuard is wanted here.
    ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(ccvel, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // ccvel carries no ghost cells, so this is exactly its valid region -- the same
        // cells the advective estimate has always visited.
        const Box& bx = mfi.tilebox();

        Array4<Real const> const& u    = ccvel.const_array(mfi);
        Array4<Real const> const& Hz   = vec_Hz[level]->const_array(mfi);

        // h, pm, pn and mskr live on the z-slab BoxArray built box-for-box from
        // grids[level] with the same DistributionMapping, so they index off this same
        // MFIter and are read at k = 0.
        Array4<Real const> const& h    = vec_h[level]->const_array(mfi);
        Array4<Real const> const& pm   = vec_pm[level]->const_array(mfi);
        Array4<Real const> const& pn   = vec_pn[level]->const_array(mfi);
        Array4<Real const> const& mskr = vec_mskr[level]->const_array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            // Land contributes nothing: h there is a fill depth, not a real one.
            if (mskr(i,j,0) <= Real(0.5)) { return {zero, zero}; }

            Real inv_adv = amrex::max(amrex::Math::abs(u(i,j,k,0)) * pm(i,j,0),
                                      amrex::Math::abs(u(i,j,k,1)) * pn(i,j,0));

            // Hz is the true thickness of this terrain-following layer; the geometry's
            // uniform dz is the thickness of a sigma level, not of a cell.
            if (Hz(i,j,k) > zero) {
                inv_adv = amrex::max(inv_adv, amrex::Math::abs(u(i,j,k,2)) / Hz(i,j,k));
            }

            // Component 0 of vec_h is the bathymetry, positive down. Component 1 is not a
            // second copy of it -- stretch_transform overwrites it with the bottom z_w.
            const Real c      = std::sqrt(grav * amrex::Math::abs(h(i,j,0,0)));
            const Real inv_bt = c * std::sqrt(pm(i,j,0)*pm(i,j,0) + pn(i,j,0)*pn(i,j,0));

            return {inv_adv, inv_bt};
        });
    } // mfi

    ReduceTuple host_tuple = reduce_data.value(reduce_op);

    // A rank owning no boxes never calls eval and leaves its local maxima at
    // std::numeric_limits<Real>::lowest(); the MPI max repairs that, and every test below
    // is against the reduced values, which are identical on every rank.
    Real inv[2] = { amrex::get<0>(host_tuple), amrex::get<1>(host_tuple) };
    ParallelDescriptor::ReduceRealMax(inv, 2);
    const Real inv_adv = inv[0];
    const Real inv_bt  = inv[1];

    // Guard both divisions rather than letting them produce inf: inv_adv is zero for any
    // quiescent start, and inv_bt is zero for a level that is entirely land.
    const Real estdt_adv = (inv_adv > zero) ? cfl / inv_adv                : bogus_large_value;
    const Real estdt_bt  = (inv_bt  > zero) ? cfl * Real(ndtfast) / inv_bt : bogus_large_value;

    const Real estdt_lowM = amrex::min(estdt_adv, estdt_bt);

    if (verbose) {
        amrex::Print() << "Using cfl = " << cfl << std::endl;
        if (inv_adv > zero) {
            amrex::Print() << "  advective    limit at level " << level << ":  " << estdt_adv << std::endl;
        } else {
            amrex::Print() << "  advective    limit at level " << level << ":  none (velocity is zero)" << std::endl;
        }
        if (inv_bt > zero) {
            amrex::Print() << "  gravity wave limit at level " << level << ":  " << estdt_bt
                           << "  (ndtfast = " << ndtfast << ")" << std::endl;
        } else {
            amrex::Print() << "  gravity wave limit at level " << level << ":  none (no wet cell has a depth)" << std::endl;
        }
        if (fixed_dt > zero) {
            if (estdt_lowM < bogus_large_value) {
                amrex::Print() << "Slow  dt at level " << level << " would be:  " << estdt_lowM << std::endl;
            } else {
                amrex::Print() << "Slow  dt at level " << level << " would be undefined " << std::endl;
            }
            amrex::Print() << "Fixed dt at level " << level << "       is:  " << fixed_dt << std::endl;
        } else {
            amrex::Print() << "Slow  dt at level " << level << ":  " << estdt_lowM << std::endl;
        }
    }

    if (fixed_dt > zero) {
        return fixed_dt;
    }

    // Water at rest still carries gravity waves, so estdt_bt is finite for any level with
    // one wet cell of nonzero depth. Reaching here means there is no such cell; returning
    // the sentinel would let ComputeDt clamp dt to stop_time - t_new[0] and run the whole
    // simulation in a single step.
    if (estdt_lowM >= bogus_large_value) {
        amrex::Abort("REMORA::estTimeStep: cannot estimate a timestep at level " +
                     std::to_string(level) + ". No wet cell (mskr > 0.5) has a nonzero "
                     "depth, so neither the advective nor the gravity wave limit is "
                     "defined. Check the bathymetry and the land mask, or set "
                     "remora.fixed_dt");
    }

    return estdt_lowM;
}

#include <iomanip>

#include "REMORA.H"
#include "AMReX_MultiFab.H"

using namespace amrex;
/**
 * @param[in   ] time    current time
 */
void
REMORA::sum_integrated_quantities(Real time)
{
    BL_PROFILE("REMORA::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    int datwidth = 14;
    // Six digits is enough to watch a run by eye but not to compare two runs: an averaged-down
    // bathymetry shifts the volume by ~1e-5 relative. Raise it when the sums are being used as
    // a diagnostic to assert on rather than to read.
    int datprecision = 6;
    amrex::ParmParse pp("remora");
    pp.queryAdd("sum_precision", datprecision);
    bool local = true;

    // One sum per cell-centered tracer past salinity: the passive scalars and the
    // biology tracers. Empty when the run carries neither.
    const int n_tracers = ncons - Tracer_comp;

    Vector<Real> tracer_ml(n_tracers, zero);
    Real kineng_ml = zero;
    Real volume_ml = zero;
    Real max_vel_ml = zero;

    Vector<Real> tracer_sl(n_tracers, zero);
    Real kineng_sl = zero;
    Real volume_sl = zero;
    Real max_vel_sl = zero;

    for (int t = 0; t < n_tracers; ++t) {
        tracer_sl[t] = volWgtSumMF(0,*cons_new[0],Tracer_comp+t,local,false);
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab kineng_mf(grids[lev], dmap[lev], 1, 0);
        MultiFab ones_mf(grids[lev], dmap[lev], 1, 0);
        ones_mf.setVal(one);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<      Real> kineng_arr = kineng_mf.array(mfi);
            const Array4<const Real> xvel_u_arr = xvel_new[lev]->const_array(mfi);
            const Array4<const Real> yvel_v_arr = yvel_new[lev]->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // This is the same expression for kinetic energy that is used in ROMS
                kineng_arr(i,j,k) = fourth * ( xvel_u_arr(i,j,k)*xvel_u_arr(i,j,k) + xvel_u_arr(i+1,j,k)*xvel_u_arr(i+1,j,k) +
                                                yvel_v_arr(i,j,k)*yvel_v_arr(i,j,k) + yvel_v_arr(i  ,j+1,k)*yvel_v_arr(i,j+1,k));
            });
        } // mfi

        const int icomp = 0;
        Real max_vel_local = std::sqrt(two * kineng_mf.max(icomp));

        if (lev==0) {
          kineng_sl = volWgtSumMF(lev,kineng_mf,0,local,false);
          volume_sl = volWgtSumMF(lev,ones_mf  ,0,local,false);
          max_vel_sl = max_vel_local;
        }

        for (int t = 0; t < n_tracers; ++t) {
            tracer_ml[t] += volWgtSumMF(lev,*cons_new[lev],Tracer_comp+t,local,true);
        }
        kineng_ml += volWgtSumMF(lev,kineng_mf     ,             0,local,true);
        volume_ml += volWgtSumMF(lev,ones_mf       ,             0,local,true);
        max_vel_ml = std::max(max_vel_ml, max_vel_local);
    }

    if (verbose > 0) {
        // Layout: every tracer's single-level sum, then kineng/volume, then the same
        // for the multi-level sums. Unpacked below in the same order.
        Vector<Real> sum_vars;
        sum_vars.reserve(2*n_tracers + 4);
        for (int t = 0; t < n_tracers; ++t) { sum_vars.push_back(tracer_sl[t]); }
        sum_vars.push_back(kineng_sl);
        sum_vars.push_back(volume_sl);
        for (int t = 0; t < n_tracers; ++t) { sum_vars.push_back(tracer_ml[t]); }
        sum_vars.push_back(kineng_ml);
        sum_vars.push_back(volume_ml);
        const int n_sum_vars = static_cast<int>(sum_vars.size());

        const int n_max_vars = 2;
        Real max_vars[n_max_vars] = {max_vel_sl, max_vel_ml};
#ifdef AMREX_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
        ParallelDescriptor::ReduceRealSum(
            sum_vars.data(), n_sum_vars, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealMax(
            max_vars, n_max_vars, ParallelDescriptor::IOProcessorNumber());

          if (ParallelDescriptor::IOProcessor()) {
            int i = 0;
            for (int t = 0; t < n_tracers; ++t) { tracer_sl[t] = sum_vars[i++]; }
            kineng_sl = sum_vars[i++];
            volume_sl = sum_vars[i++];
            for (int t = 0; t < n_tracers; ++t) { tracer_ml[t] = sum_vars[i++]; }
            kineng_ml = sum_vars[i++];
            volume_ml = sum_vars[i++];
            int j = 0;
            max_vel_sl = max_vars[j++];
            max_vel_ml = max_vars[j++];

            // Label field is as wide as the widest tracer name so the columns line up
            // whether the run carries "tracer" or "phytoplankton"
            int namewidth = 9;
            for (int t = 0; t < n_tracers; ++t) {
                namewidth = std::max(namewidth, static_cast<int>(cons_names[Tracer_comp+t].size()));
            }

            if (finest_level == 0) {
                amrex::Print() << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "TIME" << std::right
                               << " = " << std::setw(datwidth) << std::setprecision(datprecision) << time << '\n';
                for (int t = 0; t < n_tracers; ++t) {
                    amrex::Print() << std::setw(namewidth) << std::left << cons_names[Tracer_comp+t] << std::right
                                   << " = " << std::setw(datwidth) << std::setprecision(datprecision) << tracer_sl[t] << '\n';
                }
                amrex::Print() << std::setw(namewidth) << std::left << "KIN. ENG." << std::right
                               << " = " << std::setw(datwidth) << std::setprecision(datprecision) << kineng_sl  << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "VOLUME" << std::right
                               << " = " << std::setw(datwidth) << std::setprecision(datprecision) << volume_sl  << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "MAX. VEL." << std::right
                               << " = " << std::setw(datwidth) << std::setprecision(datprecision) << max_vel_sl << '\n';
            } else {
                amrex::Print() << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "TIME" << std::right
                               << "       = " << std::setw(datwidth) << std::setprecision(datprecision) << time << '\n';
                for (int t = 0; t < n_tracers; ++t) {
                    amrex::Print() << std::setw(namewidth) << std::left << cons_names[Tracer_comp+t] << std::right
                                   << " SL/ML = " << std::setw(datwidth) << std::setprecision(datprecision) << tracer_sl[t] << ' '
                                                  << std::setw(datwidth) << std::setprecision(datprecision) << tracer_ml[t] << '\n';
                }
                amrex::Print() << std::setw(namewidth) << std::left << "KIN. ENG." << std::right
                               << " SL/ML = " << std::setw(datwidth) << std::setprecision(datprecision) << kineng_sl  << ' '
                                              << std::setw(datwidth) << std::setprecision(datprecision) << kineng_ml << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "VOLUME" << std::right
                               << " SL/ML = " << std::setw(datwidth) << std::setprecision(datprecision) << volume_sl  << ' '
                                              << std::setw(datwidth) << std::setprecision(datprecision) << volume_ml << '\n';
                amrex::Print() << std::setw(namewidth) << std::left << "MAX. VEL." << std::right
                               << " SL/ML = " << std::setw(datwidth) << std::setprecision(datprecision) << max_vel_sl << ' '
                                              << std::setw(datwidth) << std::setprecision(datprecision) << max_vel_ml << '\n';
            }

            if (NumDataLogs() > 0) {
                std::ostream& data_log1 = DataLog(0);
                if (data_log1.good()) {
                    if (time == zero) {
                        data_log1 << std::setw(datwidth) << "time";
                        for (int t = 0; t < n_tracers; ++t) {
                            data_log1 << std::setw(datwidth) << cons_names[Tracer_comp+t];
                        }
                        data_log1 << std::setw(datwidth) << "kineng";
                        data_log1 << std::setw(datwidth) << "volume";
                        data_log1 << std::setw(datwidth) << "max_vel";
                        data_log1 << std::endl;
                    }

                  // Write the quantities at this time
                  data_log1 << std::setw(datwidth) << time;
                  for (int t = 0; t < n_tracers; ++t) {
                      data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                                << tracer_ml[t];
                  }
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << kineng_ml;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << volume_ml;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << max_vel_ml;
                  data_log1 << std::endl;
              }
            }
          }
#ifdef AMREX_LAZY
        });
#endif
    }
}

/**
 * @param[in   ] lev      level to calculate on
 * @param[in   ] mf       data to sum over
 * @param[in   ] comp     component on which to calculate sum
 * @param[in   ] local    whether to do the sum locally
 * @param[in   ] finemask whether to mask fine level
 */
Real
REMORA::volWgtSumMF(int lev, const MultiFab& mf, int comp, bool local, bool finemask)
{
    BL_PROFILE("REMORA::volWgtSumMF()");

    Real sum = zero;
    MultiFab tmp(grids[lev], dmap[lev], 1, 0);
    MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

    if (lev < finest_level && finemask) {
        const MultiFab& mask = build_fine_mask(lev+1);
        MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
    }

    MultiFab volume(grids[lev], dmap[lev], 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> vol_arr = volume.array(mfi);
        const Array4<const Real>      Hz = vec_Hz[lev]->const_array(mfi);
        const Array4<const Real>      pm = vec_pm[lev]->const_array(mfi);
        const Array4<const Real>      pn = vec_pn[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vol_arr(i,j,k) = Hz(i,j,k) / (pm(i,j,0) * pn(i,j,0));
        });
    } // mfi

    sum = MultiFab::Dot(tmp, 0, volume, 0, 1, 0, local);

    if (!local)
      ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

/**
 * @param[in   ] lev      level to calculate on
 */
MultiFab&
REMORA::build_fine_mask(int level)
{
    // Mask for zeroing covered cells
    AMREX_ASSERT(level > 0);

    const BoxArray& cba = grids[level-1];
    const DistributionMapping& cdm = dmap[level-1];

    // TODO -- we should make a vector of these a member of REMORA class
    fine_mask.define(cba, cdm, 1, 0, MFInfo());
    fine_mask.setVal(one);

    BoxArray fba = grids[level];
    iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, ref_ratio[level-1], 1, 0);

    const auto  fma =  fine_mask.arrays();
    const auto ifma = ifine_mask.arrays();
    ParallelFor(fine_mask, [=] AMREX_GPU_DEVICE(int bno, int i, int j, int k) noexcept
    {
        fma[bno](i,j,k) = ifma[bno](i,j,k);
    });

    Gpu::synchronize();

    return fine_mask;
}

/**
 * @param[in   ] nstep              what step we're on
 * @param[in   ] lev                level to calculate on
 * @param[in   ] dtlev              time step for this level
 * @param[in   ] action_interval    number of time steps between actions
 * @param[in   ] action_per         time interval between actions
 */
bool
REMORA::is_it_time_for_action(int nstep, Real time, Real dtlev, int action_interval, Real action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > zero) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}

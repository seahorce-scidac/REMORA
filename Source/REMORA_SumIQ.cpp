#include <iomanip>

#include "REMORA.H"
#include "AMReX_MultiFab.H"

using namespace amrex;

void
REMORA::sum_integrated_quantities(Real time)
{
    BL_PROFILE("REMORA::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    int datwidth = 14;
    int datprecision = 6;

    Real scalar = 0.0_rt;
    Real kineng = 0.0_rt;
    Real volume = 0.0_rt;
    Real max_vel = 0.0_rt;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab kineng_mf(grids[lev], dmap[lev], 1, 0);
        MultiFab ones_mf(grids[lev], dmap[lev], 1, 0);
        ones_mf.setVal(1.0_rt);

        for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<      Real> kineng_arr = kineng_mf.array(mfi);
            const Array4<const Real> xvel_u_arr = xvel_new[lev]->const_array(mfi);
            const Array4<const Real> yvel_v_arr = yvel_new[lev]->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // This is the same expression for kinetic energy that is used in ROMS
                kineng_arr(i,j,k) = 0.25_rt * ( xvel_u_arr(i,j,k)*xvel_u_arr(i,j,k) + xvel_u_arr(i+1,j,k)*xvel_u_arr(i+1,j,k) +
                                                yvel_v_arr(i,j,k)*yvel_v_arr(i,j,k) + yvel_v_arr(i  ,j+1,k)*yvel_v_arr(i,j+1,k));
            });
        } // mfi

        const int icomp = 0;
        Real max_vel_local = std::sqrt(2.0_rt * kineng_mf.max(icomp));

        scalar += volWgtSumMF(lev,*cons_new[lev],Scalar_comp,false,true);
        kineng += volWgtSumMF(lev,kineng_mf     ,             0,false,true);
        volume += volWgtSumMF(lev,ones_mf       ,             0,false,true);
        max_vel = std::max(max_vel, max_vel_local);
    }

    if (verbose > 0) {
        const int n_sum_vars = 3;
        Real sum_vars[n_sum_vars] = {scalar,kineng,volume};

        const int n_max_vars = 1;
        Real max_vars[n_max_vars] = {max_vel};
#ifdef AMREX_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
        ParallelDescriptor::ReduceRealSum(
            sum_vars, n_sum_vars, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealMax(
            max_vars, n_max_vars, ParallelDescriptor::IOProcessorNumber());

          if (ParallelDescriptor::IOProcessor()) {
            int i = 0;
            scalar = sum_vars[i++];
            kineng = sum_vars[i++];
            volume = sum_vars[i++];
            int j = 0;
            max_vel = max_vars[j++];

            amrex::Print() << '\n';
            amrex::Print() << "TIME= " << time << " SCALAR      = " << scalar  << '\n';
            amrex::Print() << "TIME= " << time << " KIN. ENG.   = " << kineng  << '\n';
            amrex::Print() << "TIME= " << time << " VOLUME      = " << volume  << '\n';
            amrex::Print() << "TIME= " << time << " MAX. VEL.   = " << max_vel << '\n';

            if (NumDataLogs() > 0) {
                std::ostream& data_log1 = DataLog(0);
                if (data_log1.good()) {
                    if (time == 0.0_rt) {
                        data_log1 << std::setw(datwidth) << "          time";
                        data_log1 << std::setw(datwidth) << "        scalar";
                        data_log1 << std::setw(datwidth) << "        kineng";
                        data_log1 << std::setw(datwidth) << "        volume";
                        data_log1 << std::setw(datwidth) << "       max_vel";
                        data_log1 << std::endl;
                    }

                  // Write the quantities at this time
                  data_log1 << std::setw(datwidth) << time;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << scalar;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << kineng;
                  data_log1 << std::endl;
              }
            }
          }
#ifdef AMREX_LAZY
        });
#endif
    }
}

Real
REMORA::volWgtSumMF(int lev, const MultiFab& mf, int comp, bool local, bool finemask)
{
    BL_PROFILE("REMORA::volWgtSumMF()");

    Real sum = 0.0_rt;
    MultiFab tmp(grids[lev], dmap[lev], 1, 0);
    MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

    if (lev < finest_level && finemask) {
        const MultiFab& mask = build_fine_mask(lev+1);
        MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
    }

    MultiFab volume(grids[lev], dmap[lev], 1, 0);
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

MultiFab&
REMORA::build_fine_mask(int level)
{
    // Mask for zeroing covered cells
    AMREX_ASSERT(level > 0);

    const BoxArray& cba = grids[level-1];
    const DistributionMapping& cdm = dmap[level-1];

    // TODO -- we should make a vector of these a member of REMORA class
    fine_mask.define(cba, cdm, 1, 0, MFInfo());
    fine_mask.setVal(1.0_rt);

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

bool
REMORA::is_it_time_for_action(int nstep, Real time, Real dtlev, int action_interval, Real action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > 0.0_rt) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}

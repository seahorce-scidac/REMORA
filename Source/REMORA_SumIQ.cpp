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

    Real scalar = 0.0;
    Real kineng = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab kineng_mf(grids[lev], dmap[lev], 1, 0);
        MultiFab cc_vel_mf(grids[lev], dmap[lev], 3, 0);
        average_face_to_cellcenter(cc_vel_mf,0,
            Array<const MultiFab*,3>{xvel_new[lev],yvel_new[lev],zvel_new[lev]});

        for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<      Real> kineng_arr = kineng_mf.array(mfi);
            const Array4<const Real>    vel_arr = cc_vel_mf.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                kineng_arr(i,j,k) = 0.5 * ( vel_arr(i,j,k,0)*vel_arr(i,j,k,0) + vel_arr(i,j,k,1)*vel_arr(i,j,k,1) +
                                            vel_arr(i,j,k,2)*vel_arr(i,j,k,2) );
            });
        } // mfi

        scalar += volWgtSumMF(lev,*cons_new[lev],Scalar_comp,false,true);
        kineng += volWgtSumMF(lev,kineng_mf     ,             0,false,true);
    }

    if (verbose > 0) {
        const int nfoo = 2;
        Real foo[nfoo] = {scalar,kineng};
#ifdef AMREX_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
        ParallelDescriptor::ReduceRealSum(
            foo, nfoo, ParallelDescriptor::IOProcessorNumber());

          if (ParallelDescriptor::IOProcessor()) {
            int i = 0;
            scalar = foo[i++];
            kineng = foo[i++];

            amrex::Print() << '\n';
            amrex::Print() << "TIME= " << time << " SCALAR      = " << scalar << '\n';
            amrex::Print() << "TIME= " << time << " KIN. ENG.   = " << kineng << '\n';

            if (NumDataLogs() > 0) {
                std::ostream& data_log1 = DataLog(0);
                if (data_log1.good()) {
                    if (time == 0.0) {
                        data_log1 << std::setw(datwidth) << "          time";
                        data_log1 << std::setw(datwidth) << "        scalar";
                        data_log1 << std::setw(datwidth) << "        kineng";
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

    Real sum = 0.0;
    MultiFab tmp(grids[lev], dmap[lev], 1, 0);
    MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

    if (lev < finest_level && finemask) {
        const MultiFab& mask = build_fine_mask(lev+1);
        MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
    }

    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    auto const& dx = geom[lev].CellSizeArray();
    Real cell_vol = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
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
    fine_mask.setVal(1.0);

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
  if (action_per > 0.0) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}

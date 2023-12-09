#include <EOS.H>
#include <ROMSX.H>

using namespace amrex;

// a wrapper for estTimeStep
void
ROMSX::ComputeDt ()
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

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

Real
ROMSX::estTimeStep(int level) const
{
  BL_PROFILE("ROMSX::estTimeStep()");

  amrex::Real estdt_lowM = 1.e20;

  auto const dxinv = geom[level].InvCellSizeArray();

  MultiFab ccvel(grids[level],dmap[level],3,0);

  average_face_to_cellcenter(ccvel,0,
      Array<const MultiFab*,3>{xvel_new[level], yvel_new[level], zvel_new[level]});

   Real estdt_lowM_inv = amrex::ReduceMax(ccvel, 0,
       [=] AMREX_GPU_HOST_DEVICE (Box const& b,
                                  Array4<Real const> const& u) -> Real
       {
           Real new_lm_dt = -1.e100;
           amrex::Loop(b, [=,&new_lm_dt] (int i, int j, int k) noexcept
           {
               new_lm_dt = amrex::max(((amrex::Math::abs(u(i,j,k,0)))*dxinv[0]),
                                      ((amrex::Math::abs(u(i,j,k,1)))*dxinv[1]),
                                      ((amrex::Math::abs(u(i,j,k,2)))*dxinv[2]), new_lm_dt);
           });
           return new_lm_dt;
       });

   amrex::ParallelDescriptor::ReduceRealMax(estdt_lowM_inv);
   if (estdt_lowM_inv > 0.0_rt)
       estdt_lowM = cfl / estdt_lowM_inv;;

  if (verbose) {
    if (fixed_dt <= 0.0) {
        amrex::Print() << "Using cfl = " << cfl << std::endl;
        if (estdt_lowM_inv > 0.0_rt) {
            amrex::Print() << "Slow  dt at level " << level << ":  " << estdt_lowM << std::endl;
        } else {
            amrex::Print() << "Slow  dt at level " << level << ": undefined " << std::endl;
        }
    }
    if (fixed_dt > 0.0) {
        amrex::Print() << "Based on cfl of 1.0 " << std::endl;
        if (estdt_lowM_inv > 0.0_rt) {
            amrex::Print() << "Slow  dt at level " << level << " would be:  " << estdt_lowM/cfl << std::endl;
        } else {
            amrex::Print() << "Slow  dt at level " << level << " would be undefined " << std::endl;
        }
        amrex::Print() << "Fixed dt at level " << level << "       is:  " << fixed_dt << std::endl;
    }
  }

  if (fixed_dt > 0.0) {
    return fixed_dt;
  } else {
        return estdt_lowM;
  }
}

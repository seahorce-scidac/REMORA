#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

using namespace amrex;

namespace derived {

void
remora_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Array4<const amrex::Real>& /*pm*/,
  const amrex::Array4<const amrex::Real>& /*pn*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine does nothing -- we use it as a placeholder.
}

void
remora_dervort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Array4<const amrex::Real>& pm,
  const amrex::Array4<const amrex::Real>& pn,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int level)
{
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real d2x = 0.5_rt / pm(i-1,j,  0) + 1.0_rt / pm(i,j,0) + 0.5_rt / pm(i+1,j,0);
        Real d2y = 0.5_rt / pn(i,  j-1,0) + 1.0_rt / pm(i,j,0) + 0.5_rt / pm(i,j+1,0);
        tfab(i,j,k,dcomp) = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (d2x)  // dv/dx
                          - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (d2y); // du/dy
    });
}
}


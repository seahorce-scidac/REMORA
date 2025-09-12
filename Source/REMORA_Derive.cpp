#include "REMORA_Derive.H"
#include "REMORA_IndexDefines.H"

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

/*
 * \brief Function to calculate vorticity for derived field output
 *
 * @param[in   ] bx         box to calculate on
 * @param[  out] derfab     derived field output
 * @param[in   ] dcomp      component to store derived field
 * @param[in   ] ncomp      number of components
 * @param[in   ] datfab     velocity to be used to calculate derived field
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] geomdata   geometry data
 * @param[in   ] time       current time
 * @param[in   ] bcrec      BC info
 * @param[in   ] level      current level
 */
void
remora_dervort(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Array4<const amrex::Real>& pm,
  const amrex::Array4<const amrex::Real>& pn,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real d2x = 0.5_rt / pm(i-1,j,  0) + 1.0_rt / pm(i,j,0) + 0.5_rt / pm(i+1,j,0);
        Real d2y = 0.5_rt / pn(i,  j-1,0) + 1.0_rt / pn(i,j,0) + 0.5_rt / pm(i,j+1,0);
        tfab(i,j,k,dcomp) = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (d2x)  // dv/dx
                          - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (d2y); // du/dy
    });
}
}


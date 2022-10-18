#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

namespace derived {

void romsx_derrhodivide(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  const amrex::FArrayBox& datfab,
  const int scalar_index)
{
  // This routine divides any cell-centered conserved quantity by density
  auto const dat = datfab.array();
  auto primitive  = dromsxab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho       = dat(i, j, k, Rho_comp);
    const amrex::Real conserved = dat(i, j, k, scalar_index);
    primitive(i,j,k) = conserved / rho;
  });
}

void
romsx_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*dromsxab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  // This routine does nothing -- we use it as a placeholder.
}

void
romsx_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto pfab      = dromsxab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    AMREX_ALWAYS_ASSERT(rhotheta > 0.);
    pfab(i,j,k) = getPgivenRTh(rhotheta);
  });
}

void
romsx_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto cfab      = dromsxab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    const amrex::Real rho      = dat(i, j, k, Rho_comp);
    AMREX_ALWAYS_ASSERT(rhotheta > 0.);
    cfab(i,j,k) = std::sqrt(Gamma * getPgivenRTh(rhotheta) / rho);
  });
}

void
romsx_derscalar(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  romsx_derrhodivide(bx, dromsxab, datfab, RhoScalar_comp);
}

void
romsx_derKE(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  romsx_derrhodivide(bx, dromsxab, datfab, RhoKE_comp);
}

void
romsx_derQKE(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  romsx_derrhodivide(bx, dromsxab, datfab, RhoQKE_comp);
}

void
romsx_deromega(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  romsx_derrhodivide(bx, dromsxab, datfab, Omega_comp);
}

void
romsx_dersalt(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  romsx_derrhodivide(bx, dromsxab, datfab, Salt_comp);
}
}

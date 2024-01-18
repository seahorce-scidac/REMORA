#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

namespace derived {

void
remora_dernull(
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*dremoraab*/,
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
}


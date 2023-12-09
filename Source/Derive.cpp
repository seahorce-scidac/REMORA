#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

namespace derived {

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
}


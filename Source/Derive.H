#ifndef _Derive_H_
#define _Derive_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

namespace derived {

void romsx_dernull(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);
}
#endif


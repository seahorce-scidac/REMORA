#ifndef _Derive_H_
#define _Derive_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

namespace derived {

void romsx_derrhodivide(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  const amrex::FArrayBox& datfab,
  const int scalar_index);

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

void romsx_derpres(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void romsx_dersoundspeed(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void romsx_derscalar(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void romsx_derKE(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void romsx_derQKE(
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

void romsx_deromega(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

void romsx_dersalt(
  const amrex::Box& bx,
  amrex::FArrayBox& dromsxab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);
#endif

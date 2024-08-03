#ifndef REMORA_PhysBCFunct_H_
#define REMORA_PhysBCFunct_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FilCC_C.H>
#include <AMReX_FilND_C.H>
#include <AMReX_FilFC_C.H>
#include "AMReX_TypeTraits.H"
#include "AMReX_Orientation.H"

#include <TimeInterpolatedData.H>
#include <IndexDefines.H>
#include <DataStruct.H>

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

struct NullFill
{
    AMREX_GPU_DEVICE
    void operator() (amrex::Array4<amrex::Real> const& /*dest*/,
                     amrex::Array4<amrex::Real> const& /*mask*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
        }
};

class REMORAPhysBCFunct
{
public:
    REMORAPhysBCFunct (const int lev,
                      const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NCONS> bc_extdir_vals
                    )
        : m_lev(lev), m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals)
          {}

    ~REMORAPhysBCFunct () {}

    //
    // mf is the multifab to be filled
    // icomp is the index into the MultiFab -- if cell-centered this can be any value
    //       from 0 to NCONS-1, if face-centered this must be 0
    // ncomp is the number of components -- if cell-centered this can be any value
    //       from 1 to NCONS as long as icomp+ncomp <= NCONS-1.  If face-centered this
    //       must be 1
    // nghost is how many ghost cells to be filled
    // time is the time at which the data should be filled
    // bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
    //     so this follows the BCVars enum
    //
    void operator() (amrex::MultiFab& mf, const amrex::MultiFab& mask, int icomp, int ncomp, amrex::IntVect const& nghost,
                     amrex::Real time, int bccomp, int n_not_fill=0);

    void impose_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& msku,
                          amrex::Real time, int bccomp);

    void impose_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskv,
                          amrex::Real time, int bccomp);


    void impose_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                          const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskr,
                          amrex::Real time, int bccomp);

    void impose_cons_bcs (const amrex::Array4<amrex::Real>& mf, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskr,
                          int icomp, int ncomp, amrex::Real time, int bccomp, int n_not_fill);

    // For backward compatibility
    // void FillBoundary (amrex::MultiFab& mf, int dcomp, int ncomp, amrex::IntVect const& nghost,
    //                    amrex::Real time, int bccomp) {
    //     this->operator()(mf,dcomp,ncomp,nghost,time,bccomp);
    // }

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NCONS> m_bc_extdir_vals;
};

#endif

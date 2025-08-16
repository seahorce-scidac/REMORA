#ifndef _REMORA_STRETCH_H_
#define _REMORA_STRETCH_H_

#include <cmath>
#include <REMORA_DataStruct.H>
#include <REMORA.H>
#include <REMORA_prob_common.H>

using namespace amrex;

void
REMORA::calc_stretch_coeffs ()
{
    BL_PROFILE("REMORA::calc_stretch_coeffs()");
    int nz = geom[0].Domain().length(2);

    auto N = nz; // Number of vertical "levels" aka, NZ
    Real ds = 1.0_rt / Real(N);

    Box bx(IntVect(0,0,0),IntVect(0,0,N));
    const auto local_theta_s = solverChoice.theta_s;
    const auto local_theta_b = solverChoice.theta_b;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int , int , int k)
    {
        Real Csur,Cbot;

        s_w[k]=ds*(k-N);
        if (local_theta_s > 0.0_rt) {
            Csur=(1.0_rt-std::cosh(local_theta_s*s_w[k]))/
                (std::cosh(local_theta_s)-1.0_rt);
        } else {
            Csur=-s_w[k]*s_w[k];
        }

        if (local_theta_b > 0.0_rt) {
            Cbot=(std::exp(local_theta_b*Csur)-1.0_rt)/
                (1.0_rt-std::exp(-local_theta_b));
            Cs_w[k]=Cbot;
        } else {
            Cs_w[k]=Csur;
        }

        if (k<N) {
            s_r[k]=ds*(k-N+0.5_rt);

            if (local_theta_s > 0.0_rt) {
                Csur=(1.0_rt-std::cosh(local_theta_s*s_r[k]))/
                    (std::cosh(local_theta_s)-1.0_rt);
            } else {
                Csur=-s_r[k]*s_r[k];
            }

            if (local_theta_b > 0.0_rt) {
                Cbot=(std::exp(local_theta_b*Csur)-1.0_rt)/
                    (1.0_rt-std::exp(-local_theta_b));
                Cs_r[k]=Cbot;
            } else {
                Cs_r[k]=Csur;
            }
        }
    });
}

/**
 * @param[in] lev    level to operate on
 */
void
REMORA::stretch_transform (int lev)
{
    BL_PROFILE("REMORA::stretch_transform()");
    std::unique_ptr<MultiFab>& mf_z_w = vec_z_w[lev];
    std::unique_ptr<MultiFab>& mf_z_r = vec_z_r[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = vec_Hz[lev];
    std::unique_ptr<MultiFab>& mf_h  = vec_hOfTheConfusingName[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];
    std::unique_ptr<MultiFab>& mf_z_phys_nd  = vec_z_phys_nd[lev];

    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& z_w = (mf_z_w)->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& Hz  = (mf_Hz)->array(mfi);
        Array4<Real> const& h  = (mf_h)->array(mfi);
        Array4<Real> const& Zt_avg1  = (mf_Zt_avg1)->array(mfi);
        Box bx = mfi.tilebox();
        Box gbx2 = bx;
        gbx2.grow(IntVect(NGROW,NGROW,0));
        Box gbx3 = bx;
        gbx3.grow(IntVect(NGROW+1,NGROW+1,0));
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);
        Box gbx3D = gbx3;
        gbx3D.makeSlab(2,0);
        Box wgbx3 = gbx3;
        wgbx3.surroundingNodes(2);

        const auto & geomdata = Geom(lev).data();

        int nz = geom[lev].Domain().length(2);

        auto N = nz; // Number of vertical "levels" aka, NZ
        //forcing tcline to be the same as probhi for now, one in DataStruct.H other in inputs
        Real hc=-min(geomdata.ProbHi(2),-solverChoice.tcline); // Do we need to enforce min here?
        Real ds = 1.0_rt / Real(N);

        amrex::ParallelFor(wgbx3, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            z_w(i,j,k) = h(i,j,0);
        });

        // ROMS Transform 2
        Gpu::streamSynchronize();

        amrex::ParallelFor(wgbx3, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real sc_r=ds*(k-N+0.5_rt);
            Real sc_w=ds*(k-N);

            Real cff_r=hc*sc_r;
            Real cff1_r=Cs_r[k];

            Real cff_w=hc*sc_w;
            Real cff1_w=Cs_w[k];

            Real hwater=h(i,j,0);

            Real hinv=1.0_rt/(hc+hwater);
            Real cff2_r=(cff_r+cff1_r*hwater)*hinv;
            Real cff2_w=(cff_w+cff1_w*hwater)*hinv;

            if(k==N) {
            // HACK: should actually be the normal expression with coeffs evaluated at k=N-1
                z_w(i,j,N)=Zt_avg1(i,j,0);

            } else if (k==0) {
                  h(i,j,0,1) = Zt_avg1(i,j,0)+(Zt_avg1(i,j,0)+hwater)*cff2_w;
                  z_w(i,j,0) = h(i,j,0,1);

            } else {
                  z_w(i,j,k)=Zt_avg1(i,j,0)+(Zt_avg1(i,j,0)+hwater)*cff2_w;

            }

            if(k!=N) {
                  z_r(i,j,k)=Zt_avg1(i,j,0)+(Zt_avg1(i,j,0)+hwater)*cff2_r;
            }
        });

        Gpu::streamSynchronize();

        amrex::ParallelFor(gbx3, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              Hz(i,j,k)=z_w(i,j,k+1)-z_w(i,j,k);
        });
    } // mfi

    vec_z_w[lev]->FillBoundary(geom[lev].periodicity());
    vec_z_r[lev]->FillBoundary(geom[lev].periodicity());
    vec_Hz[lev]->FillBoundary(geom[lev].periodicity());

    // Define nodal z as average of z on w-faces
    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& z_w = (mf_z_w)->array(mfi);
        Array4<Real> const& z_phys_nd  = (mf_z_phys_nd)->array(mfi);

        Box z_w_box = Box(z_w);
        auto const lo = amrex::lbound(z_w_box);
        auto const hi = amrex::ubound(z_w_box);

        //
        // NOTE: we assume that all boxes extend the full extent of the domain in the vertical direction
        //
        // NOTE: z_phys_nd(i,j,k) and z_w(i,j,k) both refer to the node on the LOW  side of cell (i,j,k)
        //

        // We shrink the box in the vertical since z_phys_nd has a ghost cell in the vertical
        //    which we will fill by extrapolation, not from z_w
        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));//Box(z_phys_nd); bx.grow(2,-1);

        ParallelFor(convert(bx,IntVect(0,0,1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // For now assume all boundaries are constant height --
            //     we will enforce periodicity below
            if ( i >= lo.x && i <= hi.x-1 && j >= lo.y && j <= hi.y-1 )
            {
                z_phys_nd(i,j,k)=0.25_rt*( z_w(i,j  ,k) + z_w(i+1,j  ,k) +
                                         z_w(i,j+1,k) + z_w(i+1,j+1,k) );
            } else {
                int ii = std::min(std::max(i, lo.x), hi.x);
                int jj = std::min(std::max(j, lo.y), hi.y);
                z_phys_nd(i,j,k) = z_w(ii,jj,k);
            }
        });

        // Fill nodes below the surface (to avoid out of bounds errors in particle functions)
        int klo = -1;
        ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            z_phys_nd(i,j,klo) = 2.0_rt * z_phys_nd(i,j,klo+1) - z_phys_nd(i,j,klo+2);
        });
    } // mf

    // Note that we do *not* want to do a multilevel fill here -- we have
    // already filled z_phys_nd on the grown boxes, but we enforce periodicity just in case
    vec_z_phys_nd[lev]->FillBoundary(geom[lev].periodicity());
}

#endif

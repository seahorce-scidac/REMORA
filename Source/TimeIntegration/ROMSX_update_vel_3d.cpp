#include <ROMSX.H>

using namespace amrex;

//
// update_vel_3d -- called from prestep_uv_3d or prestep_t_3d
// NOTE: "vel" here represents either u, v, or a tracer
// When updating u, ioff=1, joff=0
// When updating v, ioff=0, joff=1
// When updating tracer, ioff=0, joff=0
// The tracer update is a bit different from the u,v updates so we test
// for it, but checking if ioff=0 and joff=0. In some cases, though, we
// can recover the tracer update from the generic one by setting those indices.
// Setting icc and ntfirst identically for the tracers should be equivalent
// to setting ioff=0 and joff=0
//

void
ROMSX::update_vel_3d (const Box& vel_bx, const Box& gbx,
                      const int ioff, const int joff,
                      Array4<Real>  vel, Array4<Real> vel_old,
                      Array4<Real> rvel, Array4<Real> Hz,
                      Array4<Real>  Akv,
                      Array4<Real>   DC, Array4<Real> FC,
                      Array4<Real> sstr, Array4<Real> bstr,
                      Array4<Real> z_r,
                      Array4<Real>   pm, Array4<Real>  pn,
                      const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                      const Real lambda, const Real dt_lev)
{

    BoxArray ba_gbxvel = intersect(BoxArray(vel_bx), gbx);
    AMREX_ASSERT((ba_gbxvel.size() == 1));
    Box gbxvel = ba_gbxvel[0];

    //
    //  Weighting coefficient for the newest (implicit) time step derivatives
    //  using either a Crack-Nicolson implicit scheme (lambda=0.5) or a
    //  backward implicit scheme (lambda=1.0).
    //


    Real oml_dt = dt_lev*(1.0-lambda);
    //N is one less than ROMS

    //  Except the commented out part means lambda is always 1.0
    if (verbose > 1) {
        amrex::Print() << "in update_vel_3d with box " << vel_bx << std::endl;
        Print() << "vel old " << Box(vel_old) << std::endl;
        Print() << "Akv " << Box(Akv) << std::endl;
    }
    ParallelFor(vel_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if(k+1<=N&&k>=0)
        {
            Real cff = 1.0 / ( z_r(i,j,k+1)+z_r(i-ioff,j-joff,k+1)
                              -z_r(i,j,k  )-z_r(i-ioff,j-joff,k  ));
            FC(i,j,k) = oml_dt * cff * (vel_old(i,j,k+1,nstp)-vel_old(i,j,k,nstp)) *
                                           (Akv(i,j,k)     +Akv(i-ioff,j-joff,k));
        }
        else if (k==-1)
        {
            FC(i,j,-1) = dt_lev*bstr(i,j,0);
        }
        else if (k==N)
        {
            FC(i,j, N) = dt_lev*sstr(i,j,0);
        }

        DC(i,j,k) = 0.25 * dt_lev * (pm(i,j,0)+pm(i-ioff,j-joff,0))
                                      * (pn(i,j,0)+pn(i-ioff,j-joff,0));
    });

    ParallelFor(gbxvel,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff1 = 0.0, cff2 = 0.0, cff4;

        int indx=0; //nrhs-3

        if (iic==ntfirst)
        {
            if (k+1<=N && k>=1)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
                vel(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==0)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=FC(i,j,k)-dt_lev*bstr(i,j,0);
                vel(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==N)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=dt_lev*sstr(i,j,0)-FC(i,j,k-1); //or: -FC(i,j,k-1)+dt_lev*sstr(i,j,0);
                vel(i,j,k,nnew)=cff1+cff2;
            }

        } else if(iic==ntfirst+1) {

            if (k<N && k>0) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
            }
            else if(k==0) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=FC(i,j,k)-dt_lev*bstr(i,j,0);
            }
            else if(k==N) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff2=dt_lev*sstr(i,j,0)-FC(i,j,k-1); //or: -FC(i,j,k-1)+dt_lev*sstr(i,j,0);
            }

            cff3=0.5*DC(i,j,k);
            if(ioff==0&&joff==0)
                vel(i,j,k,nnew)=cff1 + cff2;
            else {
                indx=nrhs ? 0 : 1;
                Real r_swap= rvel(i,j,k,indx);
                rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
                rvel(i,j,k,nrhs) = r_swap;
                vel(i,j,k,nnew)=cff1- cff3*rvel(i,j,k,indx)+ cff2;
            }
        } else {
            cff1= 5.0/12.0;
            cff2=16.0/12.0;
            if (k==0) {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff4=FC(i,j,k)-dt_lev*bstr(i,j,0);
                //cff4=FC(i,j,k)-FC(i,j,k-1);//-bustr(i,j,0);

            } else if (k == N) {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff4=dt_lev*sstr(i,j,0)-FC(i,j,k-1);

            } else {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                cff4=FC(i,j,k)-FC(i,j,k-1);
            }
            if(ioff==0&&joff==0)
                vel(i,j,k,nnew)=cff3 + cff4;
            else {
                indx=nrhs ? 0 : 1;
                Real r_swap= rvel(i,j,k,indx);
                rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
                rvel(i,j,k,nrhs) = r_swap;

                vel(i,j,k,nnew)=cff3+
                    DC(i,j,k)*(cff1*rvel(i,j,k,nrhs)-
                               cff2*rvel(i,j,k,indx))+
                    cff4;
                rvel(i,j,k,nrhs) = 0.0;
            }
        }
    });
}

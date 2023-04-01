#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// update_vel_3d -- called from prestep_uv_3d
// NOTE: "vel" here represents either u or v
//
void
ROMSX::update_vel_3d (const Box& bx,
                      const int ioff, const int joff,
                      Array4<Real>  vel_arr , Array4<Real> vel_old,
                      Array4<Real> rvel_arr , Array4<Real> Hz_arr,
                      Array4<Real> DC_arr, Array4<Real> FC_arr,
                      Array4<Real> sstr_arr,
                      const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                      const Real lambda, const Real dt_lev)
{
    // Copy and grow the tilebox
    Box gbx1 = bx;
    gbx1.grow(IntVect(1,1,0));

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff1, cff2, cff4;

        int indx=0; //nrhs-3

        if (iic==ntfirst)
        {
            //Hz still might need adjusting
            if (k+1<=N && k>=1)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=FC_arr(i,j,k)-FC_arr(i,j,k-1);
                vel_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==0)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=FC_arr(i,j,k);//-bustr(i,j,0);
                vel_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==N)
            {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=-FC_arr(i,j,k-1)+dt_lev*sstr_arr(i,j,0);
                vel_arr(i,j,k,nnew)=cff1+cff2;
            }

        } else if(iic==ntfirst+1) {

            if (k<N && k>0) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=FC_arr(i,j,k)-FC_arr(i,j,k-1);
            }
            else if(k==0) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=FC_arr(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff1=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff2=-FC_arr(i,j,k-1)+dt_lev*sstr_arr(i,j,0);
            }
            cff3=0.5*DC_arr(i,j,k);
            indx=nrhs ? 0 : 1;
            Real r_swap= rvel_arr(i,j,k,indx);
            rvel_arr(i,j,k,indx) = rvel_arr(i,j,k,nrhs);
            rvel_arr(i,j,k,nrhs) = r_swap;
            vel_arr(i,j,k,nnew)=cff1- cff3*rvel_arr(i,j,k,indx)+ cff2;

        } else {

            cff1= 5.0/12.0;
            cff2=16.0/12.0;
            if (k<N && k>0) {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff4=FC_arr(i,j,k)-FC_arr(i,j,k-1);
            }
            else if(k==0) {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff4=FC_arr(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff3=vel_old(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-ioff,j-joff,k));
                cff4=-FC_arr(i,j,k-1)+dt_lev*sstr_arr(i,j,0);
            }
            indx=nrhs ? 0 : 1;
            Real r_swap= rvel_arr(i,j,k,indx);
            rvel_arr(i,j,k,indx) = rvel_arr(i,j,k,nrhs);
            rvel_arr(i,j,k,nrhs) = r_swap;
            vel_arr(i,j,k,nnew)=cff3+
                DC_arr(i,j,k)*(cff1*rvel_arr(i,j,k,nrhs)-
                           cff2*rvel_arr(i,j,k,indx))+
                cff4;
            rvel_arr(i,j,k,nrhs) = 0.0;
        }
    });
}

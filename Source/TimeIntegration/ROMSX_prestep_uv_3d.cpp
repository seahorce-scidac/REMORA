#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//
void
ROMSX::prestep_uv_3d (const Box& bx,
                      Array4<Real> uold  , Array4<Real> vold,
                      Array4<Real> u_arr , Array4<Real> v_arr,
                      Array4<Real> ru_arr, Array4<Real> rv_arr,
                      Array4<Real> Hz_arr, Array4<Real> Akv_arr,
                      Array4<Real> on_u, Array4<Real> om_v,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm  , Array4<Real> pn,
                      Array4<Real> W   , Array4<Real> DC,
                      Array4<Real> FC  , Array4<Real> z_r_arr,
                      Array4<Real> sustr_arr, Array4<Real> svstr_arr,
                      int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                      Real lambda, Real dt_lev)
{
    //copy the tilebox
    Box gbx1 = bx;
    Box gbx11 = bx;
    Box gbx2 = bx;

    Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
    Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(2,2,0));
    gbx1.grow(IntVect(1,1,0));
    gbx11.grow(IntVect(1,1,1));

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
          //-----------------------------------------------------------------------
          //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
          //-----------------------------------------------------------------------
            if (k+1<=N) {
                if (i-1>=-2)
                {
                    Huon(i,j,k)=0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k)) * uold(i,j,k,nrhs) * on_u(i,j,0);
                } else {
                    Huon(i,j,k)=(Hz_arr(i,j,k))*uold(i,j,k,nrhs) * on_u(i,j,0);
                }

                if (j-1>=-2)
                {
                    Hvom(i,j,k)=0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k))*vold(i,j,k,nrhs)* om_v(i,j,0);
                } else {
                    Hvom(i,j,k)=(Hz_arr(i,j,k))*vold(i,j,k,nrhs)* om_v(i,j,0);
                }
            }
    });

    //Should really use gbx3uneven
    amrex::ParallelFor(gbx2uneven,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //
        //------------------------------------------------------------------------
        //  Vertically integrate horizontal mass flux divergence.
        //------------------------------------------------------------------------
        //
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        if(k==0) {
            W(i,j,k)=0.0;
        } else {
            W(i,j,k) = W(i,j,k-1)- (Huon(i+1,j,k)-Huon(i,j,k)+ Hvom(i,j+1,k)-Hvom(i,j,k));
        }
    });

    //Need to include pre_step3d.F terms

    //
    //  Weighting coefficient for the newest (implicit) time step derivatives
    //  using either a Crack-Nicolson implicit scheme (lambda=0.5) or a
    //  backward implicit scheme (lambda=1.0).
    //

    //  Except the commented out part means lambda is always 1.0
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff;

        if(k+1<=N&&k>=0)
        {
            cff = 1.0 / (z_r_arr(i,j,k+1)+z_r_arr(i-1,j,k+1)- z_r_arr(i,j,k  )-z_r_arr(i-1,j,k  ));

            FC(i,j,k)=cff3*cff*(uold(i,j,k,nstp)-uold(i,j,k-1,nstp))* (Akv_arr(i,j,k)+Akv_arr(i-1,j,k));
        }
        else
        {
            //  FC(i,j,-1)=0.0;//dt_lev*bustr(i,j,0);
            //  FC(i,j,N)=0.0;//dt_lev*sustr_arr(i,j,0);
        }
        cff=dt_lev*.25;
        DC(i,j,k)=cff*(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff1, cff2, cff4;

        int indx=0; //nrhs-3

        if (iic==ntfirst)
        {
            //Hz still might need adjusting
            if(k+1<=N&&k>=1)
            {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
                u_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==0)
            {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=FC(i,j,k);//-bustr(i,j,0);
                u_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==N)
            {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=-FC(i,j,k-1)+dt_lev*sustr_arr(i,j,0);
                u_arr(i,j,k,nnew)=cff1+cff2;
            }

        } else if(iic==ntfirst+1) {

            if(k<N&&k>0) {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
            }
            else if(k==0) {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=FC(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff1=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff2=-FC(i,j,k-1)+dt_lev*sustr_arr(i,j,0);
            }
            cff3=0.5*DC(i,j,k);
            indx=nrhs ? 0 : 1;
            Real r_swap= ru_arr(i,j,k,indx);
            ru_arr(i,j,k,indx) = ru_arr(i,j,k,nrhs);
            ru_arr(i,j,k,nrhs) = r_swap;
            u_arr(i,j,k,nnew)=cff1- cff3*ru_arr(i,j,k,indx)+ cff2;

        } else {

            cff1= 5.0/12.0;
            cff2=16.0/12.0;
            if(k<N&&k>0) {
                cff3=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff4=FC(i,j,k)-FC(i,j,k-1);
            }
            else if(k==0) {
                cff3=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff4=FC(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff3=uold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i-1,j,k));
                cff4=-FC(i,j,k-1)+dt_lev*sustr_arr(i,j,0);
            }
            indx=nrhs ? 0 : 1;
            Real r_swap= ru_arr(i,j,k,indx);
            ru_arr(i,j,k,indx) = ru_arr(i,j,k,nrhs);
            ru_arr(i,j,k,nrhs) = r_swap;
            u_arr(i,j,k,nnew)=cff3+
                DC(i,j,k)*(cff1*ru_arr(i,j,k,nrhs)-
                           cff2*ru_arr(i,j,k,indx))+
                cff4;
            ru_arr(i,j,k,nrhs) = 0.0;
        }
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff;

        if(k+1<=N&&k>=0)
        {
            cff=1.0/(z_r_arr(i,j,k+1)+z_r_arr(i,j-1,k+1)-
                     z_r_arr(i,j,k  )-z_r_arr(i,j-1,k  ));
            FC(i,j,k)=cff3*cff*(vold(i,j,k,nstp)-vold(i,j,k-1,nstp))*
                (Akv_arr(i,j,k)+Akv_arr(i,j-1,k));
        }
        else
        {
            //              FC(i,j,-1)=0.0;//dt_lev*bustr(i,j,0);
            //              FC(i,j,N)=0.0;//dt_lev*sustr_arr(i,j,0);
        }
        cff=dt_lev*.25;
        DC(i,j,k)=cff*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff3=dt_lev*(1.0-lambda);
        Real cff1, cff2, cff4;

        int indx=0; //nrhs-3

        if(iic==ntfirst)
        {
            //if(j>0&&j<Mm-1)
            {
            //Hz still might need adjusting
            if(k+1<=N&&k>=1)
            {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
                v_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==0)
            {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=FC(i,j,k);//-bustr(i,j,0);
                v_arr(i,j,k,nnew)=cff1+cff2;
            }
            else if(k==N)
            {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=-FC(i,j,k-1)+dt_lev*svstr_arr(i,j,0);
                v_arr(i,j,k,nnew)=cff1+cff2;
            } }

        } else if(iic==ntfirst+1) {

            if(k<N&&k>0) {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=FC(i,j,k)-FC(i,j,k-1);
            }
            else if(k==0) {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=FC(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff1=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff2=-FC(i,j,k-1)+dt_lev*svstr_arr(i,j,0);
            }
            cff3=0.5*DC(i,j,k);
            indx=nrhs ? 0 : 1;
            Real r_swap= rv_arr(i,j,k,indx);
            rv_arr(i,j,k,indx) = rv_arr(i,j,k,nrhs);
            rv_arr(i,j,k,nrhs) = r_swap;
            //if(j>0&&j<Mm-1)
                v_arr(i,j,k,nnew) = cff1 - cff3*rv_arr(i,j,k,indx) + cff2;

        } else {
            cff1= 5.0/12.0;
            cff2=16.0/12.0;
            if(k<N&&k>0) {
                cff3=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff4=FC(i,j,k)-FC(i,j,k-1);
            }
            else if(k==0) {
                cff3=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff4=FC(i,j,k);//-bustr(i,j,0);
            }
            else if(k==N) {
                cff3=vold(i,j,k,nstp)*0.5*(Hz_arr(i,j,k)+Hz_arr(i,j-1,k));
                cff4=-FC(i,j,k-1)+dt_lev*svstr_arr(i,j,0);
            }
            indx=nrhs ? 0 : 1;
            Real r_swap= rv_arr(i,j,k,indx);
            rv_arr(i,j,k,indx) = rv_arr(i,j,k,nrhs);
            rv_arr(i,j,k,nrhs) = r_swap;
            //if(j>0&&j<Mm-1)
            v_arr(i,j,k,nnew)=cff3+
                DC(i,j,k)*(cff1*rv_arr(i,j,k,nrhs)-
                           cff2*rv_arr(i,j,k,indx))+
                cff4;

            rv_arr(i,j,k,nrhs) = 0.0;

        }
    });
}

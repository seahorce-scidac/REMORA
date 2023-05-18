#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::rho_eos (const Box& phi_bx,
                Array4<Real> rho,
                Array4<Real> rhoA, 
                Array4<Real> rhoS,
                Array4<Real> pden,
		Array4<Real> Hz,
		const int N)
{
    const int Mn = Geom(0).Domain().size()[0];
    const int Mm = Geom(0).Domain().size()[1];
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);

//
//=======================================================================
//  Linear equation of state.
//=======================================================================
//
//-----------------------------------------------------------------------
//  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
//  equation of state.
//-----------------------------------------------------------------------
//
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
	rho(i,j,k)=R0-                                          
	    R0*Tcoef*(t(i,j,k,nrhs,itemp)-T0);
	rho(i,j,k)=rho(i,j,k)+                                      
	    R0*Scoef*(t(i,j,k,nrhs,isalt)-S0);
	rho(i,j,k)=rho(i,j,k)-1000.0_r8;
	pden(i,j,k)=rho(i,j,k);
    });
//
//-----------------------------------------------------------------------
//  Compute vertical averaged density (rhoA) and density perturbation
//  used (rhoS) in barotropic pressure gradient.
//-----------------------------------------------------------------------
//
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
	cff1=rho(i,j,N)*Hz(i,j,N);
	rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N);
	rhoA(i,j)=cff1;
    });
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
	if(k!=0) {
        cff1=rho(i,j,N-1+k)*Hz(i,j,N-1+k);
	rhoS(i,j)=rhoS(i,j)+Hz(i,j,N-1+k)*(rhoA(i,j)+0.5_r8*cff1);
	rhoA(i,j)=rhoA(i,j)+cff1;
	}
    });
    amrex::ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
	cff2=1.0_r8/rho0;
        cff1=1.0_r8/(z_w(i,j,N)-z_w(i,j,0));
	rhoA(i,j)=cff2*cff1*rhoA(i,j);
	rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j);
    });

}

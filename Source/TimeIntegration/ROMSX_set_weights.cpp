#include <cmath>
#include <DataStruct.H>
#include <ROMSX.H>
#include <prob_common.H>

using namespace amrex;

//  This routine sets the weigth functions for the time averaging of
//  2D fields over all short time-steps.
void ROMSX::set_weights (int lev) {

    int i,j,iter;
    Real gamma, scale;
    Real cff1, cff2;
    Real wsum, shift, cff;

    //HACK should possibly store fixed_ndtfast elsewhere
    int ndtfast=fixed_ndtfast_ratio>0 ? fixed_ndtfast_ratio : fixed_fast_dt / fixed_dt;

    //From mod_scalars
    Real Falpha = 2.0_rt;
    Real Fbeta = 4.0_rt;
    Real Fgamma = 0.284_rt;

    vec_weight1.resize(2*ndtfast+1);
    vec_weight2.resize(2*ndtfast+1);

    auto weight1 = vec_weight1.dataPtr();
    auto weight2 = vec_weight2.dataPtr();

//
//=======================================================================
//  Compute time-averaging filter for barotropic fields.
//=======================================================================
//
//  Initialize both sets of weights to zero.
//
    nfast=0;
    for(int i=1;i<=2*ndtfast;i++) {
        weight1[i-1]=0.0_rt;
        weight2[i-1]=0.0_rt;
    }
//
//-----------------------------------------------------------------------
//  Power-law shape filters.
//-----------------------------------------------------------------------
//
//  The power-law shape filters are given by:
//
//     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
//
//  where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
//  normalization are chosen to yield the correct zeroth-order
//  (normalization), first-order (consistency), and second-order moments,
//  resulting in overall second-order temporal accuracy for time-averaged
//  barotropic motions resolved by baroclinic time step.
//
      scale=(Falpha+1.0_rt)*(Falpha+Fbeta+1.0_rt)/
      ((Falpha+2.0_rt)*(Falpha+Fbeta+2.0_rt)*Real(ndtfast));
//
//  Find center of gravity of the primary weighting shape function and
//  iteratively adjust "scale" to place the  centroid exactly at
//  "ndtfast".
//
      gamma=Fgamma*max(0.0_rt, 1.0_rt-10.0_rt/Real(ndtfast));
      for(iter=1;iter<=16;iter++) {
          nfast=0;
      for(int i=1;i<=2*ndtfast;i++) {
          cff=scale*Real(i);
          weight1[i-1]=pow(cff,Falpha)-pow(cff,(Falpha+Fbeta))-gamma*cff;
          if (weight1[i-1]>0.0_rt) nfast=i;
          if ((nfast>0)&&(weight1[i-1]<0.0_rt)) {
            weight1[i-1]=0.0_rt;
          }
    }
        wsum=0.0_rt;
        shift=0.0_rt;
        for(i=1;i<=nfast;i++) {
          wsum=wsum+weight1[i-1];
          shift=shift+weight1[i-1]*Real(i);
        }
        scale=scale*shift/(wsum*Real(ndtfast));
      }
//
//-----------------------------------------------------------------------
//  Post-processing of primary weights.
//-----------------------------------------------------------------------
//
//  Although it is assumed that the initial settings of the primary
//  weights has its center of gravity "reasonably close" to NDTFAST,
//  it may be not so according to the discrete rules of integration.
//  The following procedure is designed to put the center of gravity
//  exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
//  applying basically an upstream advection of weights to eliminate
//  the mismatch iteratively. Once this procedure is complete primary
//  weights are normalized.
//
//  Find center of gravity of the primary weights and subsequently
//  calculate the mismatch to be compensated.
//
      for(iter=1;iter<=ndtfast;iter++) {
        wsum=0.0_rt;
        shift=0.0_rt;
        for(i=1;i<=nfast;i++) {
          wsum=wsum+weight1[i-1];
          shift=shift+Real(i)*weight1[i-1];
        }
        shift=shift/wsum;
        cff=Real(ndtfast)-shift;
//
//  Apply advection step using either whole, or fractional shifts.
//  Notice that none of the four loops here is reversible.
//
        if (cff>1.0_rt) {
          nfast=nfast+1;
          for(i=nfast;i>=2;i--) {
            weight1[i-1]=weight1[i-1-1];
          }
          weight1[1-1]=0.0_rt;
        } else if (cff>0.0_rt) {
          wsum=1.0_rt-cff;
          for(i=nfast;i>=2;i--) {
          weight1[i-1]=wsum*weight1[i-1]+cff*weight1[i-1-1];
          }
          weight1[1-1]=wsum*weight1[1-1];
        } else if (cff<-1.0_rt) {
          nfast=nfast-1;
          for(i=1;i<=nfast;i++) {
          weight1[i-1]=weight1[i+1-1];
          }
          weight1[nfast+1-1]=0.0_rt;
        } else if (cff<0.0_rt) {
          wsum=1.0_rt+cff;
          for(i=1;i<=nfast-1;i++) {
          weight1[i-1]=wsum*weight1[i-1]-cff*weight1[i+1-1];
          }
          weight1[nfast-1]=wsum*weight1[nfast-1];
        }
      }
//
//  Set SECONDARY weights assuming that backward Euler time step is used
//  for free surface.  Notice that array weight2[i] is assumed to
//  have all-zero status at entry in this segment of code.
//
        for(j=1;j<=nfast;j++) {
        cff=weight1[j-1];
        for(i=1;i<=j;i++) {
            weight2[i-1]=weight2[i-1]+cff;
        }
      }
//
//  Normalize both set of weights.
//
      wsum=0.0_rt;
      cff=0.0_rt;
      for(i=1;i<=nfast;i++) {
        wsum=wsum+weight1[i-1];
        cff=cff+weight2[i-1];
      }
      wsum=1.0_rt/wsum;
      cff=1.0_rt/cff;
      for(i=1;i<=nfast;i++) {
        weight1[i-1]=wsum*weight1[i-1];
        weight2[i-1]=cff*weight2[i-1];
      }
//
//  Report weights.
//
    if (ParallelDescriptor::IOProcessor()) {
        Print().SetPrecision(18)<<ParallelDescriptor::NProcs()<<"  "<<ndtfast<<"  "<<nfast<<"  "<<wsum<<std::endl;
        cff=0.0_rt;
        cff1=0.0_rt;
        cff2=0.0_rt;
        wsum=0.0_rt;
        shift=0.0_rt;
        for(i=1;i<=nfast;i++) {
          cff=cff+weight1[i-1];
          cff1=cff1+weight1[i-1]*Real(i);
          cff2=cff2+weight1[i-1]*Real(i*i);
          wsum=wsum+weight2[i-1];
          shift=shift+weight2[i-1]*(Real(i)-0.5_rt);
          Print().SetPrecision(18)<<"i="<<i<<"  "<<weight1[i-1]<<"  "<<weight2[i-1]<<"  "<<cff<<"  "<<wsum<<std::endl;
        }
        cff1=cff1/Real(ndtfast);
        cff2=cff2/(Real(ndtfast)*Real(ndtfast));
        shift=shift/Real(ndtfast);
        Print().SetPrecision(18)<<ndtfast <<"  "<< nfast<<"  "<<Real(nfast)/Real(ndtfast)<<std::endl;
        Print().SetPrecision(18)<<cff1<<"  "<<cff2<<"  "<<shift<<"  "<<cff<<"  "<<wsum<<"  "<<Fgamma<<"  "<<gamma<<std::endl;
      if (cff2<1.0001_rt) Print()<<"\n\n\n"<<std::endl;
      }
}

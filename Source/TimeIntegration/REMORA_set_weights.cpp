#include <cmath>
#include <REMORA_DataStruct.H>
#include <REMORA.H>
#include <REMORA_prob_common.H>

using namespace amrex;

/**
 * @param[in   ] lev     level to operate on
 */
void REMORA::set_weights (int /*lev*/) {

    Real gamma, scale;
    Real wsum, shift, cff;

    //HACK should possibly store fixed_ndtfast elsewhere
    int ndtfast=fixed_ndtfast_ratio>0 ? fixed_ndtfast_ratio : static_cast<int>(fixed_fast_dt / fixed_dt);

    //From mod_scalars
    Real Falpha = two;
    Real Fbeta = Real(4.0);
    Real Fgamma = Real(0.284);

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
        weight1[i-1]=zero;
        weight2[i-1]=zero;
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
    scale=(Falpha+one)*(Falpha+Fbeta+one) /
        ((Falpha+two)*(Falpha+Fbeta+two)*Real(ndtfast));
    //
    //  Find center of gravity of the primary weighting shape function and
    //  iteratively adjust "scale" to place the  centroid exactly at
    //  "ndtfast".
    //
    gamma = Fgamma*max(zero, one-Real(10.0)/Real(ndtfast));

    for (int iter=1;iter<=16;iter++) {
        nfast=0;
        for(int i=1;i<=2*ndtfast;i++) {
            cff=scale*Real(i);

            weight1[i-1]=Real(pow(cff,Falpha)-pow(cff,(Falpha+Fbeta)))-gamma*cff;

            if (weight1[i-1] > zero) {
                nfast=i;
            }

            if ( (nfast>0) && (weight1[i-1] < zero) ) {
                weight1[i-1] = zero;
            }
        }
        wsum  = zero;
        shift = zero;
        for(int i=1;i<=nfast;i++) {
            wsum=wsum+weight1[i-1];
            shift=shift+weight1[i-1]*Real(i);
        }
        scale *= shift/(wsum*Real(ndtfast));
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
    for (int iter=1;iter<=ndtfast;iter++) {
        wsum  = zero;
        shift = zero;
        for(int i=1;i<=nfast;i++) {
            wsum=wsum+weight1[i-1];
            shift=shift+Real(i)*weight1[i-1];
        }
        shift=shift/wsum;
        cff=Real(ndtfast)-shift;
        //
        //  Apply advection step using either whole, or fractional shifts.
        //  Notice that none of the four loops here is reversible.
        //
        if (cff > one) {
            nfast=nfast+1;
            for (int i=nfast;i>=2;i--) {
                weight1[i-1]=weight1[i-1-1];
            }
            weight1[1-1] = zero;
        } else if (cff> zero) {
            wsum=one-cff;
            for (int i=nfast;i>=2;i--) {
                weight1[i-1]=wsum*weight1[i-1]+cff*weight1[i-1-1];
            }
            weight1[1-1]=wsum*weight1[1-1];
        } else if (cff < Real(-1.0)) {
            nfast=nfast-1;
            for (int i=1;i<=nfast;i++) {
                weight1[i-1]=weight1[i+1-1];
            }
            weight1[nfast+1-1] = zero;
        } else if (cff < zero) {
            wsum=one+cff;
            for (int i=1;i<=nfast-1;i++) {
                weight1[i-1]=wsum*weight1[i-1]-cff*weight1[i+1-1];
            }
            weight1[nfast-1]=wsum*weight1[nfast-1];
        }
    }

    //  Set SECONDARY weights assuming that backward Euler time step is used
    //  for free surface.  Notice that array weight2[i] is assumed to
    //  have all-zero status at entry in this segment of code.
    for(int j=1;j<=nfast;j++) {
        cff=weight1[j-1];
        for(int i=1;i<=j;i++) {
            weight2[i-1]=weight2[i-1]+cff;
        }
    }

    //
    //  Normalize both set of weights.
    //
    wsum = zero;
    cff  = zero;
    for(int i=1;i<=nfast;i++) {
        wsum=wsum+weight1[i-1];
        cff=cff+weight2[i-1];
    }

    wsum = one / wsum;
    cff  = one / cff;

    for(int i=1;i<=nfast;i++) {
        weight1[i-1]=wsum*weight1[i-1];
        weight2[i-1]=cff*weight2[i-1];
    }
}

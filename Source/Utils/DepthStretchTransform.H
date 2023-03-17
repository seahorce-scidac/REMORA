#ifndef STRETCH_H_
#define STRETCH_H_

#include <DataStruct.H>
#include <ROMSX.H>
#include <prob_common.H>

using namespace amrex;

void
ROMSX::set_2darrays (int lev)
{
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];
    std::unique_ptr<MultiFab>& mf_x_r = x_r[lev];
    std::unique_ptr<MultiFab>& mf_y_r = y_r[lev];

    for ( MFIter mfi(*(mf_x_r), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& x_r_arr = (mf_x_r)->array(mfi);
      Array4<Real> const& y_r_arr = (mf_y_r)->array(mfi);
      const Box& bx = mfi.growntilebox();
      const auto & geomdata = Geom(lev).data();
      int ncomp = 1;

      amrex::LoopConcurrentOnCpu(bx, ncomp,
      [=] (int i, int j, int k, int n)
      {
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();

        x_r_arr(i,j,0) = prob_lo[0] + (i + 0.5) * dx[0];
        y_r_arr(i,j,0) = prob_lo[1] + (j + 0.5) * dx[1];
        //        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

      });
    }
    //    x_r[lev]->FillBoundary(geom[lev].periodicity());
    //    y_r[lev]->FillBoundary(geom[lev].periodicity());
}

void ROMSX::set_vmix (int lev)
{
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];
    std::unique_ptr<MultiFab>& mf_Akv = Akv[lev];
    std::unique_ptr<MultiFab>& mf_z_w = z_w[lev];

    for ( MFIter mfi(*(mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv_arr = (mf_Akv)->array(mfi);
      Array4<Real> const& z_w_arr = (mf_z_w)->array(mfi);
      Box bx = mfi.tilebox();
      bx.grow(IntVect(2,2,0));
      const auto & geomdata = Geom(lev).data();
      int ncomp = 1;

      amrex::LoopConcurrentOnCpu(bx, ncomp,
      [=] (int i, int j, int k, int n)
      {

        Akv_arr(i,j,k) = 2.0e-03+8.0e-03*std::exp(z_w_arr(i,j,k)/150.0);

      });
    }
}

void
ROMSX::set_smflux(int lev, Real time)
{
    std::unique_ptr<MultiFab>& mf_sustr = sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = svstr[lev];
    auto geomdata = Geom(lev).data();
    bool NSPeriodic = geomdata.isPeriodic(1);
    bool EWPeriodic = geomdata.isPeriodic(0);
    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359;
    Real tdays=time/(24.0*60.0*60.0);
    amrex::Print()<<"Hacking in time offset for fixed dt=300"<<std::endl;
    //this is a hack because time is off by dt. this needs to be fixed for non-fixed dt
    Real dstart=-300.0/(24.0*60.0*60.0);
    Real rho0=parms.rho0;
    Real windamp;
    if(NSPeriodic) {
        mf_sustr->setVal(0.0);
    }
    else if(EWPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/rho0;
        else
            windamp=-0.1/rho0;
        mf_sustr->setVal(windamp);
    }
    if(NSPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1*sin(pi*(tdays-dstart)/4.0)/rho0;
        else
            windamp=-0.1/rho0;
        mf_svstr->setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr->setVal(0.0);
    }

}
void
ROMSX::set_depth (int lev)
{
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];
    std::unique_ptr<MultiFab>& mf_z_w = z_w[lev];
    std::unique_ptr<MultiFab>& mf_z_r = z_r[lev];
    std::unique_ptr<MultiFab>& mf_s_r = s_r[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = Hz[lev];
    std::unique_ptr<MultiFab>& mf_h  = hOfTheConfusingName[lev];
    mf_h->setVal(geom[lev].ProbHi(2));
    Real depth = geom[lev].ProbHi(2);
    const int Lm = Geom(lev).Domain().size()[0];
    const int Mm = Geom(lev).Domain().size()[1];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = Zt_avg1[lev];
    mf_Zt_avg1->setVal(0.0);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& z_w_arr = (mf_z_w)->array(mfi);
      Array4<Real> const& z_r_arr = (mf_z_r)->array(mfi);
      Array4<Real> const& s_r_arr = (mf_s_r)->array(mfi);
      Array4<Real> const& Hz_arr  = (mf_Hz)->array(mfi);
      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(2,2,0));
      Box gbx21 = bx;
      gbx21.grow(IntVect(2,2,1));
      Box gbx22 = bx;
      gbx22.grow(IntVect(2,2,2));
      Array4<Real> const& h_arr  = (mf_h)->array(mfi);
      Array4<Real> const& Zt_avg1_arr  = (mf_Zt_avg1)->array(mfi);
      const auto & geomdata = Geom(lev).data();

      int ncomp = 1;
      Box subdomain;
      if (lev == 0) {
          subdomain = geom[lev].Domain();
      } else {
          amrex::Abort("Geometry information needs to be updated for multilevel");
      }

      int nx = subdomain.length(0);
      int ny = subdomain.length(1);
      int nz = subdomain.length(2);

      auto N = nz; // Number of vertical "levels" aka, NZ
      //forcing tcline to be the same as probhi for now, one in DataStruct.H other in inputs
      Real hc=-min(geomdata.ProbHi(2),-solverChoice.tcline); // Do we need to enforce min here?
      bool NSPeriodic = geomdata.isPeriodic(1);
      bool EWPeriodic = geomdata.isPeriodic(0);
      amrex::LoopConcurrentOnCpu(gbx22, ncomp,
      [=] (int i, int j, int k, int n)
      {
#if 0
          Real val1, val2;
          int iFort = i+1;
          int jFort = j+1;
          if(NSPeriodic) {
              if (iFort<=Lm/2.0)
                  val1=iFort;
              else
                  val1=Lm+1-iFort;
              val2=-min(-geomdata.ProbLo(2),(84.5+66.526*std::tanh((val1-10.0)/7.0)));
              h_arr(i,j,0) = val2;
          }
          else if(EWPeriodic) {
              if (jFort<=Mm/2.0)
                  val1=jFort;
              else
                  val1=Mm+1-jFort;
              val2=-min(-geomdata.ProbLo(2),(84.5+66.526*std::tanh((val1-10.0)/7.0)));
              h_arr(i,j,0) = val2;
          }
#else
          h_arr(i,j,0) = geomdata.ProbLo(2);
#endif
          z_w_arr(i,j,k) = -h_arr(i,j,0);
      });
      // ROMS Transform 2
      amrex::LoopConcurrentOnCpu(gbx21, ncomp,
      [=] (int i, int j, int k, int n)
      {
        //        const Real z = prob_lo[2] + (k + 0.5) * dx[2];
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();
        //This is the z for the bottom of the cell this k corresponds to
        //   if we weren't stretching and transforming
        const Real z = prob_lo[2] + (k) * dx[2];
        //        h_arr(i,j,0) = -prob_lo[2]; // conceptually

        ////////////////////////////////////////////////////////////////////
        //ROMS Stretching 4
        // Move this block to it's own function for maintainability if needed
        // Information about the problem dimension would need to be added
        // This file would need a k dependent function to return the 
        // stretching scalars, or access to 4 vectors of length prob_length(2)
        /////////////////////////////////////////////////////////////////////
        Real ds = 1.0_rt / Real(N);
        
        Real cff_r, cff_w, cff1_r, cff1_w, cff2_r, cff2_w, Csur, Cbot;
        Real sc_r,sc_w,Cs_r,Cs_w;

        if(k==N) // end of array // pretend we're storing 0?
        {
        sc_w=0.0; //sc_w / hc
        Cs_w=0.0; //Cs_w
        }
        else if (k==0) // beginning of array
        {
        sc_w=-1.0; //sc_w / hc
        Cs_w=-1.0; //Cs_w
        }
        else
        {
          sc_w=ds*(k-N);
          if (solverChoice.theta_s > 0.0_rt)
            Csur=(1.0_rt-std::cosh(solverChoice.theta_s*sc_w))/
              (std::cosh(solverChoice.theta_s)-1.0_rt);
          else
            Csur=-sc_w*sc_w;
          if (solverChoice.theta_b > 0.0_rt)
          {
            Cbot=(std::exp(solverChoice.theta_b*Csur)-1.0_rt)/
              (1.0_rt-std::exp(-solverChoice.theta_b));
            Cs_w=Cbot;
          }
          else
            Cs_w=Csur;
        }

        cff_w=hc*sc_w;
        cff1_w=Cs_w;

        //cff_r => sc_r *hc
        //cff1_r => Cs_r
        //Don't do anything special for first/last index
        {
          sc_r=ds*(k-N+0.5_rt);
          if (solverChoice.theta_s > 0.0_rt)
            Csur=(1.0_rt-std::cosh(solverChoice.theta_s*sc_r))/
              (std::cosh(solverChoice.theta_s)-1.0_rt);
          else
            Csur=-sc_r*sc_r;
          if (solverChoice.theta_b > 0.0_rt)
          {
            Cbot=(std::exp(solverChoice.theta_b*Csur)-1.0_rt)/
              (1.0_rt-std::exp(-solverChoice.theta_b));
            Cs_r=Cbot;
          }
          else
            Cs_r=Csur;
        }
        if(i==0&&j==0&&k<N&&k>=0)
        s_r_arr(0,0,k) = sc_r;
        cff_r=hc*sc_r;
        cff1_r=Cs_r;

        ////////////////////////////////////////////////////////////////////
        Real hwater=-h_arr(i,j,0);
        if(k==0)
         {
          //extra guess added (maybe not actually defined in ROMS)
        Real hinv=1.0_rt/(hc+hwater);
         cff2_r=(cff_r+cff1_r*hwater)*hinv;
         z_w_arr(i,j,k-2) = hwater;
         z_w_arr(i,j,k-1)= -hwater;

         z_r_arr(i,j,k) = Zt_avg1_arr(i,j,0)+(Zt_avg1_arr(i,j,0)+hwater)*cff2_r;
         Hz_arr(i,j,k)=0.0;//z_w_arr(i,j,k)-z_w_arr(i,j,k-1);
	 if((i==0||i==-1)&&j==0&&k==0)
	      printf("%d %d %d %15.15g %15.15g %15.15g %15.15g %15.15g\n",i,j,k,Hz_arr(i,j,k),Hz_arr(i-1,j,k),z_w_arr(i,j,k),z_w_arr(i,j,k-1));
         }
        else
        {
        //Note, we are not supporting ICESHELF flag
          Real hinv=1.0_rt/(hc+hwater);
          cff2_r=(cff_r+cff1_r*hwater)*hinv;
          cff2_w=(cff_w+cff1_w*hwater)*hinv;

          z_w_arr(i,j,k-1)=Zt_avg1_arr(i,j,0)+(Zt_avg1_arr(i,j,0)+hwater)*cff2_w;
          z_r_arr(i,j,k)=Zt_avg1_arr(i,j,0)+(Zt_avg1_arr(i,j,0)+hwater)*cff2_r;
        //Note, we are not supporting ICESHELF flag
        }
      });

      amrex::LoopConcurrentOnCpu(bx, ncomp,
      [=] (int i, int j, int k, int n)
      {
          Hz_arr(i,j,k)=z_w_arr(i,j,k)-z_w_arr(i,j,k-1);
	 if((i==0||i==-1)&&j==0&&k==0)
	      printf("%d %d %d %15.15g %15.15g %15.15g %15.15g %15.15g\n",i,j,k,Hz_arr(i,j,k),Hz_arr(i-1,j,k),z_w_arr(i,j,k),z_w_arr(i,j,k-1));
      });
    } // mf
    mf_Hz->FillBoundary(geom[lev].periodicity());
    /*
    z_w[lev]->FillBoundary(geom[lev].periodicity());
    z_r[lev]->FillBoundary(geom[lev].periodicity());
    s_r[lev]->FillBoundary(geom[lev].periodicity());
    Hz[lev]->FillBoundary(geom[lev].periodicity());
    hOfTheConfusingName[lev]->FillBoundary(geom[lev].periodicity());
    Zt_avg1[lev]->FillBoundary(geom[lev].periodicity());
    */      
}

#endif

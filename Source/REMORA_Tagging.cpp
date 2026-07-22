#include <REMORA.H>
#include <REMORA_Derive.H>

using namespace amrex;

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in]  levc    level of refinement (0 is coarsest level)
 * @param[out] tags    array of tagged cells
 * @param[in]  time    current time
 * @param[in]  ngrow   number of grow cells
*/
void
REMORA::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    //
    // This mf must have ghost cells because we may take differences between adjacent values
    //
    std::unique_ptr<MultiFab> mf = std::make_unique<MultiFab>(grids[levc], dmap[levc], 1, 1);

    for (int j=0; j < ref_tags.size(); ++j)
    {
        if (ref_tags[j].Field() == "tracer" || ref_tags[j].Field() == "temp" ||
            ref_tags[j].Field() == "salt") {
            FillPatch(levc, time, *cons_new[levc], cons_new, BCVars::cons_bc, BdyVars::t,
                0,true,false);
        }
        // This allows dynamic refinement based on the value of the tracer
        if (ref_tags[j].Field() == "tracer")
        {
            MultiFab::Copy(*mf,*cons_new[levc],Tracer_comp,0,1,1);
        } else if (ref_tags[j].Field() == "temp") {
            MultiFab::Copy(*mf,*cons_new[levc],Temp_comp,0,1,1);
        } else if (ref_tags[j].Field() == "salt") {
            MultiFab::Copy(*mf,*cons_new[levc],Salt_comp,0,1,1);
        } else if (ref_tags[j].Field() == "x_velocity") {
            FillPatch(levc, time, *xvel_new[levc], xvel_new, xvel_bc(), BdyVars::u,0,true,true);
            MultiFab::Copy(*mf,*xvel_new[levc],0,0,1,1);
        } else if (ref_tags[j].Field() == "y_velocity") {
            FillPatch(levc, time, *yvel_new[levc], yvel_new, yvel_bc(), BdyVars::v,0,true,true);
            MultiFab::Copy(*mf,*yvel_new[levc],0,0,1,1);
        } else if (ref_tags[j].Field() == "z_velocity") {
            FillPatch(levc, time, *zvel_new[levc], zvel_new, zvel_bc(), BdyVars::null,0,true,true);
            MultiFab::Copy(*mf,*zvel_new[levc],0,0,1,1);
        } else if (ref_tags[j].Field() == "vorticity") {
            MultiFab mf_cc_vel(grids[levc],dmap[levc],3,1);
            average_face_to_cellcenter(mf_cc_vel,0,
                                       Array<const MultiFab*,3>{xvel_new[levc],
                                                                yvel_new[levc],
                                                                zvel_new[levc]});
            // Impose bc's at domain boundaries at all levels
            FillBdyCCVels(levc, mf_cc_vel);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                auto& dfab = (*mf)[mfi];
                auto& sfab = mf_cc_vel[mfi];
                auto pm = vec_pm[levc]->const_array(mfi);
                auto pn = vec_pn[levc]->const_array(mfi);
                auto maskr = vec_mskr[levc]->const_array(mfi);
                derived::remora_dervort(bx, dfab, 0, 1, sfab, pm, pn, maskr, Geom(levc), time, nullptr, levc);
            } // mfi

          mf->FillBoundary(geom[levc].periodicity());
          //
          // TODO: we may need to fill physical boundaries here before tagging criteria are imposed
          //

        } else if (ref_tags[j].Field() == "mask") {
            MultiFab::Copy(*mf,*vec_mskr3d[levc],0,0,1,IntVect(1,1,0));
#ifdef REMORA_USE_PARTICLES
        } else {
            //
            // This allows dynamic refinement based on the number of particles per cell
            //
            // Note that we must count all the particles in levels both at and above the current,
            //      since otherwise, e.g., if the particles are all at level 1, counting particles at
            //      level 0 will not trigger refinement when regridding so level 1 will disappear,
            //      then come back at the next regridding
            //
            const auto& particles_namelist( particleData.getNames() );
            mf->setVal(zero);
            for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++)
            {
                std::string tmp_string(particles_namelist[i]+"_count");
                IntVect rr = IntVect::TheUnitVector();
                if (ref_tags[j].Field() == tmp_string) {
                    for (int lev = levc; lev <= finest_level; lev++)
                    {
                        MultiFab temp_dat(grids[lev], dmap[lev], 1, 0); temp_dat.setVal(0);
                        particleData[particles_namelist[i]]->IncrementWithTotal(temp_dat, lev);

                        MultiFab temp_dat_crse(grids[levc], dmap[levc], 1, 0); temp_dat_crse.setVal(0);

                        if (lev == levc) {
                            MultiFab::Copy(*mf, temp_dat, 0, 0, 1, 0);
                        } else {
                            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                                rr[d] *= ref_ratio[levc][d];
                            }
                            average_down(temp_dat, temp_dat_crse, 0, 1, rr);
                            MultiFab::Add(*mf, temp_dat_crse, 0, 0, 1, 0);
                        }
                    }
                }
            }


#endif
        }

        ref_tags[j](tags,mf.get(),clearval,tagval,time,levc,geom[levc]);
    }
}

/**
 * Function to define the refinement criteria based on user input
*/
void
REMORA::refinement_criteria_setup ()
{
    if (max_level > 0)
    {
        ParmParse pp(pp_prefix);
        Vector<std::string> refinement_indicators;
        pp.queryarr("refinement_indicators",refinement_indicators,0,pp.countval("refinement_indicators"));

        for (int i=0; i<refinement_indicators.size(); ++i)
        {
            std::string ref_prefix = pp_prefix + "." + refinement_indicators[i];

            ParmParse ppr(ref_prefix);
            RealBox realbox;
            int lev_for_box;

            int num_real_lo      = ppr.countval("in_box_lo");
            int num_indx_lo      = ppr.countval("in_box_lo_indices");
            int num_indx_lo_crse = ppr.countval("in_box_lo_indices_crse");

            int num_real_hi      = ppr.countval("in_box_hi");
            int num_indx_hi      = ppr.countval("in_box_hi_indices");
            int num_indx_hi_crse = ppr.countval("in_box_hi_indices_crse");

            AMREX_ALWAYS_ASSERT( (num_real_lo      == num_real_hi)      && (num_real_lo      == 0 || num_real_lo      >= 2) );
            AMREX_ALWAYS_ASSERT( (num_indx_lo      == num_indx_hi)      && (num_indx_lo      == 0 || num_indx_lo      >= 2) );
            AMREX_ALWAYS_ASSERT( (num_indx_lo_crse == num_indx_hi_crse) && (num_indx_lo_crse == 0 || num_indx_lo_crse >= 2) );

            // Problem low and high (in real not index space) are the same at all levels
            if ( !((num_real_lo >= AMREX_SPACEDIM-1 && num_indx_lo == 0 && num_indx_lo_crse == 0) ||
                   (num_indx_lo >= AMREX_SPACEDIM-1 && num_real_lo == 0 && num_indx_lo_crse == 0) ||
                   (num_indx_lo ==              0   && num_real_lo == 0 && num_indx_lo_crse == 0) ||
                   (num_indx_lo_crse >= AMREX_SPACEDIM-1 && num_real_lo == 0 && num_indx_lo == 0)
                ) )
            {
                amrex::Abort("Must only specify box for refinement using real OR index space with fine/coarse grid indices");
            }

            if (num_real_lo > 0) {
                std::vector<Real> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box > 0 && lev_for_box <= max_level)
                {
                    ppr.getarr("in_box_lo",box_lo,0,2);
                    ppr.getarr("in_box_hi",box_hi,0,2);
                    box_lo[2] = geom[0].ProbLo(2);
                    box_hi[2] = geom[0].ProbHi(2);
                    realbox = RealBox(&(box_lo[0]),&(box_hi[0]));

                    amrex::Print() << "Reading " << realbox << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    const auto* dx  = geom[lev_for_box].CellSize();
                    const Real* plo = geom[lev_for_box].ProbLo();

                    int ilo = static_cast<int>((box_lo[0] - plo[0])/dx[0]);
                    int jlo = static_cast<int>((box_lo[1] - plo[1])/dx[1]);
                    int klo = static_cast<int>((box_lo[2] - plo[2])/dx[2]);
                    int ihi = static_cast<int>((box_hi[0] - plo[0])/dx[0]-1);
                    int jhi = static_cast<int>((box_hi[1] - plo[1])/dx[1]-1);
                    int khi = static_cast<int>((box_hi[2] - plo[2])/dx[2]-1);

                    Box bx_old(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));

                    int mod_ilo = ilo%ref_ratio[lev_for_box-1][0];
                    int mod_jlo = jlo%ref_ratio[lev_for_box-1][1];

                    int mod_ihi = (ihi+1)%ref_ratio[lev_for_box-1][0];
                    int mod_jhi = (jhi+1)%ref_ratio[lev_for_box-1][1];

                    if (mod_ilo != 0) {
                        ilo -= mod_ilo;
                    }
                    if (mod_jlo != 0) {
                        jlo -= mod_jlo;
                    }
                    if (mod_ihi != 0) {
                        ihi += ref_ratio[lev_for_box-1][0] - mod_ihi;
                    }
                    if (mod_jhi != 0) {
                        jhi += ref_ratio[lev_for_box-1][1] - mod_jhi;
                    }
                    Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));
                    if (mod_ilo !=0 || mod_jlo !=0 || mod_ihi != 0 || mod_jhi != 0) {
                        amrex::Print() << "Fine box on level " << lev_for_box << " adjusted from " << bx_old << " to " << bx << " to make it valid for refinement." << std::endl;
                    }
                    boxes_at_level[lev_for_box].push_back(bx);
                    amrex::Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev

            } else if (num_indx_lo > 0) {

                std::vector<int> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box > 0 && lev_for_box <= max_level)
                {
                    if (n_error_buf[0] != IntVect::TheZeroVector()) {
                        amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
                    }

                    ppr.getarr("in_box_lo_indices",box_lo,0,num_indx_lo);
                    ppr.getarr("in_box_hi_indices",box_hi,0,num_indx_hi);

                    if (num_indx_lo < AMREX_SPACEDIM) {
                        box_lo[2] = geom[lev_for_box].Domain().smallEnd(2);
                        box_hi[2] = geom[lev_for_box].Domain().bigEnd(2);
                    }

                    Box bx(IntVect(box_lo[0],box_lo[1],box_lo[2]),IntVect(box_hi[0],box_hi[1],box_hi[2]));
                    const Box& domain = geom[lev_for_box].Domain();

                    if (!domain.contains(bx)) {
                        amrex::Print() << "\n";
                        amrex::Print() << "Box specified       is " << bx << std::endl;
                        amrex::Print() << "But domain at level is " << domain << std::endl;
                        amrex::Error("Specified box doesn't fit in the domain");
                    }

                    const auto* dx  = geom[lev_for_box].CellSize();
                    const Real* plo = geom[lev_for_box].ProbLo();
                    realbox = RealBox(plo[0]+ box_lo[0]   *dx[0], plo[1]+ box_lo[1]   *dx[1], plo[2]+ box_lo[2]   *dx[2],
                                      plo[0]+(box_hi[0]+1)*dx[0], plo[1]+(box_hi[1]+1)*dx[1], plo[2]+(box_hi[2]+1)*dx[2]);

                    Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    if(box_lo[0]%ref_ratio[lev_for_box-1][0] != 0){
                        amrex::Print()<< "Requested ilo in x-direction : " << box_lo[0] << std::endl;
                        amrex::Print() << "ilo = " << box_lo[0] << " is not divisible by ref_ratio in x direction = " <<
                                          ref_ratio[lev_for_box-1][0] << std::endl;
                        amrex::Error("Adjust in_box_lo_indices in x-direction to be divisible by ref_ratio and try again");
                    }
                    if((box_hi[0]+1)%ref_ratio[lev_for_box-1][0] != 0){
                        amrex::Print()<< "Requested ihi in x-direction : " << box_hi[0] << std::endl;
                        amrex::Print() << "ihi+1 = " << box_hi[0]+1 << " is not divisible by ref_ratio in x direction = " <<
                                          ref_ratio[lev_for_box-1][0] << std::endl;
                        amrex::Error("Adjust in_box_hi_indices in x-direction to be divisible by ref_ratio and try again");
                    }
                     if(box_lo[1]%ref_ratio[lev_for_box-1][1] != 0){
                        amrex::Print()<< "Requested jlo in y-direction : " << box_lo[1] << std::endl;
                        amrex::Print() << "jlo = " << box_lo[1] << " is not divisible by ref_ratio in y direction = " <<
                                          ref_ratio[lev_for_box-1][1] << std::endl;
                        amrex::Error("Adjust in_box_lo_indices in y-direction to be divisible by ref_ratio and try again");
                    }
                    if((box_hi[1]+1)%ref_ratio[lev_for_box-1][1] != 0){
                        amrex::Print()<< "Requested jhi in y-direction : " << box_hi[1] << std::endl;
                        amrex::Print() << "jhi+1 = " << box_hi[1]+1 << " is not divisible by ref_ratio in y direction = " <<
                                          ref_ratio[lev_for_box-1][1] << std::endl;
                        amrex::Error("Adjust in_box_hi_indices in y-direction to be divisible by ref_ratio and try again");
                    }
                    if(box_lo[2]%ref_ratio[lev_for_box-1][2] != 0){
                        amrex::Print()<< "Requested klo in z-direction : " << box_lo[2] << std::endl;
                        amrex::Print() << "klo = " << box_lo[2] << " is not   divisible by ref_ratio in z direction = " <<
                                          ref_ratio[lev_for_box-1][2] << std::endl;
                        amrex::Error("Adjust in_box_lo_indices in z-direction to be divisible by ref_ratio and try again");
                    }
                    if((box_hi[2]+1)%ref_ratio[lev_for_box-1][2] != 0){
                        amrex::Print()<< "Requested khi in z-direction : " << box_hi[2] << std::endl;
                        amrex::Print() << "khi+1 = " << box_hi[2]+1 << " is not divisible by ref_ratio in z direction = " <<
                                          ref_ratio[lev_for_box-1][2] << std::endl;
                        amrex::Error("Adjust in_box_hi_indices in z-direction to be divisible by ref_ratio and try again");
                    }

                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev

            } else if (num_indx_lo_crse > 0) {

                std::vector<int> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box > 0 && lev_for_box <= max_level)
                {
                    if (n_error_buf[0] != IntVect::TheZeroVector()) {
                        amrex::Abort("Don't use n_error_buf > 0 when setting the box explicitly");
                    }

                    ppr.getarr("in_box_lo_indices_crse",box_lo,0,num_indx_lo_crse);
                    ppr.getarr("in_box_hi_indices_crse",box_hi,0,num_indx_hi_crse);

                    if (num_indx_lo_crse < AMREX_SPACEDIM) {
                        box_lo[2] = geom[lev_for_box-1].Domain().smallEnd(2);
                        box_hi[2] = geom[lev_for_box-1].Domain().bigEnd(2);
                    }

                    Box bx(IntVect(box_lo[0],box_lo[1],box_lo[2]),IntVect(box_hi[0],box_hi[1],box_hi[2]));

                    if (!geom[lev_for_box-1].Domain().contains(bx)) {
                        amrex::Print() << "\n";
                        amrex::Print() << "(Coarse) Box specified       is " << bx << std::endl;
                        amrex::Print() << "But (coarse) domain at level is " << geom[lev_for_box-1].Domain() << std::endl;
                        amrex::Error("Specified box doesn't fit in the domain");
                    }

                    bx.refine(ref_ratio[lev_for_box-1]);

                    const auto* dx  = geom[lev_for_box-1].CellSize();

                    const Real* plo = geom[lev_for_box].ProbLo();
                    realbox = RealBox(plo[0]+ box_lo[0]   *dx[0], plo[1]+ box_lo[1]   *dx[1], plo[2]+ box_lo[2]   *dx[2],
                                      plo[0]+(box_hi[0]+1)*dx[0], plo[1]+(box_hi[1]+1)*dx[1], plo[2]+(box_hi[2]+1)*dx[2]);

                    Print() << "Reading " << bx << " at level " << lev_for_box << std::endl;
                    num_boxes_at_level[lev_for_box] += 1;

                    boxes_at_level[lev_for_box].push_back(bx);
                    Print() << "Saving in 'boxes at level' as " << bx << std::endl;
                } // lev
            }

            AMRErrorTagInfo info;

            if (realbox.ok()) {
                info.SetRealBox(realbox);
            }
            if (ppr.countval("start_time") > 0) {
                Real ref_min_time; ppr.get("start_time",ref_min_time);
                info.SetMinTime(ref_min_time);
            }
            if (ppr.countval("end_time") > 0) {
                Real ref_max_time; ppr.get("end_time",ref_max_time);
                info.SetMaxTime(ref_max_time);
            }
            if (ppr.countval("max_level") > 0) {
                int ref_max_level; ppr.get("max_level",ref_max_level);
                info.SetMaxLevel(ref_max_level);
            }

            if (ppr.countval("value_greater")) {
                int num_val = ppr.countval("value_greater");
                Vector<Real> value(num_val);
                ppr.getarr("value_greater",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
            }
            else if (ppr.countval("value_less")) {
                int num_val = ppr.countval("value_less");
                Vector<Real> value(num_val);
                ppr.getarr("value_less",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
            }
            else if (ppr.countval("adjacent_difference_greater")) {
                int num_val = ppr.countval("adjacent_difference_greater");
                Vector<Real> value(num_val);
                ppr.getarr("adjacent_difference_greater",value,0,num_val);
                std::string field; ppr.get("field_name",field);
                ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
            }
            else if (realbox.ok())
            {
                ref_tags.push_back(AMRErrorTag(info));
            } else {
                Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
            }
        } // loop over criteria
        // Untag anywhere we have masks
        AMRErrorTagInfo info;
        info.SetDerefine(1);
        Real value = half;
        ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,"mask",info));
    } // if max_level > 0
}

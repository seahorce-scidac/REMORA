#include <REMORA.H>
#include <Derive.H>

using namespace amrex;

/**
 * Function to tag cells for refinement -- this overrides the pure virtual function in AmrCore
 *
 * @param[in] levc level of refinement (0 is coarsest level)
 * @param[out] tags array of tagged cells
 * @param[in] time current time
*/

void
REMORA::ErrorEst (int levc, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    bool tag_found = false;

    for (int j=0; j < ref_tags.size(); ++j)
    {
        std::unique_ptr<MultiFab> mf = std::make_unique<MultiFab>(grids[levc], dmap[levc], 1, 0);

        // This allows dynamic refinement based on the value of the scalar
        if (ref_tags[j].Field() == "scalar")
        {
            MultiFab::Copy(*mf,*cons_new[levc],Scalar_comp,0,1,0);
            tag_found = true;

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
            mf->setVal(0.0);
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
                    tag_found = true;
                }
            }


#endif
        }

        if (tag_found) {
            ref_tags[j](tags,mf.get(),clearval,tagval,time,levc,geom[levc]);
        }
    } // loop over j
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
            if (ppr.countval("in_box_lo")) {
                std::vector<Real> box_lo(3), box_hi(3);
                ppr.get("max_level",lev_for_box);
                if (lev_for_box <= max_level)
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
                    Box bx(IntVect(ilo,jlo,klo),IntVect(ihi,jhi,khi));
                    if ( (ilo%ref_ratio[lev_for_box-1][0] != 0) || ((ihi+1)%ref_ratio[lev_for_box-1][0] != 0) ||
                         (jlo%ref_ratio[lev_for_box-1][1] != 0) || ((jhi+1)%ref_ratio[lev_for_box-1][1] != 0) )
                         amrex::Error("Fine box is not legit with this ref_ratio");
                    boxes_at_level[lev_for_box].push_back(bx);
                    amrex::Print() << "Saving in 'boxes at level' as " << bx << std::endl;
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
    } // if max_level > 0
}


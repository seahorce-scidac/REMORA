#include <REMORA.H>
#include <REMORA_Derive.H>

using namespace amrex;

namespace {
/**
 * Copy mf's topmost and bottommost valid z-planes into its z-ghost planes.
 *
 * Several of the fields ErrorEst tags on live in MultiFabs with no ghost cells in z
 * (zvel_new, vec_mskr3d), and the vorticity branch only ever computes valid cells, so
 * mf's k = klo-1 and k = khi+1 planes would otherwise be read uninitialized by the GRAD
 * (adjacent_difference_greater) test, which unconditionally differences k+-1. There is
 * no data below k=0 or above k=N to copy, so impose a zero vertical gradient: GRAD then
 * sees no vertical difference at the surface and bottom cells instead of garbage.
 *
 * @param[inout] mf   single-component MultiFab with at least one ghost cell in z
 */
void
fill_z_ghost_planes (MultiFab& mf)
{
    AMREX_ALWAYS_ASSERT(mf.nComp() == 1 && mf.nGrowVect()[2] >= 1);

    // Not tiled: each iteration writes the two z-planes of a whole column, so tiling in
    // z would have several tiles writing the same cell.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const int klo = vbx.smallEnd(2);
        const int khi = vbx.bigEnd(2);

        // Laterally grown as well, so the ghost columns get their z-planes too --
        // GRAD reads i+-1 and j+-1 in those columns as well as k+-1.
        const Box& gbx = mfi.fabbox();
        auto arr = mf.array(mfi);

        ParallelFor(makeSlab(gbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            arr(i,j,klo-1) = arr(i,j,klo);
            arr(i,j,khi+1) = arr(i,j,khi);
        });
    }
}

/**
 * Fill mf's lateral ghost cells from the nearest valid cell of the same box.
 *
 * Call this *before* FillBoundary: FillBoundary reads only valid regions, so it
 * overwrites every ghost cell backed by a real same-level or periodic neighbor and
 * leaves only the physical-boundary and coarse/fine ghosts holding the zero-gradient
 * value written here. That keeps the GRAD test's i+-1 / j+-1 reads defined and
 * deterministic everywhere without tagging on a jump we invented.
 *
 * @param[inout] mf   single-component MultiFab
 */
void
fill_lateral_ghosts_zero_grad (MultiFab& mf)
{
    AMREX_ALWAYS_ASSERT(mf.nComp() == 1);
    const IntVect ng = mf.nGrowVect();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box  gbx = mfi.growntilebox(IntVect(ng[0],ng[1],0));
        auto arr = mf.array(mfi);

        const int ilo = vbx.smallEnd(0); const int ihi = vbx.bigEnd(0);
        const int jlo = vbx.smallEnd(1); const int jhi = vbx.bigEnd(1);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i < ilo || i > ihi || j < jlo || j > jhi) {
                arr(i,j,k) = arr(amrex::min(amrex::max(i,ilo),ihi),
                                 amrex::min(amrex::max(j,jlo),jhi), k);
            }
        });
    }
}
} // namespace

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

    // Any cell-centered tracer may drive refinement, named as it is elsewhere: "temp",
    // "salt", "tracer", "tracer_1", or a biology tracer such as "NO3".
    auto cons_comp_for_field = [this] (const std::string& field) {
        for (int icomp = 0; icomp < ncons; ++icomp) {
            if (cons_names[icomp] == field) { return icomp; }
        }
        return -1;
    };

    for (int j=0; j < ref_tags.size(); ++j)
    {
        const int cons_comp = cons_comp_for_field(ref_tags[j].Field());

        if (cons_comp >= 0) {
            FillPatch(levc, time, *cons_new[levc], cons_new, BCVars::cons_bc, BdyVars::t,
                0,true,false);
        }
        // This allows dynamic refinement based on the value of a tracer
        if (cons_comp >= 0)
        {
            MultiFab::Copy(*mf,*cons_new[levc],cons_comp,0,1,1);
        } else if (ref_tags[j].Field() == "x_velocity") {
            FillPatch(levc, time, *xvel_new[levc], xvel_new, xvel_bc(), BdyVars::u,0,true,true);
            MultiFab::Copy(*mf,*xvel_new[levc],0,0,1,1);
        } else if (ref_tags[j].Field() == "y_velocity") {
            FillPatch(levc, time, *yvel_new[levc], yvel_new, yvel_bc(), BdyVars::v,0,true,true);
            MultiFab::Copy(*mf,*yvel_new[levc],0,0,1,1);
        } else if (ref_tags[j].Field() == "z_velocity") {
            FillPatch(levc, time, *zvel_new[levc], zvel_new, zvel_bc(), BdyVars::null,0,true,true);
            // zvel_new has no ghost cells in z, so we can only ask the copy for lateral ones
            MultiFab::Copy(*mf,*zvel_new[levc],0,0,1,IntVect(1,1,0));
            fill_z_ghost_planes(*mf);
        } else if (ref_tags[j].Field() == "vorticity") {
            // Fill the ghost cells of the face-based velocities -- including at
            // coarse/fine boundaries, which is what FillPatch's FillPatchTwoLevels
            // interpolates -- since the cell-centered velocities in those ghost cells
            // are read when computing vorticity below
            FillPatch(levc, time, *xvel_new[levc], xvel_new, xvel_bc(), BdyVars::u,0,true,true);
            FillPatch(levc, time, *yvel_new[levc], yvel_new, yvel_bc(), BdyVars::v,0,true,true);
            FillPatch(levc, time, *zvel_new[levc], zvel_new, zvel_bc(), BdyVars::null,0,true,true);

            // Vorticity needs no vertical neighbor, and zvel_new has no z ghost cells
            // to average from, so only grow laterally
            MultiFab mf_cc_vel(grids[levc],dmap[levc],3,IntVect(1,1,0));
            average_face_to_cellcenter(mf_cc_vel,0,
                                       Array<const MultiFab*,3>{xvel_new[levc],
                                                                yvel_new[levc],
                                                                zvel_new[levc]},
                                       IntVect(1,1,0));
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

            // remora_dervort only writes valid cells, so fill mf's own ghosts before the
            // tagging criteria difference across them. The zero-gradient pass goes first;
            // FillBoundary then replaces every ghost that has a real same-level or periodic
            // neighbor, leaving only the physical-boundary and coarse/fine ghosts extrapolated.
            fill_lateral_ghosts_zero_grad(*mf);
            mf->FillBoundary(geom[levc].periodicity());
            fill_z_ghost_planes(*mf);

        } else if (ref_tags[j].Field() == "mask") {
            // vec_mskr3d has no z ghost cells, so this copy leaves mf's top and bottom
            // ghost planes unwritten while the GRAD test differences in z. The mask is
            // constant down a column, so the zero vertical gradient imposed here is exact.
            MultiFab::Copy(*mf,*vec_mskr3d[levc],0,0,1,IntVect(1,1,0));
            fill_z_ghost_planes(*mf);
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
            bool matched_particle_count = false;
            for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++)
            {
                std::string tmp_string(particles_namelist[i]+"_count");
                IntVect rr = IntVect::TheUnitVector();
                if (ref_tags[j].Field() == tmp_string) {
                    matched_particle_count = true;
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

            // A box-only indicator carries no field name and tags geometrically, so it
            // never reads mf; only a named field that matched nothing is an error.
            if (!matched_particle_count && !ref_tags[j].Field().empty()) {
                amrex::Abort("Unknown refinement field '" + ref_tags[j].Field() +
                             "'. Use a tracer name, x_velocity, y_velocity, z_velocity, "
                             "vorticity, mask, or <particle>_count.");
            }
#else
        } else if (!ref_tags[j].Field().empty()) {
            // mf is uninitialized until something writes it, so an unrecognized named
            // field would otherwise tag on garbage. A box-only indicator has no field
            // name, tags geometrically, and never reads mf.
            amrex::Abort("Unknown refinement field '" + ref_tags[j].Field() +
                         "'. Use a tracer name, x_velocity, y_velocity, z_velocity, "
                         "vorticity, or mask.");
#endif
        }

        ref_tags[j](tags,mf.get(),clearval,tagval,time,levc,geom[levc]);
    }

    // Promote any tagged cell to a full local z-column.
    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const& tag = tags.array(mfi);

        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);

        amrex::ParallelFor(makeSlab(bx, 2, 0),
        [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            bool refine_col = false;
            for (int k = klo; k <= khi; ++k) {
                refine_col = refine_col || (tag(i,j,k) != TagBox::CLEAR);
            }

            if (refine_col) {
                for (int k = klo; k <= khi; ++k) {
                    tag(i,j,k) = TagBox::SET;
                }
            }
        });
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
        {
            // Untag anywhere we have masks
            AMRErrorTagInfo info;
            info.SetDerefine(1);
            Real value = Real(0.5);
            ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,"mask",info));
        }
        {
            // Also untag at mask-water boundaries
            AMRErrorTagInfo info;
            info.SetDerefine(1);
            Real value = Real(0.5);
            ref_tags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,"mask",info));
        }
    } // if max_level > 0
}

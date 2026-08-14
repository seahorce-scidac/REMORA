#include "remora_bio_cppdefs.h"
!
!=======================================================================
!  REMORA Fennel Fortran bridge (Path A oracle).
!
!  This single translation unit contains, in order:
!
!    1. Minimal stand-ins for the ROMS modules that fennel.h and
!       fennel_mod.h reference (mod_kinds, mod_param, mod_scalars,
!       mod_iounits, mod_ncparam, mod_grid, mod_forces, mod_ocean,
!       mod_stepping, dateclock_mod).  These exist only to let the tracked
!       ROMS sources compile and run outside ROMS.  They are NOT faithful
!       reproductions of the ROMS modules; they carry exactly the symbols
!       the two included files touch, and nothing more.
!
!    2. MODULE mod_biology, whose body is the tracked copy
!       REMORA_fennel_mod.h, included verbatim, plus the two bridge
!       diagnostic emitters appended inside its CONTAINS section.
!
!    3. MODULE biology_mod, whose body is the tracked copy
!       REMORA_fennel.h, included verbatim.
!
!    4. MODULE remora_fennel_isohelper, the iso_c_binding wrapper.  All
!       adaptation logic -- unit handling, index-base translation, packing
!       REMORA's split old/new state into the ROMS three-slot tracer array,
!       and runtime validation of the compiled option set -- lives here, so
!       the ROMS sources stay unmodified apart from the tag emitters.
!
!  Everything is kept in one file because AMReX's GNUmake dependency scanner
!  resolves module dependencies for .f90/.F90 sources only; .F (fixed-form)
!  sources are compiled without inter-file module ordering.  A single
!  translation unit removes the ordering problem entirely.
!
!  Build:  USE_FENNEL_FORT=TRUE (GNUmake).  See Exec/Make.REMORA.
!=======================================================================
!
      MODULE mod_kinds
!
!  ROMS numeric kinds.  Values match ROMS mod_kinds.F.
!
      implicit none
      PUBLIC
!
      integer, parameter :: i8 = selected_int_kind(8)
      integer, parameter :: r4 = selected_real_kind(6,30)
      integer, parameter :: r8 = selected_real_kind(12,300)
      integer, parameter :: dp = selected_real_kind(12,300)
!
      END MODULE mod_kinds
!
!=======================================================================
!
      MODULE mod_param
!
!  Grid, tile and tracer dimensioning.  Ngrids is fixed at 1: the bridge
!  validates one REMORA level per call and has no ROMS nesting.
!
      USE mod_kinds
!
      implicit none
      PUBLIC
!
      integer, parameter :: Ngrids = 1
      integer, parameter :: iNLM = 1
!
!  Tracer counts.  NAT is temperature and salinity; REMORA has no passive,
!  cohesive or non-cohesive sediment tracers, so NPT/NCS/NNS are zero.
!  NBT is assigned by initialize_biology from the compiled option set.
!
      integer, parameter :: NAT = 2
      integer, parameter :: NPT = 0
      integer, parameter :: NCS = 0
      integer, parameter :: NNS = 0
!
      integer :: NBT = 0
      integer :: NDbio2d = 0
      integer :: NDbio3d = 0
!
      integer :: N(Ngrids)
      integer :: NT(Ngrids)
!
!  Dynamic-memory accounting.  fennel_mod.h increments Dmem(1) on every
!  allocation; the bridge never reports it.
!
      real(r8) :: Dmem(Ngrids) = 0.0_r8
!
!  Header-name bookkeeping referenced by SUBROUTINE biology.  In ROMS these
!  live in mod_scalars and mod_iounits; they are hosted here because
!  SUBROUTINE biology does not USE those modules.
!
      logical :: Lbiofile(4) = .TRUE.
      character (len=256) :: BIONAME(4) = ' '
!
!  Tile bounds.  The bridge always presents a single tile, index 0, whose
!  extents are the REMORA tilebox.  Indices are REMORA's own 0-based i and j,
!  so no horizontal index translation happens anywhere in the bridge.
!
      TYPE T_BOUNDS
        integer :: Istr(0:0)
        integer :: Iend(0:0)
        integer :: Jstr(0:0)
        integer :: Jend(0:0)
        integer :: LBi(0:0)
        integer :: UBi(0:0)
        integer :: LBj(0:0)
        integer :: UBj(0:0)
      END TYPE T_BOUNDS
!
      TYPE (T_BOUNDS) :: BOUNDS(Ngrids)
!
      END MODULE mod_param
!
!=======================================================================
!
      MODULE mod_scalars
!
!  Time-stepping and physical scalars used by biology_tile, plus the
!  bridge's runtime diagnostic state.
!
      USE mod_kinds
      USE mod_param, ONLY : Ngrids
!
      implicit none
      PUBLIC
!
      real(r8), parameter :: sec2day = 1.0_r8/86400.0_r8
!
      real(r8) :: dt(Ngrids) = 0.0_r8
      real(dp) :: tdays(Ngrids) = 0.0_dp
      integer  :: iic(Ngrids) = 0
!
!  Reference density and specific heat.  fennel.h uses these to convert
!  ROMS kinematic srflx (degC m/s) into W/m2 before applying PARfrac.
!
      real(r8) :: rho0 = 1025.0_r8
      real(r8) :: Cp = 3985.0_r8
!
!-----------------------------------------------------------------------
!  Bridge diagnostic state (not part of ROMS).
!
!  Two-level gating.  REMORA_USE_BIOLOGY_DIAG decides at compile time whether
!  the emitters exist at all (off by default, because the C++ half of the
!  contract sits inside a GPU device lambda); these variables then select
!  verbosity at run time, so a validation campaign never rebuilds to change
!  paths or levels.  Set by the isohelper on every call from the values the
!  C++ side parsed out of remora.biology_debug*.  They are declared
!  unconditionally so the bind(C) signature stays fixed regardless of the
!  flag.
!
!    bio_dbg_level  0 off, 1 target column only, 2 every column
!    bio_dbg_i/j    target column, REMORA 0-based indices
!-----------------------------------------------------------------------
!
      integer :: bio_dbg_level = 0
      integer :: bio_dbg_i = 0
      integer :: bio_dbg_j = 0
!
      END MODULE mod_scalars
!
!=======================================================================
!
      MODULE mod_iounits
!
      implicit none
      PUBLIC
!
      integer :: stdout = 6
!
      END MODULE mod_iounits
!
!=======================================================================
!
      MODULE mod_ncparam
!
!  Tracer metadata indices.  Only the two active-tracer indices are needed;
!  biology_tile uses them to read temperature and salinity out of t.
!
      implicit none
      PUBLIC
!
      integer, parameter :: itemp = 1
      integer, parameter :: isalt = 2
!
      END MODULE mod_ncparam
!
!=======================================================================
!
      MODULE mod_biology
!
!  ROMS mod_biology.  Body is the tracked copy of fennel_mod.h, included
!  verbatim; the bridge diagnostic emitters are appended to its CONTAINS
!  section below the include.
!
#include "REMORA_fennel_mod.h"
!
!-----------------------------------------------------------------------
!  Bridge diagnostic emitters (not part of ROMS).
!
!  Format contract is frozen in Source/Biology/Fortran/tag_map.md and must
!  stay byte-identical to remora_fennel_tag() on the C++ side:
!
!    FENNEL-FORT <TAG> iter=%3d i=%5d j=%5d k=%4d <VAR> <ES26.17E2>
!
!  k is printed in REMORA's 0-based convention (Fortran k minus one) so the
!  two logs diff without an index map.  Values are printed unmodified.
!-----------------------------------------------------------------------
!
      SUBROUTINE fennel_tag (tag, iter, Istr, Iend, j,
     &                       IminS, ImaxS, UBk, UBt, Bio)
!
!  Body is gated, not the twelve call sites in REMORA_fennel.h.  That file is
!  the tracked copy of ROMS fennel.h and its diff against the pristine
!  ../fennel.h is exactly 36 added lines with zero deletions; wrapping each
!  call site would double that and make the copy harder to audit.  The
!  residual cost of a no-op call per block seam per j is negligible, and the
!  Fortran path only exists in bridge builds anyway.
!
#ifdef REMORA_USE_BIOLOGY_DIAG
      USE mod_scalars, ONLY : bio_dbg_level, bio_dbg_i, bio_dbg_j
#endif
!
      character (len=*), intent(in) :: tag
      integer, intent(in) :: iter, Istr, Iend, j
      integer, intent(in) :: IminS, ImaxS, UBk, UBt
      real(r8), intent(in) :: Bio(IminS:ImaxS,UBk,UBt)
!
#ifdef REMORA_USE_BIOLOGY_DIAG
      integer :: i, k
!
      IF (bio_dbg_level.le.0) RETURN
      IF ((bio_dbg_level.eq.1).and.(j.ne.bio_dbg_j)) RETURN
!
      DO i=Istr,Iend
        IF ((bio_dbg_level.eq.1).and.(i.ne.bio_dbg_i)) CYCLE
        DO k=1,UBk
          CALL fennel_tag_line (tag, iter, i, j, k, 'NO3 ',
     &                          Bio(i,k,iNO3_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'NH4 ',
     &                          Bio(i,k,iNH4_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'CHLO',
     &                          Bio(i,k,iChlo))
          CALL fennel_tag_line (tag, iter, i, j, k, 'PHYT',
     &                          Bio(i,k,iPhyt))
          CALL fennel_tag_line (tag, iter, i, j, k, 'ZOOP',
     &                          Bio(i,k,iZoop))
          CALL fennel_tag_line (tag, iter, i, j, k, 'LDEN',
     &                          Bio(i,k,iLDeN))
          CALL fennel_tag_line (tag, iter, i, j, k, 'SDEN',
     &                          Bio(i,k,iSDeN))
#ifdef PO4
          CALL fennel_tag_line (tag, iter, i, j, k, 'PO4 ',
     &                          Bio(i,k,iPO4_))
#endif
#ifdef CARBON
          CALL fennel_tag_line (tag, iter, i, j, k, 'LDEC',
     &                          Bio(i,k,iLDeC))
          CALL fennel_tag_line (tag, iter, i, j, k, 'SDEC',
     &                          Bio(i,k,iSDeC))
          CALL fennel_tag_line (tag, iter, i, j, k, 'TIC ',
     &                          Bio(i,k,iTIC_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'TALK',
     &                          Bio(i,k,iTAlk))
#endif
#ifdef OXYGEN
          CALL fennel_tag_line (tag, iter, i, j, k, 'OXYG',
     &                          Bio(i,k,iOxyg))
#endif
#ifdef ODU
          CALL fennel_tag_line (tag, iter, i, j, k, 'ODU ',
     &                          Bio(i,k,iODU_))
#endif
        END DO
      END DO
#endif
!
      RETURN
      END SUBROUTINE fennel_tag
!
      SUBROUTINE fennel_tag_t (tag, iter, Istr, Iend, j,
     &                         LBi, UBi, LBj, UBj, UBk, UBt, nnew, t)
!
!  As fennel_tag, but reads the updated global tracer array rather than the
!  Bio scratch array.  Used for G11_UPDATE, the only tag emitted after the
!  Hz-weighted increment has been applied to t(nnew).
!
#ifdef REMORA_USE_BIOLOGY_DIAG
      USE mod_scalars, ONLY : bio_dbg_level, bio_dbg_i, bio_dbg_j
#endif
!
      character (len=*), intent(in) :: tag
      integer, intent(in) :: iter, Istr, Iend, j
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt, nnew
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
!
#ifdef REMORA_USE_BIOLOGY_DIAG
      integer :: i, k
!
      IF (bio_dbg_level.le.0) RETURN
      IF ((bio_dbg_level.eq.1).and.(j.ne.bio_dbg_j)) RETURN
!
      DO i=Istr,Iend
        IF ((bio_dbg_level.eq.1).and.(i.ne.bio_dbg_i)) CYCLE
        DO k=1,UBk
          CALL fennel_tag_line (tag, iter, i, j, k, 'NO3 ',
     &                          t(i,j,k,nnew,iNO3_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'NH4 ',
     &                          t(i,j,k,nnew,iNH4_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'CHLO',
     &                          t(i,j,k,nnew,iChlo))
          CALL fennel_tag_line (tag, iter, i, j, k, 'PHYT',
     &                          t(i,j,k,nnew,iPhyt))
          CALL fennel_tag_line (tag, iter, i, j, k, 'ZOOP',
     &                          t(i,j,k,nnew,iZoop))
          CALL fennel_tag_line (tag, iter, i, j, k, 'LDEN',
     &                          t(i,j,k,nnew,iLDeN))
          CALL fennel_tag_line (tag, iter, i, j, k, 'SDEN',
     &                          t(i,j,k,nnew,iSDeN))
#ifdef PO4
          CALL fennel_tag_line (tag, iter, i, j, k, 'PO4 ',
     &                          t(i,j,k,nnew,iPO4_))
#endif
#ifdef CARBON
          CALL fennel_tag_line (tag, iter, i, j, k, 'LDEC',
     &                          t(i,j,k,nnew,iLDeC))
          CALL fennel_tag_line (tag, iter, i, j, k, 'SDEC',
     &                          t(i,j,k,nnew,iSDeC))
          CALL fennel_tag_line (tag, iter, i, j, k, 'TIC ',
     &                          t(i,j,k,nnew,iTIC_))
          CALL fennel_tag_line (tag, iter, i, j, k, 'TALK',
     &                          t(i,j,k,nnew,iTAlk))
#endif
#ifdef OXYGEN
          CALL fennel_tag_line (tag, iter, i, j, k, 'OXYG',
     &                          t(i,j,k,nnew,iOxyg))
#endif
#ifdef ODU
          CALL fennel_tag_line (tag, iter, i, j, k, 'ODU ',
     &                          t(i,j,k,nnew,iODU_))
#endif
        END DO
      END DO
#endif
!
      RETURN
      END SUBROUTINE fennel_tag_t
!
      SUBROUTINE fennel_tag_line (tag, iter, i, j, k, var, val)
!
!  Single record.  k is converted from the ROMS 1-based bottom-up index to
!  REMORA's 0-based bottom-up index here, and nowhere else.
!
      character (len=*), intent(in) :: tag, var
      integer, intent(in) :: iter, i, j, k
      real(r8), intent(in) :: val
!
#ifdef REMORA_USE_BIOLOGY_DIAG
!  Re-pad the tag and variable names into fixed-width locals before writing.
!
!  This is load-bearing, not cosmetic.  AMReX compiles .F sources through
!  "cpp -E -traditional-cpp -P", which collapses runs of two or more spaces
!  even inside Fortran character literals.  The 10-character tag literals in
!  REMORA_fennel.h therefore reach the compiler at whatever length survives
!  that collapse ('G00_PRE   ' arrives as 'G00_PRE ', while 'G03_ZMETAB'
!  arrives intact).  Edit descriptor Aw right-justifies a string shorter than
!  w, so writing the dummy directly would place the tag at a column that
!  varies from tag to tag and never matches the C++ "%-10s".
!
!  Character assignment blank-pads on the right, so assigning into a
!  fixed-length local restores left justification without depending on any
!  literal surviving preprocessing.
!
      character (len=10) :: tagf
      character (len=4) :: varf
!
      tagf = tag
      varf = var
!
      write (6,10) tagf, iter, i, j, k-1, varf, val
  10  format ('FENNEL-FORT ',a10,' iter=',i3,' i=',i5,' j=',i5,
     &        ' k=',i4,' ',a4,' ',ES26.17E2)
#endif
!
      RETURN
      END SUBROUTINE fennel_tag_line
!
      END MODULE mod_biology
!
!=======================================================================
!
      MODULE mod_grid
!
!  Referenced only by SUBROUTINE biology, which the bridge never calls.
!
      USE mod_kinds
      USE mod_param, ONLY : Ngrids
!
      implicit none
      PUBLIC
!
      TYPE T_GRID
        real(r8), pointer :: rmask(:,:)
        real(r8), pointer :: Hz(:,:,:)
        real(r8), pointer :: z_r(:,:,:)
        real(r8), pointer :: z_w(:,:,:)
      END TYPE T_GRID
!
      TYPE (T_GRID) :: GRID(Ngrids)
!
      END MODULE mod_grid
!
!=======================================================================
!
      MODULE mod_forces
!
      USE mod_kinds
      USE mod_param, ONLY : Ngrids
!
      implicit none
      PUBLIC
!
      TYPE T_FORCES
        real(r8), pointer :: srflx(:,:)
        real(r8), pointer :: sustr(:,:)
        real(r8), pointer :: svstr(:,:)
        real(r8), pointer :: Uwind(:,:)
        real(r8), pointer :: Vwind(:,:)
      END TYPE T_FORCES
!
      TYPE (T_FORCES) :: FORCES(Ngrids)
!
      END MODULE mod_forces
!
!=======================================================================
!
      MODULE mod_ocean
!
      USE mod_kinds
      USE mod_param, ONLY : Ngrids
!
      implicit none
      PUBLIC
!
      TYPE T_OCEAN
        real(r8), pointer :: t(:,:,:,:,:)
        real(r8), pointer :: pH(:,:)
      END TYPE T_OCEAN
!
      TYPE (T_OCEAN) :: OCEAN(Ngrids)
!
      END MODULE mod_ocean
!
!=======================================================================
!
      MODULE mod_stepping
!
      USE mod_param, ONLY : Ngrids
!
      implicit none
      PUBLIC
!
      integer :: nstp(Ngrids) = 1
      integer :: nnew(Ngrids) = 2
!
      END MODULE mod_stepping
!
!=======================================================================
!
      MODULE dateclock_mod
!
!  Stub.  caldate is referenced by fennel.h only under PCO2AIR_SECULAR,
!  which the bridge does not enable.  Present so the tracked copy compiles
!  unmodified if that option is ever turned on for a separate lane.
!
      USE mod_kinds
!
      implicit none
      PUBLIC :: caldate
!
      CONTAINS
!
      SUBROUTINE caldate (CurrentTime, yy_i, yd_i, mm_i, dd_i,
     &                    h_i, m_i, s_r, yd_dp)
!
      real(dp), intent(in) :: CurrentTime
      integer, intent(out), optional :: yy_i, yd_i, mm_i, dd_i
      integer, intent(out), optional :: h_i, m_i
      real(r8), intent(out), optional :: s_r
      real(dp), intent(out), optional :: yd_dp
!
      IF (PRESENT(yy_i)) yy_i=1951
      IF (PRESENT(yd_i)) yd_i=1
      IF (PRESENT(mm_i)) mm_i=1
      IF (PRESENT(dd_i)) dd_i=1
      IF (PRESENT(h_i)) h_i=0
      IF (PRESENT(m_i)) m_i=0
      IF (PRESENT(s_r)) s_r=0.0_r8
      IF (PRESENT(yd_dp)) yd_dp=1.0_dp
!
      RETURN
      END SUBROUTINE caldate
!
      END MODULE dateclock_mod
!
!=======================================================================
!
      MODULE biology_mod
!
!  ROMS biology.F.  Body is the tracked copy of fennel.h, included verbatim
!  apart from the fennel_tag calls at the frozen block seams (tag_map.md).
!
      implicit none
!
      CONTAINS
!
#include "REMORA_fennel.h"
!
      END MODULE biology_mod
!
!=======================================================================
!
      MODULE remora_fennel_isohelper
!
!  iso_c_binding wrapper.  Everything that adapts REMORA's data model to the
!  ROMS calling convention lives here:
!
!    * REMORA splits the tracer state across two MultiFabs (cons_old at nstp,
!      cons_new at nnew); ROMS carries one array with a three-slot time index.
!      The wrapper packs slot 1 from t_nstp and slot 2 from t_nnew, calls with
!      nstp=1/nnew=2, and copies slot 2 back out.
!    * REMORA k runs 0..N bottom-up; ROMS k runs 1..N(ng) bottom-up.  The C++
!      packer applies the +1 shift when filling the buffers, so inside this
!      module every k is already ROMS-based.
!    * i and j are passed through as REMORA's own 0-based indices.  Fortran
!      lower bounds are set to match, so no horizontal translation occurs.
!    * The compiled ROMS option set is validated against the C++ runtime
!      toggles on every call.
!
      USE ISO_C_BINDING
      USE mod_kinds
      USE mod_param
      USE mod_scalars
      USE mod_biology
      USE mod_ncparam, ONLY : itemp, isalt
      USE biology_mod, ONLY : biology_tile
!
      implicit none
!
      CONTAINS
!
      subroutine fennel_bridge_advance_c (
     &    ilo, ihi, jlo, jhi, nz, nbt_in,
     &    istr_in, iend_in, jstr_in, jend_in,
     &    use_po4, use_carbon, use_oxygen, use_odu,
     &    use_denit, use_biosed, use_bulk,
     &    bioiter_in, dt_sec, rho0_in, cp_in,
     &    p_AttSW, p_AttChl, p_PARfrac, p_Vp0, p_I_thNH4, p_D_p5NH4,
     &    p_NitriR, p_K_NO3, p_K_NH4, p_K_PO4, p_K_Phy, p_Chl2C_m,
     &    p_ChlMin, p_PhyCN, p_R_P2N, p_PhyIP, p_PhyIS, p_PhyMin,
     &    p_PhyMR, p_ZooAE_N, p_ZooCN, p_ZooBM, p_ZooER, p_ZooGR,
     &    p_ZooMin, p_ZooMR, p_LDeRRN, p_LDeRRC, p_CoagR, p_SDeRRN,
     &    p_SDeRRC, p_wPhy, p_wLDet, p_wSDet, p_pCO2air,
     &    rmask, Hz, z_r, z_w, srflx, sustr, svstr, uwind, vwind, pH,
     &    t_nstp, t_nnew,
     &    dbg_level, dbg_i, dbg_j)
     &    bind(C, name="fennel_bridge_advance_c")
!
!  Storage bounds (ilo..ihi, jlo..jhi, nz) and active bounds
!  (istr_in..iend_in, jstr_in..jend_in) are passed separately.  The bridge
!  currently packs exactly the tilebox, so they coincide, but keeping them
!  distinct preserves the ROMS memory-vs-compute distinction and leaves room
!  for a haloed pack later without an interface change.
!
      integer(C_INT), VALUE, intent(in) :: ilo, ihi, jlo, jhi, nz
      integer(C_INT), VALUE, intent(in) :: nbt_in
      integer(C_INT), VALUE, intent(in) :: istr_in, iend_in
      integer(C_INT), VALUE, intent(in) :: jstr_in, jend_in
      integer(C_INT), VALUE, intent(in) :: use_po4, use_carbon
      integer(C_INT), VALUE, intent(in) :: use_oxygen, use_odu
      integer(C_INT), VALUE, intent(in) :: use_denit, use_biosed
      integer(C_INT), VALUE, intent(in) :: use_bulk
      integer(C_INT), VALUE, intent(in) :: bioiter_in
      real(C_DOUBLE), VALUE, intent(in) :: dt_sec, rho0_in, cp_in
      real(C_DOUBLE), VALUE, intent(in) :: p_AttSW, p_AttChl, p_PARfrac
      real(C_DOUBLE), VALUE, intent(in) :: p_Vp0, p_I_thNH4, p_D_p5NH4
      real(C_DOUBLE), VALUE, intent(in) :: p_NitriR, p_K_NO3, p_K_NH4
      real(C_DOUBLE), VALUE, intent(in) :: p_K_PO4, p_K_Phy, p_Chl2C_m
      real(C_DOUBLE), VALUE, intent(in) :: p_ChlMin, p_PhyCN, p_R_P2N
      real(C_DOUBLE), VALUE, intent(in) :: p_PhyIP, p_PhyIS, p_PhyMin
      real(C_DOUBLE), VALUE, intent(in) :: p_PhyMR, p_ZooAE_N, p_ZooCN
      real(C_DOUBLE), VALUE, intent(in) :: p_ZooBM, p_ZooER, p_ZooGR
      real(C_DOUBLE), VALUE, intent(in) :: p_ZooMin, p_ZooMR, p_LDeRRN
      real(C_DOUBLE), VALUE, intent(in) :: p_LDeRRC, p_CoagR, p_SDeRRN
      real(C_DOUBLE), VALUE, intent(in) :: p_SDeRRC, p_wPhy, p_wLDet
      real(C_DOUBLE), VALUE, intent(in) :: p_wSDet, p_pCO2air
      integer(C_INT), VALUE, intent(in) :: dbg_level, dbg_i, dbg_j
!
      real(C_DOUBLE), intent(in) :: rmask(ilo:ihi,jlo:jhi)
      real(C_DOUBLE), intent(in) :: Hz(ilo:ihi,jlo:jhi,nz)
      real(C_DOUBLE), intent(in) :: z_r(ilo:ihi,jlo:jhi,nz)
      real(C_DOUBLE), intent(in) :: z_w(ilo:ihi,jlo:jhi,0:nz)
      real(C_DOUBLE), intent(in) :: srflx(ilo:ihi,jlo:jhi)
!
!  fennel.h reads sustr(i+1,j) and svstr(i,j+1) when averaging the surface
!  stress to cell centres for the gas-exchange transfer velocity, so both
!  stress buffers are packed one cell wider in each direction.
!
      real(C_DOUBLE), intent(in) :: sustr(ilo:ihi+1,jlo:jhi+1)
      real(C_DOUBLE), intent(in) :: svstr(ilo:ihi+1,jlo:jhi+1)
      real(C_DOUBLE), intent(in) :: uwind(ilo:ihi,jlo:jhi)
      real(C_DOUBLE), intent(in) :: vwind(ilo:ihi,jlo:jhi)
      real(C_DOUBLE), intent(inout) :: pH(ilo:ihi,jlo:jhi)
      real(C_DOUBLE), intent(in) ::
     &   t_nstp(ilo:ihi,jlo:jhi,nz,nbt_in+2)
      real(C_DOUBLE), intent(inout) ::
     &   t_nnew(ilo:ihi,jlo:jhi,nz,nbt_in+2)
!
      integer, parameter :: ng = 1
      integer, parameter :: tile = 0
      integer, parameter :: kstp = 1
      integer, parameter :: knew = 2
!
      integer :: i, j, k, itrc, ntrc
!
      real(r8), allocatable :: t(:,:,:,:,:)
!
      ntrc=nbt_in+NAT
!
!-----------------------------------------------------------------------
!  Publish grid, time-step and tile geometry into the stub modules.
!-----------------------------------------------------------------------
!
      N(ng)=nz
      NT(ng)=ntrc
      dt(ng)=REAL(dt_sec,r8)
      rho0=REAL(rho0_in,r8)
      Cp=REAL(cp_in,r8)
!
      BOUNDS(ng)%Istr(tile)=istr_in
      BOUNDS(ng)%Iend(tile)=iend_in
      BOUNDS(ng)%Jstr(tile)=jstr_in
      BOUNDS(ng)%Jend(tile)=jend_in
      BOUNDS(ng)%LBi(tile)=ilo
      BOUNDS(ng)%UBi(tile)=ihi
      BOUNDS(ng)%LBj(tile)=jlo
      BOUNDS(ng)%UBj(tile)=jhi
!
      bio_dbg_level=dbg_level
      bio_dbg_i=dbg_i
      bio_dbg_j=dbg_j
!
!-----------------------------------------------------------------------
!  Allocate and index the biology parameter arrays, then validate the
!  compiled option set against the caller's runtime toggles.
!-----------------------------------------------------------------------
!
      CALL initialize_biology
!
      CALL fennel_check_config (nbt_in, use_po4, use_carbon,
     &                          use_oxygen, use_odu, use_denit,
     &                          use_biosed, use_bulk)
!
!-----------------------------------------------------------------------
!  Copy the parameter set.  Order here mirrors the FennelParameters struct
!  declaration order in REMORA_Biology.H so the two lists review side by
!  side.
!-----------------------------------------------------------------------
!
      BioIter(ng)=bioiter_in
      AttSW(ng)=REAL(p_AttSW,r8)
      AttChl(ng)=REAL(p_AttChl,r8)
      PARfrac(ng)=REAL(p_PARfrac,r8)
      Vp0(ng)=REAL(p_Vp0,r8)
      I_thNH4(ng)=REAL(p_I_thNH4,r8)
      D_p5NH4(ng)=REAL(p_D_p5NH4,r8)
      NitriR(ng)=REAL(p_NitriR,r8)
      K_NO3(ng)=REAL(p_K_NO3,r8)
      K_NH4(ng)=REAL(p_K_NH4,r8)
      K_PO4(ng)=REAL(p_K_PO4,r8)
      K_Phy(ng)=REAL(p_K_Phy,r8)
      Chl2C_m(ng)=REAL(p_Chl2C_m,r8)
      ChlMin(ng)=REAL(p_ChlMin,r8)
      PhyCN(ng)=REAL(p_PhyCN,r8)
      R_P2N(ng)=REAL(p_R_P2N,r8)
      PhyIP(ng)=REAL(p_PhyIP,r8)
      PhyIS(ng)=REAL(p_PhyIS,r8)
      PhyMin(ng)=REAL(p_PhyMin,r8)
      PhyMR(ng)=REAL(p_PhyMR,r8)
      ZooAE_N(ng)=REAL(p_ZooAE_N,r8)
      ZooCN(ng)=REAL(p_ZooCN,r8)
      ZooBM(ng)=REAL(p_ZooBM,r8)
      ZooER(ng)=REAL(p_ZooER,r8)
      ZooGR(ng)=REAL(p_ZooGR,r8)
      ZooMin(ng)=REAL(p_ZooMin,r8)
      ZooMR(ng)=REAL(p_ZooMR,r8)
      LDeRRN(ng)=REAL(p_LDeRRN,r8)
      LDeRRC(ng)=REAL(p_LDeRRC,r8)
      CoagR(ng)=REAL(p_CoagR,r8)
      SDeRRN(ng)=REAL(p_SDeRRN,r8)
      SDeRRC(ng)=REAL(p_SDeRRC,r8)
      wPhy(ng)=REAL(p_wPhy,r8)
      wLDet(ng)=REAL(p_wLDet,r8)
      wSDet(ng)=REAL(p_wSDet,r8)
      pCO2air(ng)=REAL(p_pCO2air,r8)
!
!-----------------------------------------------------------------------
!  Pack REMORA's split old/new state into the ROMS three-slot tracer array.
!  Slot 3 is unused by biology_tile but is allocated for shape fidelity.
!-----------------------------------------------------------------------
!
      allocate ( t(ilo:ihi,jlo:jhi,nz,3,ntrc) )
!
      DO itrc=1,ntrc
        DO k=1,nz
          DO j=jlo,jhi
            DO i=ilo,ihi
              t(i,j,k,kstp,itrc)=REAL(t_nstp(i,j,k,itrc),r8)
              t(i,j,k,knew,itrc)=REAL(t_nnew(i,j,k,itrc),r8)
              t(i,j,k,3,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Call the unmodified ROMS kernel.
!-----------------------------------------------------------------------
!
      CALL biology_tile (ng, tile,
     &                   ilo, ihi, jlo, jhi, nz, ntrc,
     &                   ilo, ihi, jlo, jhi,
     &                   kstp, knew,
#ifdef MASKING
     &                   rmask,
#endif
     &                   Hz, z_r, z_w, srflx,
#if defined CARBON || defined OXYGEN
# ifdef BULK_FLUXES
     &                   uwind, vwind,
# else
     &                   sustr, svstr,
# endif
#endif
#ifdef CARBON
     &                   pH,
#endif
     &                   t)
!
!-----------------------------------------------------------------------
!  Return the updated nnew slot.  Temperature and salinity are untouched by
!  biology_tile, so copying the whole array back is safe and keeps the
!  unpack symmetric with the pack.
!-----------------------------------------------------------------------
!
      DO itrc=1,ntrc
        DO k=1,nz
          DO j=jlo,jhi
            DO i=ilo,ihi
              t_nnew(i,j,k,itrc)=REAL(t(i,j,k,knew,itrc),C_DOUBLE)
            END DO
          END DO
        END DO
      END DO
!
      deallocate ( t )
!
      RETURN
      end subroutine fennel_bridge_advance_c
!
!=======================================================================
!
      SUBROUTINE fennel_check_config (nbt_in, use_po4, use_carbon,
     &                                use_oxygen, use_odu, use_denit,
     &                                use_biosed, use_bulk)
!
!  ROMS resolves biology options at compile time; REMORA resolves them at
!  run time.  A bridge built for one option set and run against a case
!  configured for another would produce a plausible-looking but meaningless
!  comparison, so every call checks the two agree and stops if they do not.
!
      integer, intent(in) :: nbt_in
      integer, intent(in) :: use_po4, use_carbon, use_oxygen, use_odu
      integer, intent(in) :: use_denit, use_biosed, use_bulk
!
      logical :: bad
!
      integer :: want_po4, want_carbon, want_oxygen, want_odu
      integer :: want_denit, want_biosed, want_bulk
!
      want_po4=0
      want_carbon=0
      want_oxygen=0
      want_odu=0
      want_denit=0
      want_biosed=0
      want_bulk=0
#ifdef PO4
      want_po4=1
#endif
#ifdef CARBON
      want_carbon=1
#endif
#ifdef OXYGEN
      want_oxygen=1
#endif
#ifdef ODU
      want_odu=1
#endif
#ifdef DENITRIFICATION
      want_denit=1
#endif
#ifdef BIO_SEDIMENT
      want_biosed=1
#endif
#ifdef BULK_FLUXES
      want_bulk=1
#endif
!
      bad=.FALSE.
      IF (use_po4.ne.want_po4) bad=.TRUE.
      IF (use_carbon.ne.want_carbon) bad=.TRUE.
      IF (use_oxygen.ne.want_oxygen) bad=.TRUE.
      IF (use_odu.ne.want_odu) bad=.TRUE.
      IF (use_denit.ne.want_denit) bad=.TRUE.
      IF (use_biosed.ne.want_biosed) bad=.TRUE.
      IF (use_bulk.ne.want_bulk) bad=.TRUE.
      IF (nbt_in.ne.NBT) bad=.TRUE.
!
      IF (bad) THEN
        write (6,10)
        CALL fennel_config_row ('po4', use_po4, want_po4)
        CALL fennel_config_row ('carbon', use_carbon, want_carbon)
        CALL fennel_config_row ('oxygen', use_oxygen, want_oxygen)
        CALL fennel_config_row ('odu', use_odu, want_odu)
        CALL fennel_config_row ('denitrification', use_denit,
     &                          want_denit)
        CALL fennel_config_row ('bio_sediment', use_biosed, want_biosed)
        CALL fennel_config_row ('bulk_fluxes', use_bulk, want_bulk)
        CALL fennel_config_row ('tracer count', nbt_in, NBT)
        write (6,30)
!
!  ERROR STOP, not STOP.  gfortran gives a character STOP code an exit status
!  of zero, so a plain STOP here would abort the run while still reporting
!  success to the shell, and any exit-status gate (CTest, a campaign runner)
!  would score the mismatch as a pass.
!
        ERROR STOP 1
      END IF
!
!  Format literals use nX edit descriptors rather than runs of spaces: the
!  .F preprocessing pass collapses two or more consecutive spaces inside
!  character literals, which silently destroys column alignment.
!
  10  format (/,' FENNEL BRIDGE: runtime option set does not match the',
     &          ' compiled ROMS CPP set.',/,
     &          ' Comparing these two paths would be invalid.',//,
     &        3x,'option',11x,'runtime',3x,'compiled')
  30  format (/,' Edit Source/Biology/Fortran/remora_bio_cppdefs.h to',
     &          ' match the case, then rebuild.',/)
!
      RETURN
      END SUBROUTINE fennel_check_config
!
      SUBROUTINE fennel_config_row (name, runtime_val, compiled_val)
!
!  One row of the option-mismatch table.  The name is re-padded into a
!  fixed-length local for the same reason the tag emitter does it: trailing
!  blanks in the caller's literals do not survive preprocessing, and Aw
!  right-justifies anything shorter than w.
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: runtime_val, compiled_val
!
      character (len=15) :: nm
!
      nm = name
      write (6,10) nm, runtime_val, compiled_val
  10  format (3x,a15,2x,i7,3x,i8)
!
      RETURN
      END SUBROUTINE fennel_config_row
!
      END MODULE remora_fennel_isohelper

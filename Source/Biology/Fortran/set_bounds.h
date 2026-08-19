!
!  REMORA bridge replacement for the ROMS "set_bounds.h" include.
!
!  ROMS ships a large set_bounds.h that derives every tile-bound alias
!  (IstrB, IstrP, IstrR, IstrT, IstrU, IstrM, and the J equivalents) from
!  BOUNDS(ng).  biology_tile uses only Istr, Iend, Jstr and Jend -- verified
!  by inspecting the routine body -- so only those four are declared here.
!  Declaring the unused aliases would compile, but leaving them out keeps the
!  stub surface honest about what the bridge actually depends on.
!
!  This file is included at the end of the declaration section of
!  biology_tile, so the declarations below precede the assignments, which are
!  the first executable statements of the routine.  That is the same ordering
!  ROMS itself relies on.
!
      integer :: Istr, Iend, Jstr, Jend
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)

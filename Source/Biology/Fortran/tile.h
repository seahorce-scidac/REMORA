!
!  REMORA bridge replacement for the ROMS "tile.h" include.
!
!  Included only by SUBROUTINE biology, the ROMS-level driver that the bridge
!  never calls -- the isohelper calls biology_tile directly, because REMORA
!  owns tiling through MFIter.  SUBROUTINE biology is nonetheless part of the
!  tracked fennel.h copy and must compile, so this provides the bound aliases
!  it references.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: IminS, ImaxS, JminS, JmaxS
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
      IminS=LBi
      ImaxS=UBi
      JminS=LBj
      JmaxS=UBj

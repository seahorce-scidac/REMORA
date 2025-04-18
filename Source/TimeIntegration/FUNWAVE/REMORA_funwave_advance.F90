
MODULE FUNWAVE_ADVANCE_MODULE

   use iso_c_binding

   implicit none

   contains

       subroutine funwave_advance_2d(salt, ims, ime, jms, jme, kms, kme)

           INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
           REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme), INTENT(INOUT) :: salt

       end subroutine funwave_advance_2d

END MODULE FUNWAVE_ADVANCE_MODULE

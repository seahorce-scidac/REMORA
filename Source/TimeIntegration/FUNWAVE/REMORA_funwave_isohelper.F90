
MODULE remora_funwave_isohelper
  use ISO_C_BINDING
  use FUNWAVE_ADVANCE_MODULE, ONLY: funwave_advance_2d
  implicit none

  contains

     subroutine funwave_advance_c(salt,ims,ime,jms,jme,kms,kme)  bind(C, name="funwave_advance_c")
         INTEGER(C_INT), VALUE, intent(in)                                   :: ims,ime,kms,kme,jms,jme
         REAL(C_DOUBLE), intent(inout), dimension(ims:ime, jms:jme, kms:kme) :: salt
         call funwave_advance_2d(salt,ims,ime,jms,jme,kms,kme)
     end subroutine funwave_advance_c

END MODULE remora_funwave_isohelper

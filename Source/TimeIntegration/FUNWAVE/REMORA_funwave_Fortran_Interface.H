#ifndef FUNWAVE_FORTRAN_INTERFACE_H
#define FUNWAVE_FORTRAN_INTERFACE_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * C interface to the FUNWAVE Fortran advance
 *
 * @note: All array dimensions must match between C++ and Fortran calls.
 * @note: Arrays must be allocated in C++ before passing to this routine.
 */
void funwave_advance_c(double* salt, int ims, int ime, int jms, int jme, int kms, int kme);

#ifdef __cplusplus
}
#endif

#endif // FUNWAVE_FORTRAN_INTERFACE_H

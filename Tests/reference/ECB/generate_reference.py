#!/usr/bin/env python3
"""Regenerate native ECB column fixtures using the pinned, clean ROMS sources.

Requires a ROMS git clone, NumPy, cpp and gfortran. No working-tree source is
read. Only the adapter's fixed-form continuation is repaired, as documented.
Run: python3 generate_reference.py --roms /path/to/roms
"""
import argparse
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile
import numpy as np

REVISION = '31cabd318c9379429c365c308b9e0f2ee166a6fc'
HERE = Path(__file__).resolve().parent
PARAMETERS = '''AttSW AttChl Chl2C_m ChlMin CoagR D_p5NH4 I_thNH4 K_NH4
K_NO3 K_PO4 K_Phy LDeRRN LDeRRC NitriR PARfrac PhyCN R_P2N PhyIS PhyMin
PhyMR SDeRRN SDeRRC RDeRRN RDeRRC Vp0 wLDet wPhy wSDet ZooAE_N ZooBM ZooCN
ZooER ZooGR ZooMin ZooMR ElDON EsDON deltN deltC gammaC slCexc pCO2air'''.split()
# Synthetic, explicitly specified test inputs, not a calibrated estuary setup.
VALUES = [0.04, 0.02486, 0.0535, 0.001, 0.015, 0.1, 0.0095, 6.7,
          2., 32., 2., 0.09, 0.09, 0.05, 0.43, 6.625, 0.0625, 0.08, 0.001,
          0.05, 0.09, 0.09, 0.01, 0.01, 0.55, 5., 0.1, 0.1, 0.75, 0.1,
          6.625, 0.1, 0.03, 0.001, 0.025, 0.2, 0.1, 0.15, 0.275, 0.15, 0.45, 400.]


def run(args, **kwargs):
    return subprocess.check_output(args, text=True, **kwargs)


def adapter(kernel, module, bulk):
    # Retain the entire reference parameter module and numerical subroutines.
    # Replace only ROMS's host driver and tile-bound include with a one-column
    # test driver. The original header misses a fixed-form continuation here.
    kernel = kernel[kernel.index('      SUBROUTINE ecb_tile'):]
    kernel = kernel.replace('                           bustr, bvstr,',
                            '     &                     bustr, bvstr,')
    kernel = kernel.replace('#include "set_bounds.h"',
                            '      Istr=0\n      Iend=0\n      Jstr=0\n      Jend=0')
    kernel = ('      MODULE biology_mod\n      USE mod_param\n'
              '      USE mod_biology\n      USE mod_ncparam\n'
              '      USE mod_scalars\n      USE dateclock_mod\n'
              '      implicit none\n      CONTAINS\n' + kernel)
    # These bounds normally arrive through set_bounds.h.
    kernel = kernel.replace('      integer, parameter :: Nsink = 6',
                            '      integer :: Istr,Iend,Jstr,Jend\n'
                            '      integer, parameter :: Nsink = 6')
    setup = '\n'.join(f'      {name}(1)=params({i+1})'
                      for i, name in enumerate(PARAMETERS))
    return f'''
#define PO4
#define TALK_NONCONSERV
#define RW14_OXYGEN_SC
#define RW14_CO2_SC
#define MASKING
{'#define BULK_FLUXES' if bulk else ''}
      MODULE mod_param
      implicit none
      integer,parameter :: r8=kind(1.d0),dp=kind(1.d0)
      integer,parameter :: Ngrids=1,NAT=2,NPT=0,NCS=0,NNS=0
      integer :: N(1),NT(1),NBT,NDbio2d,NDbio3d
      real(r8) :: Dmem(1)=0.d0
      END MODULE
      MODULE mod_scalars
      USE mod_param
      real(r8) :: dt(1),tdays(1),rho0,Cp=3985.d0
      real(r8),parameter :: sec2day=1.d0/86400.d0
      END MODULE
      MODULE mod_ncparam
      integer,parameter :: itemp=1,isalt=2
      END MODULE
      MODULE dateclock_mod
      USE mod_param
      CONTAINS
      SUBROUTINE caldate(t,yy_i,yd_dp)
      real(r8),intent(in) :: t
      integer,intent(out) :: yy_i
      real(r8),intent(out) :: yd_dp
      yy_i=2000
      yd_dp=1.d0
      END SUBROUTINE
      END MODULE
{module}
{kernel}
      SUBROUTINE reference(nz,iterations,seconds,params,old,new,
     &                     thickness,interfaces,shortwave,stress,
     &                     bottom,wind,density,mask) bind(C)
      USE iso_c_binding
      USE biology_mod
      USE mod_biology
      USE mod_scalars
      implicit none
      integer(c_int),value :: nz,iterations
      real(c_double),value :: seconds,shortwave,stress,bottom,wind
      real(c_double),value :: density,mask
      real(c_double),intent(in) :: params({len(PARAMETERS)})
      real(c_double),intent(in) :: old(nz,20),thickness(nz)
      real(c_double),intent(in) :: interfaces(0:nz)
      real(c_double),intent(inout) :: new(nz,20)
      real(r8) :: t(0:1,0:1,nz,3,20)
      real(r8) :: hz(0:1,0:1,nz),zr(0:1,0:1,nz)
      real(r8) :: zw(0:1,0:1,0:nz),srf(0:1,0:1)
      real(r8) :: su(0:1,0:1),sv(0:1,0:1),bu(0:1,0:1),bv(0:1,0:1)
      real(r8) :: rm(0:1,0:1),ph(0:1,0:1)
      integer :: i,j,k,c
      if (.not.allocated(idbio)) CALL initialize_biology
      N(1)=nz
      NT(1)=20
      dt(1)=seconds
      tdays(1)=0.d0
      rho0=density
      BioIter(1)=iterations
{setup}
      do c=1,20
        do k=1,nz
          t(:,:,k,1,c)=old(k,c)
          t(:,:,k,2,c)=new(k,c)
        enddo
      enddo
      do k=1,nz
        hz(:,:,k)=thickness(k)
        zr(:,:,k)=0.d0
      enddo
      do k=0,nz
        zw(:,:,k)=interfaces(k)
      enddo
      srf=shortwave
      su={'wind' if bulk else 'stress'}
      sv={'0.d0' if bulk else 'stress*0.5d0'}
      bu=bottom
      bv=bottom*0.5d0
      rm=mask
      ph=8.d0
      CALL ecb_tile(1,0,0,1,0,1,nz,20,0,1,0,1,1,2,rm,
     &              hz,zr,zw,srf,su,sv,bu,bv,ph,t)
      do c=1,20
        do k=1,nz
          new(k,c)=t(0,0,k,2,c)
        enddo
      enddo
      END SUBROUTINE
'''


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--roms', required=True)
    parser.add_argument('--compiler', default='gfortran')
    args = parser.parse_args()
    sources = [run(['git', '-C', args.roms, 'show',
                    REVISION + ':ROMS/Nonlinear/Biology/' + name])
               for name in ['ecb.h', 'ecb_mod.h']]
    for name, source in zip(['ecb.h', 'ecb_mod.h'], sources):
        print(name, hashlib.sha256(source.encode()).hexdigest())
    parameters = np.array(VALUES, dtype=np.float64)
    cases = []
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        for bulk in [False, True]:
            raw = adapter(*sources, bulk)
            preprocessed = subprocess.run(['cpp', '-P', '-traditional-cpp'],
                                         input=raw, text=True, capture_output=True, check=True).stdout
            # ROMS preprocessing normally converts its dual-form continuations.
            statements = []
            for line in preprocessed.splitlines():
                if not line.strip() or line[0] in '!cC*':
                    continue
                content = line.split('!')[0][6:].strip().rstrip('&').strip()
                if not content:
                    continue
                if len(line)>5 and line[5] not in ' 0':
                    statements[-1] += ' ' + content
                else:
                    statements.append(content)
            (work/'reference.f90').write_text('\n'.join(statements)+'\n')
            library = work/f'reference_{int(bulk)}.dylib'
            run([args.compiler, '-shared', '-fPIC', '-O0', '-fcheck=all',
                 '-ffree-line-length-none', 'reference.f90', '-o', str(library)], cwd=work)
            function = ctypes.CDLL(str(library)).reference
            array = np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS')
            function.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double,
                                 array, array, array, array, array] + [ctypes.c_double]*6
            for case in range(12):
                nz = [2,3,8,16,8,8,8,8,8,8,8,2][case]
                iterations = [1,3,1,2,3,1,2,3,1,0,1,1][case]
                seconds = [300.,3600.,60.,900.,1800.,300.,300.,300.,300.,300.,86400.,172800.][case]
                light = 0. if case in [1,4,10] else 150./(1025.*3985.)
                stress = 0.0001 if case%2 else 0.
                bottom = [0.,1.e-6,1.e-4][case%3]
                wind = 5. if case%2 else 0.
                mask = 0. if case==8 else (0.5 if case==7 else 1.)
                hz = np.linspace(1.,4.,nz)
                zw = np.concatenate(([-hz.sum()], -hz.sum()+np.cumsum(hz)))
                old = np.empty((nz,20),order='F')
                initial = [14.,30.,8.,0.5,1.,0.3,0.1,0.2,0.1,3.,0.5,
                           1.2,0.6,2100.,2300.,50.,200.,1.,40.,10000.]
                for c in range(20):
                    old[:,c] = initial[c]*(1.+0.03*np.arange(nz))
                if case==5:
                    old[:,16] = 0.01  # Low oxygen.
                if case==6:
                    old[:,2:13] *= 0.01
                    old[:,13] = 3100. # Preserve reference TIC constraints.
                new = np.asfortranarray(old*hz[:,None]+7.)
                expected = new.copy(order='F')
                function(nz,iterations,seconds,parameters,old,expected,hz,zw,
                         light,stress,bottom,wind,1025.,mask)
                assert np.isfinite(expected).all()
                cases.append((int(bulk),nz,iterations,seconds,light,stress,bottom,wind,mask,
                              hz,zw,old,new,expected))
    output = [str(len(PARAMETERS)), ' '.join(PARAMETERS),
              ' '.join(format(x,'.17e') for x in VALUES), str(len(cases))]
    for bulk,nz,iterations,seconds,light,stress,bottom,wind,mask,hz,zw,old,new,expected in cases:
        output.append(' '.join(map(str,[bulk,nz,iterations,seconds,light,stress,bottom,wind,mask])))
        for values in [hz,zw,old.ravel(order='F'),new.ravel(order='F'),expected.ravel(order='F')]:
            output.append(' '.join(format(x,'.17e') for x in values))
    (HERE/'columns.txt').write_text('\n'.join(output)+'\n')
    print(f'Wrote {len(cases)} ROMS column fixtures')


if __name__ == '__main__':
    main()

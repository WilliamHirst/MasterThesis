*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsmhcom.h                              ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsmhcom.h            ***
c----------------------------------------------------------------------c
c by Torsten Bringmann (troms@physto.se), 2010-01-23


* MH switches
      integer mhtype
      real*8 mheps
      common /MHswitches/ mheps, mhtype

* useful constants
      real*8 mpl, zeta(8)
      parameter (mpl=1.2209D19)
      data zeta/1D99,1.64493D0,1.20206D0,1.20206D0,1.03693D0,1.0173D0,
     -           1.00835D0,1.00408D0/

* needed for the integration of the Boltzmann equation
      real*8 resk(9),Tint
      integer nres
      common /MHresonances/ resk,Tint,nres

* degrees of freedom
      integer nfmax
      parameter (nfmax=1000)
      integer nf ! number of table entries in dof.dat
      integer khi,klo
      real*8 tmev(nfmax), sqrtgtab(nfmax), sqrtgttab(nfmax)
      common /MHdof/ tmev,sqrtgtab,sqrtgttab,nf,khi,klo

* QCD phase transition
      real*8 tqcdmin,tqcdmax
      parameter (tqcdmin=0.154,tqcdmax=4*tqcdmin) ! in GeV

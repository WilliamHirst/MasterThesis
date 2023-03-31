***********************************************************************
*** This routine calculates the capture rates in the Sun.
*** It does the same thing as dsntcapsunnum (i.e. dsntcapsunnumi),
*** except that it uses tabulted versions of the results intead
*** of performing a numerical integration every time. 
*** Inputs: mx - neutralino mass in GeV
***         sigsi - spin-independent capture rate in cm^2
***         type - type of distribution (same as in dsntcapsunnum)
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      real*8 function dsntcapsuntab(mx,sigsi,sigsd)

      implicit none
      real*8 mx,sigsi,sigsd,dsntctabget
      integer type
      include 'dshmcom.h'

      dsntcapsuntab=dsntctabget('su',1,mx)*(sigsi/1.d-40)
     &  *(rhox/0.3d0)
     &  +dsntctabget('su',2,mx)*(sigsd/1.d-40)*(rhox/0.3d0)

      return
      end

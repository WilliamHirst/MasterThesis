***********************************************************************
*** This routine calculates the capture rates in the Earth.
*** It does the same thing as dsntcapearthnum (i.e. dsntcapearthnumi),
*** except that it uses tabulted versions of the results intead
*** of performing a numerical integration every time. 
*** Inputs: mx - neutralino mass in GeV
***         sigsi - spin-independent capture rate in cm^2
***         type - type of distribution (same as in dsntcapearthnum)
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      real*8 function dsntcapearthtab(mx,sigsi)

      implicit none
      real*8 mx,sigsi,dsntctabget
      include 'dshmcom.h'

      dsntcapearthtab=dsntctabget('ea',1,mx)*(sigsi/1d-40)
     &  *(rhox/0.3d0)

      return
      end

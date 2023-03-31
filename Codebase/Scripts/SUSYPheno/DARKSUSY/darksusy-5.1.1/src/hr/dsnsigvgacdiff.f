**********************************************************************
*** subroutine dsnsigvgacdiff gives the number of contiuous gammas
*** per GeV at the energy ega (in GeV) times the annihilation cross
*** section.
*** The result given is the number
***     N_gamma * (sigma * v) / (10-29 cm^3 s^-1)
*** in units of GeV^-1.
***
*** author: joakim edsjo, edsjo@physto.se
*** date: 00-09-03
**********************************************************************

      subroutine dsnsigvgacdiff(ega,nsigvgacdiff)
      implicit none
      real*8 nsigvgacdiff,dshaloyield,yieldga,ega
      integer istat
      include 'dsprep.h'
      include 'dshacom.h'

c-----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsnsigvgacdiff:',
     &    ' dshasetup must be called',
     &    ' before any halo rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      yieldga=dshaloyield(ega,152,istat)
      nsigvgacdiff=yieldga*(hasv/1.0d-29)

      end






**********************************************************************
*** subroutine dsnsigvgaline gives the number of photons times
*** the annihilation cross section into gamma gamma and Z gamma
*** respectively. The result given is the dimensionless number
***     N_gamma * (sigma * v) / (10-29 cm^3 s^-1)
***
*** author: joakim edsjo, edsjo@physto.se
*** date: 00-09-03
**********************************************************************

      subroutine dsnsigvgaline(nsigvgaga,nsigvgaz)
      implicit none
      real*8 nsigvgaga,nsigvgaz,dssigmav
      include 'dsprep.h'
      include 'dshacom.h'

c-----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsnsigvgaline:',
     &    ' dshasetup must be called',
     &    ' before any halo rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

c... gamma-gamma
c      nsigvgaga=2.0d0*(dssigmav(28)/1.0d-29)
      nsigvgaga=2.0d0*(habr(28)*hasv/1.0d-29)
c... z gamma
c      nsigvgaz=(dssigmav(29)/1.0d-29)
      nsigvgaz=(habr(29)*hasv/1.0d-29)
      end






      subroutine dsmssmzero
c----------------------------------------------------------------------
c
c     Routine to zero all EW SUSY parameters in DarkSUSY. This is
c     done to make sure we don't inherit old values from previous
c     models when reading e.g. SLHA files.
c     Author: Joakim Edsjo, edsjo@physto.se
c     Date: 2008-07-01
c----------------------------------------------------------------------

      implicit none
      include 'dsmssm.h'

      integer i
      mu = 0.d0
      m2 = 0.d0
      ma = 0.d0
      tanbe = 0.d0
      m1=0.0d0
      m2=0.0d0
      m3=0.0d0
      do i=1,3
         mass2q(i) = 0.d0
         mass2u(i) = 0.d0
         mass2d(i) = 0.d0
         mass2l(i) = 0.d0
         mass2e(i) = 0.d0
         asofte(i) = 0.d0
         asoftu(i) = 0.d0
         asoftd(i) = 0.d0
      enddo

      return
      end

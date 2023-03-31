**********************************************************************
*** function dsepktdiff calculates the differential flux of
*** positrons for the energy egev as a result of
*** neutralino annihilation in the halo.
*** NOTE. This routine uses the Kamionkowski and Turner expressions.
*** from prd 43(1991)1774. Only included for comparison.
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-02-10
**********************************************************************

      real*8 function dsepktdiff(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsprep.h'

c----------------------------------------------------------------------
      real*8 egev,dsepktig,dsf_int2,eint(0:5),flx,
     &  dsepktig2,dsepktline
      integer i,eplint
      external dsepktig,dsepktig2
      parameter(eplint=2) ! way to integrate

c----------------------------------------------------------------------

c...setup common block for integration
      ehere=egev
      hasmooth=2

      if (egev.lt.hamwimp) then

c...continuum positrons
        if (eplint.eq.1.or.egev.eq.0.0d0) then
          eint(0)=egev
          do i=1,5
            eint(i)=(hamwimp-egev)/5.0d0*dble(i)+egev
          enddo
          flx=0.0d0
          do i=1,5
            flx=flx+dsf_int2(dsepktig,eint(i-1),eint(i),1.d-3)
          enddo
          dsepktdiff=flx
        elseif (eplint.eq.2) then
          flx=dsf_int2(dsepktig2,1.0d0/(hamwimp**2),
     &          1.0d0/(egev**2),1.d-3)
          dsepktdiff=flx
        endif

c...monochromatic positrons
        dsepktdiff=dsepktdiff+dsepktline(egev)
      else
        dsepktdiff=0.0d0
      endif
      hasmooth=0

      return
      end





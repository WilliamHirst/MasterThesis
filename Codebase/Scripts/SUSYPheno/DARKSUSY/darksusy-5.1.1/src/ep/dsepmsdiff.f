**********************************************************************
*** function dsepmsdiff calculates the differential flux of
*** positrons for the energy egev as a result of
*** neutralino annihilation in the halo.
*** NOTE. This routine uses the Moskaleno and Strong expressions
*** from PRD 60 (1999) 063003 for z_h=4 kpc, isothermal halo
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: edward baltz (eabaltz@alum.mit.edu), joakim edsjo
*** date: 10/18/2001
**********************************************************************

      real*8 function dsepmsdiff(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsprep.h'

c----------------------------------------------------------------------
      real*8 egev,dsepmsig,dsf_int2,eint(0:5),flx,
     &     dsepmsig2,dsepmsline
      integer i,eplint
      external dsepmsig,dsepmsig2
      parameter(eplint=2) ! way to integrate

c----------------------------------------------------------------------

c...setup common block for integration
      ehere=egev
      hasmooth=2

c...neglects reacceleration at maximum energies,
c...but that is neglible at interesting neutralino masses
      if (egev.lt.hamwimp) then

c...continuum positrons
         if (eplint.eq.1.or.egev.eq.0.0d0) then
            eint(0)=egev
            do i=1,5
               eint(i)=(hamwimp-egev)/5.0d0*dble(i)+egev
            enddo
            flx=0.0d0
            do i=1,5
               flx=flx+dsf_int2(dsepmsig,eint(i-1),eint(i),1.d-3)
            enddo
            dsepmsdiff=flx
         elseif (eplint.eq.2) then
            flx=dsf_int2(dsepmsig2,1.0d0/(hamwimp**2),
     &           1.0d0/(egev**2),1.d-3)
            dsepmsdiff=flx
         endif

c...monochromatic positrons
         dsepmsdiff=dsepmsdiff+dsepmsline(egev)
      else
         dsepmsdiff=0.0d0
      endif
      hasmooth=0
      
      return
      end

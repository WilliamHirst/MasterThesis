**********************************************************************
*** function dsepgalpropdiff calculates the differential flux of
*** positrons for the energy egev as a result of
*** neutralino annihilation in the halo.
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: edward baltz (eabaltz@alum.mit.edu), joakim edsjo
*** date: 4/28/2006
**********************************************************************

      real*8 function dsepgalpropdiff(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dsgalpropcom.h'

c----------------------------------------------------------------------
      real*8 egev,dsepgalpropig,dsf_int2,eint(0:5),flx,
     &     dsepgalpropig2,dsepgalpropline,emax
      integer i,eplint
      external dsepgalpropig,dsepgalpropig2
      parameter(eplint=2) ! way to integrate

c----------------------------------------------------------------------

c...setup common block for integration
      ehere=egev
      hasmooth=2

      emax=hamwimp+10.d0
c...neglects reacceleration at maximum energies,
c...but that is neglible at interesting neutralino masses
      if (egev.lt.emax) then

c...continuum positrons
         if (eplint.eq.1.or.egev.eq.0.0d0) then
            eint(0)=egev
            do i=1,5
               eint(i)=(emax-egev)/5.0d0*dble(i)+egev
            enddo
            flx=0.0d0
            do i=1,5
               flx=flx+dsf_int2(dsepgalpropig,eint(i-1),eint(i),1.d-3)
            enddo
            dsepgalpropdiff=flx
         elseif (eplint.eq.2) then
            flx=dsf_int2(dsepgalpropig2,
     +           1.0d0/(emax**2),1.0d0/(egev**2),1.d-3)
            dsepgalpropdiff=flx
         endif

c...monochromatic positrons
         dsepgalpropdiff=dsepgalpropdiff+dsepgalpropline(egev)
      else
         dsepgalpropdiff=0.0d0
      endif
c...  galprop Green's function are in MeV
      dsepgalpropdiff=dsepgalpropdiff*1.d-3
      hasmooth=0
      
      return
      end

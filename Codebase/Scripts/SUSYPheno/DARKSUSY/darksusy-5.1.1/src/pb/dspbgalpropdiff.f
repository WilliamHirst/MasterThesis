**********************************************************************
*** function dspbgalpropdiff calculates the differential flux of
*** antiprotons for the energy egev as a result of
*** neutralino annihilation in the halo.
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: edward baltz (eabaltz@alum.mit.edu), joakim edsjo
*** date: 4/28/2006
**********************************************************************

      real*8 function dspbgalpropdiff(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dsgalpropcom.h'

c----------------------------------------------------------------------
      real*8 egev,dspbgalpropig,dsf_int2,eint(0:5),flx,
     &     dspbgalpropig2,emax
      integer i,eplint
      external dspbgalpropig,dspbgalpropig2
      parameter(eplint=2) ! way to integrate

c----------------------------------------------------------------------

      if (.not.gpgfread) then
         write (6,*)
     +        'You must call dsgalprop_gettable(0) to use GALPROP'
         stop
      endif

c...setup common block for integration
      ehere=egev
      hasmooth=2

      emax=mx+10.d0
c...neglects reacceleration at maximum energies,
c...but that is neglible at interesting neutralino masses
      if (egev.lt.emax) then

c...continuum antiprotons
         if (eplint.eq.1.or.egev.eq.0.0d0) then
            eint(0)=egev
            do i=1,5
               eint(i)=(emax-egev)/5.0d0*dble(i)+egev
            enddo
            flx=0.0d0
            do i=1,5
               flx=flx+dsf_int2(dspbgalpropig,eint(i-1),eint(i),1.d-3)
            enddo
            dspbgalpropdiff=flx
         elseif (eplint.eq.2) then
            flx=dsf_int2(dspbgalpropig2,
     +           1.0d0/(emax**2),1.0d0/(egev**2),1.d-3)
            dspbgalpropdiff=flx
         endif
      else
         dspbgalpropdiff=0.0d0
      endif
c...  galprop Green's function are in MeV
      dspbgalpropdiff=dspbgalpropdiff*1.d-3
      hasmooth=0
      
      return
      end

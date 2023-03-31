**********************************************************************
*** wimpexample. This program is an example how to use the halo
*** and neutrino telescope routines without relying on a specific
*** SUSY model, i.e. instead using a general WIMP model.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: December 5, 2008
**********************************************************************

      program wimpexample
      implicit none

      real*8 dshrpbardiff,dshrdbardiff,dsepdiff,dsepspecm,
     &  dshrgacontdiff
      real*8 pbflux,ekin,rateea,ratesu,phidb,phiep,phigac
      real*8 suy,mww,suyint,enu,theta,dth,thmax,suytab
      integer unphys,warning,istat,i,n,yk
      real*8 dswayield

c...Initialize DarkSUSY
      call dsinit

c...Set up halo model and diffusion model (not necessary for defaults)
      call dshmset('default')
      call dspbset('default')

c...Now set up the routines for a WIMP. This is usually done with dshasetup,
c...so we have here created a copy of that routine that we now set
c...up the model with

      call dshasetup_private  ! routine is further down this code

c...Do the same for neutrino telescope in case we need it
      mww=100.d0
      call dswasetup_private(mww)  ! routine is further down this code

c...Gamma flux at 1 GeV for a J*delta-Omega of .1 sr
      phigac=dshrgacontdiff(1.d0,0.1d0,istat)
      write(*,*) '  phigac = ',phigac,' ph/(cm^2 s GeV)'


c...Now calculate the antiproton flux
      ekin=1.d0 ! kinetic energy
      pbflux=dshrpbardiff(ekin,1,4) ! use tables, faster after slower startup
      write(*,*) '  solar modulated pbar flux = ',pbflux,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

c...Antideutron flux
      phidb=dshrdbardiff(1.d0,1,4)   ! use tables, faster after slower startup
      write(*,*) '  solar modulated dbar flux at 1.00 GeV = ',phidb,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

c...Positrons
c      phiep=dsepspecm(1.0d0,4)  ! use tables
c      write (*,*) '  phiep=',phiep,
c     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

      phiep=dsepdiff(1.0d0,1,4) ! use tables
      write (*,*) '  phiep=',phiep,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'


      call dsntrates(1.d0,30.d0,3,3,rateea,ratesu,istat)

      write(*,*) '  Flux from the Earth = ',rateea,
     &  ' km^-2 yr^-1'
      write(*,*) '  Flux from the Sun =   ',ratesu,
     &  ' km^-2 yr^-1'

      enu=10.d0
      n=1000
      thmax=30.0d0
      dth=thmax/n
      suyint=0.d0
      yk=13
      do i=1,n
         theta=(dble(i)-0.5d0)*dth
         suy=dswayield(mww,enu,theta,'su',2,yk,istat)
         suyint=suyint+suy*dth
         suytab=dswayield(mww,enu,theta+dth*0.5d0,'su',3,yk,istat)
         write(47,'(1x,i5,4(1x,e12.6))') i,theta,suy,suyint,suytab
      enddo

      write(*,*) 'Solar yields: '
      write(*,*) 'Direct integration: suyint = ',suyint,' / GeV'
      suy=dswayield(mww,enu,thmax,'su',3,yk,istat)
      write(*,*) 'Tabulated yields: ',
     &  suy,' / GeV'

      suy=dswayield(mww,enu,thmax,'su',1,13,istat)
      write(*,*) 'Int (2D) yields (just for test): ',
     &  suy,' / GeV'

      end
      

**********************************************************************
*** SUBROUTINES AND FUNCTIONS
**********************************************************************

*****************************************************************************
***   subroutine dshasetup prepares the common blocks with annihilation
***   channel branching ratios and Higgs decay widths for halo yield
***   calculations.
***   Note: This routine needs to be called for each each model before
***   dshaloyield is called.
***   This routine is the interface between SUSY and the halo annihilation
***   routines.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 08-01-15, modified 08-12-05 to show how to do things
*** for general WIMP models
*****************************************************************************

      subroutine dshasetup_private
      implicit none
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      real*8 dssigmav
      integer i,j

c----------------------------------------------- set-up common variables

c...In this example, we will set up a general WIMP of mass 200 GeV
c...annihilating to W+ W- only and with an annihilation cross section
c...of 3e-26 cm^3/s.


      dshasetupcalled=.true.

c...make sure we computed the branching ratios
      hasv=3.0d-26 ! annihilation cross section

c...Now we need to set up the branching fraction, I set all to zero first
c...and then rest those we want to include
c...Transfer to Halo annihilation common blocks
c...For a description of the channels, see dshayield.f
      do i=1,29
         habr(i)=0.d0
      enddo
      habr(13)=1.d0 ! W+W-, see dshayield.f for an explaination

c...If we have scalars (Higgses) as final states, we need to set up
c...the Higgs decay widths as well. Here I just put them to zero (the 
c...values will actually don't matter unless we have Higgses produced).

c...Transfer Higgs widths
c...Neutral Higgses
      do i=1,3
         do j=1,29
            has0br(j,i)=0.d0 ! decay branching fraction
         enddo
      enddo

c...Masses of Higgses
      do i=1,3
         has0m(i)=0.d0
      enddo

c...Charged Higgses
      do j=1,15
         hascbr(j)=0.d0 ! the same for charged Higgses
      enddo

      hascm=0.d0 ! mass of H+-

c...Internal Bremsstrahlung
c...This is the only model-dependent part that the halo routines need
c...to know about. For SUSY models, set to SUSY. For other models, set
c...to 'none' not to include FSR. You can also provide your own routines,
c...see dshaib for more info on viable options.
c      haib='susy'
      haib='none' ! for non-SUSY models, 


c...Mass of WIMP
      hamwimp=200.0d0 ! WIMP mass in GeV

      return

      end


*****************************************************************************
***   subroutine dswasetup prepares the common blocks with annihilation
***   channel branching ratios and Higgs decay widths for WIMP annihilation
***   yield calculations from the Sun and the Earth.
***   Note: This routine needs to be called for each each model before
***   dswayield is called.
***   This routine is the interface between SUSY and the halo annihilation
***   routines.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 08-04-02
*****************************************************************************

      subroutine dswasetup_private(mww)
      implicit none
      include 'dswacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

      real*8 dssigmav,tmp1,tmp2,mww
      integer i,j,kh(4)

c----------------------------------------------- set-up common variables

c...In this example, we will set up a general WIMP of mass 200 GeV
c...annihilating to W+ W- only and with an annihilation cross section
c...of 3e-26 cm^3/s and spin-dependent and spin-independent scattering
c...cross sections of 1 pb.

      dswasetupcalled=.true.

c...make sure we computed the branching ratios
      wasv=3.d-26 ! annihilation cross section

c...Transfer to WIMP annihilation common blocks
c...For a description of the channels, see dswayieldone.f
      do i=1,29
         wabr(i)=0.d0
      enddo

c      wabr(13)=1.d0 ! W+W-, see dswayieldone.f for an explaination
      wabr(25)=1.d0 ! b b-bar, see dswayieldone.f for an explaination

c...Transfer Higgs widths
c...Neutral Higgses
      do i=1,3
         do j=1,29
            was0br(j,i)=0.d0
         enddo
      enddo

      do i=1,3
         was0m(i)=0.d0
      enddo

c...Charged Higgses
      do j=1,15
         wascbr(j)=0.d0
      enddo

      wascm=0.d0

c...Mass of WIMP
      wamwimp=mww

c...Scattering cross sections
      wasigsip=1.d-36*1.d-5 ! cm^2 SI scattering cross section
      wasigsdp=1.d-36*1.d-5 ! cm^2 SD scattering cross section

      return
      end

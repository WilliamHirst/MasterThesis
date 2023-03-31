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

      subroutine dswasetup
      implicit none
      include 'dsmssm.h'
      include 'dswacom.h'
      include 'dshrcom.h'
      include 'dsprep.h'

      real*8 dssigmav,tmp1,tmp2,br
      integer i,j,kh(4)

c----------------------------------------------- set-up common variables

c...Set up Higgses
      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3
      kh(4)=khc


      dswasetupcalled=.true.

c...make sure we computed the branching ratios
      wasv=dssigmav(0)

c...Transfer to WIMP annihilation common blocks
c...For a description of the channels, see dswayieldone.f
      do i=1,29
         wabr(i)=sigv(i)/wasv
      enddo

c...Transfer Higgs widths
c...Neutral Higgses
      do i=1,3
         do j=1,29
            br=hdwidth(j,i)/width(kh(i))
            if (br.gt.wasbrmin) then ! to gain some speed
               was0br(j,i)=br
            else
               was0br(j,i)=0.d0
            endif
         enddo
      enddo

      do i=1,3
         was0m(i)=mass(kh(i))
      enddo

c...Charged Higgses
      do j=1,15
         br=hdwidth(j,4)/width(kh(4))
         if (br.gt.wasbrmin) then
            wascbr(j)=br
         else
            wascbr(j)=0.d0
         endif
      enddo

      wascm=mass(kh(4))

c...Mass of WIMP
      wamwimp=mass(kn(1))

c...Scattering cross sections
      call dsddneunuc(wasigsip,tmp1,wasigsdp,tmp2)

      end




















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
*** Date: 08-01-15
*****************************************************************************

      subroutine dshasetup
      implicit none
      include 'dsmssm.h'
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      real*8 dssigmav,br
      integer i,j,kh(4)

c----------------------------------------------- set-up common variables

c...Set up Higgses
      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3
      kh(4)=khc


      dshasetupcalled=.true.

c...make sure we computed the branching ratios
      hasv=dssigmav(0)

c...Transfer to Halo annihilation common blocks
c...For a description of the channels, see dshayield.f
      do i=1,29
         habr(i)=sigv(i)/hasv
      enddo

c...Transfer Higgs widths
c...Neutral Higgses
      do i=1,3
         do j=1,29
            br=hdwidth(j,i)/width(kh(i))
            if (br.gt.hasbrmin) then
               has0br(j,i)=br
            else
               has0br(j,i)=0.d0
            endif
         enddo
      enddo

      do i=1,3
         has0m(i)=mass(kh(i))
      enddo

c...Charged Higgses
      do j=1,15
         br=hdwidth(j,4)/width(kh(4))
         if (br.gt.hasbrmin) then
            hascbr(j)=br
         else
            hascbr(j)=0.d0
         endif
      enddo

      hascm=mass(kh(4))

c...Internal Bremsstrahlung
c...See dshaib for more info on viable options.
      haib='susy'


c...Mass of WIMP
      hamwimp=mass(kn(1))

      end




















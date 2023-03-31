*****************************************************************************
***   subroutine dsprep calculates different frequently used cross
***   sections, annihilation branching ratios and sets different common
***   block variable when a new model has been defined. this routine
***   hasto be called before the relic density or any rates are
***   calculated and is called from dssusy.
***   
*** author: joakim edsjo  edsjo@physics.berkeley.edu
*** date: 98-02-27
*** modified: 99-02-23 pg, 00-08-09 je, 08-01-15 JE (moved halo yield
*** setup to dshasetup), 10-08-04 pg (removed mh)
*****************************************************************************

      subroutine dsprep
      implicit none
      include 'dsprep.h'
      include 'dsidtag.h'


c------------------------ functions ------------------------------------

      real*8 dsandwdcosnn

c------------------------ variables ------------------------------------

c      integer i,j,kh(4)

c----------------------------------------------- set-up common variables

c...flags for subroutines. with the help of these flags, some cpu
c...consuming calculations are only done once per model.
      newmodelanwx=.true.
      newmodelandwdcosnn=.true.
      newmodelep=.true.
      newmodelntdk=.true.
      newmodelsigmav=.true.
      dsprepcalled=.true.


c...Set up Higgses
c      kh(1)=kh1
c      kh(2)=kh2
c      kh(3)=kh3
c      kh(4)=khc
c      mh(1)=mass(kh1)
c      mh(2)=mass(kh2)
c      mh(3)=mass(kh3)
c      mh(4)=mass(khc)


c...Note: the sigma v setup is now moved to dssigmav

c...Note, Higgs decay widths setup is now moved to dshasetup
      call dshasetup

c...For the WIMP annihilation routines in Sun/Earth, setup
c...is moved to dswasetup
      call dswasetup

      end

***********************************************************************
*** Routine dsIBselect selects which channels to include based on the
*** mass degeneracies. Used when ibhow=2.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2007-10-19
***********************************************************************

      subroutine dsibselect

      implicit none
      include 'dsmssm.h'
      include 'dsibcom.h'

      real*8 mn0,zg1,dsabsq
      integer i

      mn0=mass(kn1)
      zg1=dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))

c...Reset all channels
      do i=1,12
            IBflag(i)=0
      enddo 

c...Hard to make 100% safe criterion, uncomment with care
c      if (zg1.gt.0.9d0.and.m0.gt.100) then
c        IBflag(1)=1 ! W+W-
c      endif
      IBflag(1)=1 ! W+W- ! always on for safety

c...Hard to make 100% safe criterion, uncomment with care
c      if (zg1.gt.0.0001d0.and.mn0.gt.800.0d0) then
c        IBflag(2)=1 ! W+H- and W-H+ ! always on for safety
c      endif
      IBflag(2)=1 ! W+H- and W-H+ ! always on for safety

c...Never important, don't include
      IBflag(3)=0 ! H+H-

      if (min(mass(kse1),mass(kse2)).lt.ibmfr*mn0) then
         IBflag(4)=1 ! e+ e-
      endif     

      if (min(mass(ksmu1),mass(ksmu2)).lt.ibmfr*mn0) then
         IBflag(5)=1 ! mu+ mu-
      endif     

      if (min(mass(kstau1),mass(kstau2)).lt.ibmfr*mn0) then
         IBflag(6)=1 ! tau+ tau-
      endif     

      if (min(mass(ksu1),mass(ksu2)).lt.ibmfr*mn0) then
         IBflag(7)=1 ! u u-bar
      endif     

      if (min(mass(ksd1),mass(ksd2)).lt.ibmfr*mn0) then
         IBflag(8)=1 ! d d-bar
      endif     

      if (min(mass(ksc1),mass(ksc2)).lt.ibmfr*mn0) then
         IBflag(9)=1 ! c c-bar
      endif     

      if (min(mass(kss1),mass(kss2)).lt.ibmfr*mn0) then
         IBflag(10)=1 ! s s-bar
      endif     

      if (min(mass(kst1),mass(kst2)).lt.ibmfr*mn0) then
         IBflag(11)=1 ! t t-bar
      endif     

      if (min(mass(ksb1),mass(ksb2)).lt.ibmfr*mn0) then
         IBflag(12)=1 ! b b-bar
      endif     


      return

      end





        


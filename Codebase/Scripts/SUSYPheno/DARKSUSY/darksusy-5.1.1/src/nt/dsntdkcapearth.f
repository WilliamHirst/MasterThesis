

***********************************************************************
*** dsntdkcapearth calculates the capture rate at present from the
*** damour-krauss distribution of wimps. a numerical integration
*** has to be performed instead of the convenient expressions in jkg.
*** author: joakim edsjo (edsjo@physto.se)
*** date: march 18, 1999
***********************************************************************
      real*8 function dsntdkcapearth(mx,sigsi,sigsd)
      implicit none
      include 'dsntdkcom.h'
      include 'dshmcom.h'
      include 'dsntcom.h'
      real*8 mx,sigsi,sigsd
      real*8 dsntcapearthnum,dsntdkgtot10,
     &  dsntcapearthtab
      character*12 vdftmp

c...integrate over the earth according to gould apj 321 (1987) 571.
c...set up things for radial and velocity integration

      gtot10=dsntdkgtot10(mx,sigsi,sigsd) ! d-k eq (5.16)

c...perform integration
      vdftmp=veldfearth
      veldfearth='dk'
      if (nttab.eq.0) then ! use full numerical integration
        dsntdkcapearth=gtot10*dsntcapearthnum(mx,sigsi)
      else  ! use tabulated results
        dsntdkcapearth=gtot10*dsntcapearthtab(mx,sigsi)
      endif
      veldfearth=vdftmp

c      write(*,*) mx,sigsi,dsntdkcapearth

      return
      end

c  dsntcapsunnumi.f:
***********************************************************************
*** dsntcapsunnumi gives the capture rate of neutralinos in the sun
*** given a specified velocity distribution. the integrations over
*** the sun's radius and over the velocity distribution are
*** performed numerically
*** input: mx [ gev ]
***        sigsi [ cm^2 ]
***        sgisd [ cm^2 ]
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
*** Author: Joakim Edsjo
*** Date: 2003-11-26
***********************************************************************

      real*8 function dsntcapsunnumi(mx,sigsi,sigsd,foveru)
      implicit none

      include 'dssun.h'
      real*8 mx,sigsi,sigsd,ma
      real*8 mxx,max,sigsix,rmin,rmax,res,dsntcsint,rx,foveru,
     &  sigsdx
      integer nel,i,ix,vtype,vt
      external dsntcsint,foveru
      common/ntdkint/mxx,max,sigsix,sigsdx,rx,vtype
      common/ntdkint2/ix
c...perform integration as given in Gould, ApJ 521 (1987) 571.
      mxx=mx
      sigsix=sigsi
      sigsdx=sigsd
      dsntcapsunnumi=0.0d0

c...sum over elements
      nel=16 ! include 16 most important elements
      do 10 i=1,nel
        ma=sdma(i) ! nucleus mass
        max=ma
        ix=i
c...begin radial integration
        rmin=0.d0
        rmax=r_sun*100.0d0  ! sun radius in centimeters
        call dshiprecint(dsntcsint,foveru,rmin,rmax,res)
        dsntcapsunnumi=dsntcapsunnumi+res
c        write(*,*) 'Element i=',i,'  C =',res,'  sum=',dsntcapsunnumi 
 10   continue

      return
      end

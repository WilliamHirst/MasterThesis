c  dsntcapearthnumi.f:
***********************************************************************
*** dsntcapearthnumi gives the capture rate of neutralinos in the earth
*** given a specified velocity distribution. the integrations over
*** the earth's radius and over the velocity distribution are
*** performed numerically
*** input: mx [ gev ]
***        sigsi [ cm^2 ]
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
***        vt (velocity type): 1=general type, 2=DK-type (i.e. only
***          include non-zero parts )
*** output: capture rate [ s^-1 ]
*** l.b. and j.e. 1999-04-06
*** Modified by Joakim Edsjo 2003-07-10 to allow for arbitrary external
*** velocity distributions foveru.
***********************************************************************

      real*8 function dsntcapearthnumi(mx,sigsi,foveru,vt)
      implicit none
      include 'dsmpconst.h'

      real*8 mx,sigsi,ma
      real*8 mxx,max,sigsix,rmin,rmax,res,dsntceint,rx,foveru,sigsdx
      integer imass(11),nmass,i,ix,vtype,vt
      external dsntceint,foveru
      common/ntdkint/mxx,max,sigsix,sigsdx,rx,vtype
      common/ntdkint2/ix
c...perform integration as given in gould, apj 521 (1987) 571.
      mxx=mx
      sigsix=sigsi
      dsntcapearthnumi=0.0d0
      vtype=vt

c...sum over elements
      nmass=11 ! include 11 most important elements
      imass(1)=16  ! Oxygen, O
      imass(2)=28  ! Silicon, Si
      imass(3)=24  ! Magnesium, Mg
      imass(4)=56  ! Iron, Fe
      imass(5)=40  ! Calcium, Ca
      imass(6)=30  ! Phosphor, P
      imass(7)=23  ! Sodium, Na
      imass(8)=32  ! Sulfur, S
      imass(9)=59  ! Nickel, Ni
      imass(10)=27 ! Aluminum, Al
      imass(11)=52 ! Chromium, Cr
      do 10 i=1,nmass
        ma=imass(i)*(m_p+m_n)/2.0d0 ! nucleus mass
        max=ma
        ix=imass(i)
c...begin radial integration
        rmin=0.d0
        rmax=6.3782d8  ! earth radius in centimeters
        call dshiprecint(dsntceint,foveru,rmin,rmax,res)
        dsntcapearthnumi=dsntcapearthnumi+res
 10   continue

      return
      end

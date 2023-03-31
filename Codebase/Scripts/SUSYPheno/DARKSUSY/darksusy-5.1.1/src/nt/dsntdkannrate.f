      subroutine dsntdkannrate(m_x,sigsip,sigsdp,sigma_v,arateea,
     &  aratesu)
c_______________________________________________________________________
c
c     wimp annihilation rate in the sun and in the earth
c     in units of 10^24 annihilations per year
c
c     also gives the capture rate and the annih/capt equilibration time
c
c    november, 1995
c    uses routines by p. gondolo and j. edsjo
c    modified by l. bergstrom and j. edsjo and p. gondolo
c    capture rate routines are written by l. bergstrom
c    input:  m_x      - wimp mass
c            sigsip  - spin-indep wimp-proton cross section in cm^2
c            sigsdp  - spin-dep wimp-proton cross section in cm^2
c            sigma_v  - wimp self-annihilation cross section in cm^3/s
c            rescale - rescale factor for local density
c    output: arateea  - 10^24 annihilations per year, earth
c            aratesu  - 10^24 annihilations per year, sun
c            aratedk  - 10^24 annihilations per year, earth including dk
c    slightly modified by j. edsjo.
c    modified by j. edsjo 97-05-15 to match new inv. rate convention
c    modified by j. edsjo 97-12-03 to match muflux3.21 routines.
c    modified by p. gondolo 98-03-04 to detach it from susy routines.
c    added damour-krauss population l.bergstrom 98-09-15
c    added better damour-krauss velocity distribution and numerical
c      integration of the capture rate for these non-gaussian
c      distributions. j. edsjo 99-03-17.
c=======================================================================
      implicit none
      include 'dsprep.h'
      include 'dsntcom.h'
      include 'dswacom.h'

      real*8 arateea,aratesu,m_x,sigsip,sigsdp,sigma_v,
     &  ca_sun,ca_earth,tt_sun,
     &  tt_earth,fluxs,fluxe,
     &  cap_sun,cap_earth,
     &  c_old,c_new,c_new_dot,t_sun,dkalpha,xi,xf,yf
      real*8 dsntcapsun,dsntcapearth2,
     &  dsntdkyf,dsntdkcapearth
      real*8 sarateea,saratesu

      save sarateea,saratesu
      parameter(t_sun=1.5d17)      ! age of the solar system

c --------------------------------------------- start calculations

c=============================================== standard capture rates

c...check if called before for this model, this is used to save
c...cpu time when the rate is calculated several times with different
c...thresholds

      if (.not.newmodelntdk) then
        arateea=sarateea
        aratesu=saratesu
        return
      endif


c...first we take the standard capture
      cap_sun=dsntcapsun(m_x,sigsip,sigsdp)
c      cap_earth=dsntcapearth(m_x,sigsip)  ! old routines
      cap_earth=dsntcapearth2(m_x,sigsip)  ! new full routines
      csu=cap_sun
      cea=cap_earth

c **************************************************************

c...for the sun, the standard calculation applies
      ca_sun=sigma_v/6.6d28*(m_x/20.d0)**(3./2.)
      tausu=1.0d0/dsqrt(cap_sun*ca_sun)
      tt_sun=1.5d17*dsqrt(cap_sun*ca_sun)

      fluxs=cap_sun*0.5d0*tanh(tt_sun)**2
      aratesu = fluxs*1.d-24*3.15d7  ! 10^24 ann. per year

c...for the Earth, we need some of the standard values as well
      ca_earth=wasv/2.3d25*(wamwimp/20.d0)**(3./2.)


c===================================== capture rates with damour-krauss
c...it now gets more complicated due to a non-constant capture rate
c...see astro-ph/99mmnnn for details.

      c_old=cap_earth   ! capture without damour-krauss population
      c_new=dsntdkcapearth(m_x,sigsip,sigsdp)

c...this is the analytical routine
      if (c_new/c_old.gt.1.0d-2) then
        c_new_dot=c_new/t_sun
        ceadk=c_new
        dkalpha=c_new_dot**(2.0/3.0)*ca_earth**(-1.0/3.0)
        xi=c_old/dkalpha
        xf=(c_old+c_new)/dkalpha
        yf=dsntdkyf(xi,xf)
        arateea=0.5d0*dkalpha*yf**2
     &    *1.0d-24*3.15d7     ! 10^24 ann. per year

c...use new numerical routine instead, change? forget about dk in this case?
c      elseif (c_new/c_old.le.1.0d-2.and.c_new.gt.0.0d0) then
c        arateea=dsntratesolve(c_old,c_new,ca_earth,t_sun,1.0d0)
c     &    *1.0d-24*3.15d7       ! 10^24 ann. per year

c...c_new=0, no dk contribution
      else  ! standard calculation for the earth as well
        ca_earth=sigma_v/2.3d25*(m_x/20.d0)**(3./2.)
        tauea=1.0d0/dsqrt(cap_earth*ca_earth)
        tt_earth=1.5d17*dsqrt(cap_earth*ca_earth)

        fluxe=cap_earth*0.5d0*tanh(tt_earth)**2
        arateea = fluxe*1.d-24*3.15d7  ! 10^24 ann. per year
      endif

      sarateea=arateea
      saratesu=aratesu
      newmodelntdk=.false.
 
      return

      end

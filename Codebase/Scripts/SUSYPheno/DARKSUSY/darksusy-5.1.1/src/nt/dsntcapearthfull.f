***********************************************************************
*** full capture rate routines for the earth.
*** this set of routines use the full expressions for the capture rate
*** in the earth as given in gould, apj 321 (1987) 571.
*** dsntcapearthfull thus replaces dsntcapearth which use the approximations
*** given in the jkg review. these routines assume a maxwell-
*** boltmann velocity distribution. for the damour-krauss population
*** of wimps, the routine dsntcapearthnumi should be used instead.
***********************************************************************

      real*8 function dsntcapearthfull(mx,sigsi,v_star,v_bar,rho_x)
c----------------------------------------------------------------------
c         capture rate in the earth
c *** full: use formulas by gould ap.j. 321 (1987) 571
c       mass fractions and phi_i from jkg review
c       mx: neutralino mass
c       sigsi: spin independent cross section in units of cm^2
c       v_star: solar system velocity through halo (220 km/s in standard case)
c       v_bar: 3D velocity dispersion of wimps (270 km/s in standard case)
c       rho_x: local wimp density (units of gev/cm**3)
c       lars bergstrom 1998-09-21
c       Modified by J. Edsjo, 2003-11-22
c References: gould: Gould ApJ 321 (1987) 571
c----------------------------------------------------------------------
       implicit none
       include 'dsmpconst.h'
       real*8 sigsi,v_star,v_bar,rho_x
       real*8 mx,aa(11),sigsinorm
       real*8 m(11),f(11),phi(11),mn,sla1,fact,sig_a,dsntsefull
       integer i
       mn=(m_p+m_n)/2.0d0
       m(1)=16.d0*mn
       m(2)=28.d0*mn
       m(3)=24.d0*mn
       m(4)=56.d0*mn
       m(5)=40.d0*mn
       m(6)=30.d0*mn
       m(7)=23.d0*mn
       m(8)=32.d0*mn
       m(9)=59.d0*mn
       m(10)=27.d0*mn
       m(11)=52.d0*mn

c...These are the values in jkg
       f(1)=0.3d0
       f(2)=0.15d0
       f(3)=0.14d0
       f(4)=0.3d0
       f(5)=0.015d0
       f(6)=0.011d0
       f(7)=0.004d0
       f(8)=0.05d0
       f(9)=0.03d0
       f(10)=0.0d0
       f(11)=0.0d0

       phi(1)=1.2d0
       phi(2)=1.2d0
       phi(3)=1.2d0
       phi(4)=1.6d0
       phi(5)=1.2d0
       phi(6)=1.2d0
       phi(7)=1.2d0
       phi(8)=1.6d0
       phi(9)=1.6d0
       phi(10)=1.0d0
       phi(11)=1.0d0

c...These are the new values based on McDonough, Treatise on Geochemistry,
c...Vol 2, Elsevier, 2003.
       f(1)=0.298d0
       f(2)=0.162d0
       f(3)=0.154d0
       f(4)=0.319d0
       f(5)=0.0171d0
       f(6)=0.00071d0
       f(7)=0.00183d0
       f(8)=0.0063d0
       f(9)=0.0181d0
       f(10)=0.0159d0
       f(11)=0.0047d0

c...old 1.6 core, 1.2 mantle averages
       phi(1)=1.2d0
       phi(2)=1.24d0
       phi(3)=1.2d0
       phi(4)=1.546d0
       phi(5)=1.2d0
       phi(6)=1.56d0
       phi(7)=1.2d0
       phi(8)=1.59d0
       phi(9)=1.57d0
       phi(10)=1.20d0
       phi(11)=1.44d0

c...new full average
       phi(1)=1.20d0
       phi(2)=1.25d0
       phi(3)=1.20d0
       phi(4)=1.54d0
       phi(5)=1.2d0
       phi(6)=1.56d0
       phi(7)=1.20d0
       phi(8)=1.58d0
       phi(9)=1.56d0
       phi(10)=1.20d0
       phi(11)=1.45d0

       aa(1)=16.d0
       aa(2)=28.d0
       aa(3)=24.d0
       aa(4)=56.d0
       aa(5)=40.d0
       aa(6)=30.d0
       aa(7)=23.d0
       aa(8)=32.d0
       aa(9)=59.d0
       aa(10)=27.d0
       aa(11)=52.d0

       sla1=0.d0
       sigsinorm=sigsi*1.d40 ! units of 10**-40 cm**2
       fact=1.d5*3.35d11 ! velocity in cm/sec; earth mass in 10^40 gev
       do 20 i=1,11
       sig_a=1.*aa(i)**2*(mx*m(i))**2/(mx+m(i))**2/(mx*mn)**2*(mx+mn)**2
              ! to go from sigsi to sigma_a
              !  see d-k (7.1), or Eq. (9-25) in jkg w/ A=m(i)/mn
              ! sig_a * sigsi is cross section on nucleus i
c...We will now sum up the contributions from the different elements in
c...the Earth. We use Eq. (A10) in gould. dsntsefull returns everything in
c...(A10) except for sigma n n_w v_bar. Multiply by these and perform
c...the integral over the Earth's volume by multiplyting by the
c...number of atoms of the given type.
         sla1=sla1+
     &   f(i)/m(i)*sigsinorm*sig_a*
     &   dsntsefull(mx,m(i),v_star,v_bar,phi(i))*(rho_x/mx)*v_bar
     &          *fact
c         write(*,*) 'caprate: ',sla1
c      capture rate on earth per second
 20    continue
        dsntcapearthfull=sla1
        return
        end

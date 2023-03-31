      subroutine dsbsgpre2007full(ratio,flag)

***********************************************************************
* Routine that calculates the b-->s+gamma branching rate              *
* The Standard Model contribution is taken from                       *
* Gambino and Misiak,Nucl. Phys. B611 (2001) 338                      *
* with new 'magic numbers' from Buras et al., hep-ph/0203135          *
* Input: flag:  0 = only standard model                               *
*               1 = standard model plus SUSY corrections              *
* Output: ratio = BR(b -> s gamma)                                    *
* author:Mia Schelke, schelke@physto.se, 2003-03-27                   *
***********************************************************************
      implicit none
      include 'dsmssm.h'
      integer flag
      real*8 ratio
      real*8 mu0,mt,ckmratio
      real*8 alpha3mu
      real*8 pofe
      real*8 dsbsgbofe,dsbsgckm,dsbsgalpha3,dsabsq
      complex*16 dsbsgkt, dsbsgkc
      real*8 bsgbr
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     We here set mt=m_t(m_t)

      mt=mtmt   ! Use DarkSUSY value at scale mt


c     We here set mu0. In this routine, the scale mu0 enters in
c     P(E0) together with the term Kt, hence we should use the scale
c     mu0=mt

      mu0=mt
 
c     value of alpha3mu=alpha3(\mu_0)
c     has here been set to the value at the scale mu0
c     Note: With mu0=mt, the terms in pofe below vanishes, so this is only
c     included for consistency if other scales would be used in the future
      alpha3mu=dsbsgalpha3(mu0)


c     We first calculate P(E_0) in eq.(3.5),
c     Nucl. Phys. B611 (2001) 338 
c     For r(\mu_0) we use the value in eq. (4.1) with mu0=mt
c     For e_ew we use the value of footnote 8 p.12


c     We also include both the real and imaginary parts of kt and kc.
c     Ignoring the imaginary parts (as they do in hep-ph/0104034) gives
c     a result 0.5% lower.
      
      pofe=dsbsgbofe(flag)
     &     +dsabsq(dsbsgkc()+(1.d0+(alpha3mu/pi)*log(mu0**2/mt**2))
     &        *0.578d0*dsbsgkt(flag)
     &       +dcmplx(0.0071d0,0.d0))

c    We here calculate the branching ratio for b-->s+gamma
c    by the formula in eq(3.1) of Nucl. Phys. B611 (2001) 338
c    The semileptonic br. ratio is taken from app. A
c    The ratio of the ckm-elements is taken from eq.(4.5)

c    When the flag is set, we instead call the function
c    dsbsgckm.f to get the susy coreected value of the 
c    ratio of the relevant ckm elements.
c    The value of \alpha_em=\alpha_em^(on shell) is 
c    taken from app A. 
c    The value of C is taken from eq.(4.2)
c    The value of N(E_0) is taken from eq.(4.7), but this is
c    for another value of E_0 (1.6 GeV) instead of E_0 = mb/20 (0.23 GeV)
c    that we use. However, the energy dependence (if any...) of N(E_0) seems
c    to be very weak, so this approximation should be OK.

      ckmratio=0.971d0

      if (flag.eq.1) then
       ckmratio=dsbsgckm()
c       write(*,*) 'ckmratio=',ckmratio
      endif     

      bsgbr=0.1045d0*ckmratio*6.d0*(1.d0/137.036d0)/(pi*0.575d0)
     &      *(pofe+0.0036d0)

      ratio=bsgbr

c      write (*,*) 'BR(bsg)',bsgbr

      return
      end


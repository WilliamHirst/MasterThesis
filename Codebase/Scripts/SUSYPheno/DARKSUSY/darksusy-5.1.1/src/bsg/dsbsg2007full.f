      subroutine dsbsg2007full(ratio,flag)

***********************************************************************
* Routine that calculates the b-->s+gamma branching rate              *
* This routine was written after a big shift in the theoretical       *
* calculation of the SM Branching ratio in 2006.                      *
* From 3.7e-4 to 3.15e-4.                                             *
* This routine implements an approximate solution which works when    *
* there are only small Beyond-SM (BSM) corrections                    *
* For the implemented solution we refer to:                           *
*   M. Misiak et al: hep-ph/0609232                                   *
* + M. Misiak and M. Steinhauser: hep-ph/0609241                      *
* + M. Misiak and M. Steinhauser, private communication               *
* + P. Gambino, private communication.                                *
* (+ref for the SUSY contributions, see the individual files)         *
* The implemented solution is:                                        *
* BR * 10^4  = (3.15 +_ 0.23) - 8.0 * deltaBSM[C7]                    *
*               - 1.9 * deltaBSM[C8] + O(deltaBSM^2)                  *
* with a matching scale m_0 = 2 m_W = 160 GeV                         *
* for other parameters see hep-ph/0609241 appendix A.                 *
* but we use the old matching scale mu_0 = m_t(m_t) = 166.6 GeV       *
* Input: flag:  0 = only standard model                               *
*               1 = standard model plus SUSY corrections              *
* Output: ratio = BR(b -> s gamma)                                    *
* author:Mia Schelke, schelke@physto.se, 2007                         *
***********************************************************************
      implicit none
      include 'dsmssm.h'
      integer flag
      real*8 ratio
      real*8 mu0,mt
      real*8 alpha3mu 
      real*8 dsbsgalpha3
      real*8 deltac7,deltac70,deltac71
      real*8 deltac8,deltac80,deltac81
      real*8 dsbsgc70h2,dsbsgc70susy
      real*8 dsbsgc71h2,dsbsgc71chisusy,dsbsgc71wsusy
      real*8 dsbsgc71phi1susy,dsbsgc71phi2susy
      real*8 dsbsgc80h2,dsbsgc80susy
      real*8 dsbsgc81h2,dsbsgc81chisusy,dsbsgc81wsusy
      real*8 dsbsgc81phi1susy,dsbsgc81phi2susy
      real*8 bsgbr
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c     We here set mt=m_t(m_t)

      mt=mtmt   ! Use DarkSUSY value at scale mt


c     We here set mu0. 

      mu0=mt
 
c     value of alpha3mu=alpha3(\mu_0)
c     has here been set to the value at the scale mu0

      alpha3mu=dsbsgalpha3(mu0)


      bsgbr=3.15d0


      if (flag.eq.1) then
        deltac70=dsbsgc70h2()+dsbsgc70susy()
        deltac71=dsbsgc71h2()+dsbsgc71chisusy()
     &              +dsbsgc71wsusy()+dsbsgc71phi1susy()
     &              +dsbsgc71phi2susy()
        deltac80=dsbsgc80h2()+dsbsgc80susy()
        deltac81=dsbsgc81h2()+dsbsgc81chisusy()
     &              +dsbsgc81wsusy()+dsbsgc81phi1susy()
     &              +dsbsgc81phi2susy()
        deltac7=deltac70+(alpha3mu/(4.0d0*pi))*deltac71
        deltac8=deltac80+(alpha3mu/(4.0d0*pi))*deltac81
        bsgbr=bsgbr-8.0d0*deltac7-1.9d0*deltac8
      endif

c To test just the two-higgs-doublet II model:

c      if (flag.eq.1) then
c        deltac70=dsbsgc70h2()
c        deltac71=dsbsgc71h2()
c        deltac80=dsbsgc80h2()
c        deltac81=dsbsgc81h2()
c        deltac7=deltac70+(alpha3mu/(4.0d0*pi))*deltac71
c        deltac8=deltac80+(alpha3mu/(4.0d0*pi))*deltac81
c        bsgbr=bsgbr-8.0d0*deltac7-1.9d0*deltac8
c      endif

c Might want to insert here a test for the smallness 
c of the susy corrections.

      ratio=bsgbr*1.0d-4

c      write (*,*) 'BR(bsg)',bsgbr

      return
      end


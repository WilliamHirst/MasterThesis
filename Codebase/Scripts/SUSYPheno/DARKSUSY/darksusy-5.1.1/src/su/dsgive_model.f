      subroutine dsgive_model(amu,am2,ama,atanbe,amsq,atm,abm)
c----------------------------------------------------------------------
c
c     To specify the supersymmetric parameters of a model.
c     Inputs:
c        amu - mu parameter (GeV)
c        am2 - M2 parameter (GeV)
c        ama - Mass of the CP-odd Higgs boson A (or H3)
c        atanbe - ratio of Higgs vacuum expecation values, tan(beta)
c        amsq - common sfermion mass scale, M_sq_tilde (GeV)
c        atm - trilinear term in units of amsq, top sector 
c        abm - trilinear term in units of amsq, bottom sector
c     Outputs:
c        The common blocks are set corresponding to the values above
c     Author: Paolo Gondolo, gondolo@mppmu.mpg.de
c     Date: 2000
c     Modified: Joakim Edsjo, edsjo@physto.se
c        2001-02-13 - setting of idtag taken away
c----------------------------------------------------------------------

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 amu,am2,ama,atanbe,amsq,atm,abm,msq
      real*8 dsgf2s2thw
      mu = amu
      m2 = am2
      ma = ama
      tanbe = atanbe
      msq = amsq
c... further assumptions
c      mqtilde=amsq   ! Still needed by FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      ! GUT relation on gaugino mass parameters
c...The following is set in dsinit and dssuconst_couplings
c      s2wmz=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),3)
      m3 =  m2 * (alph3mz*s2wmz)/alphem    ! alph3->alph3mz 020912 (JE)
      m1 = m2 * 5.0d0/3.0d0 * s2wmz/(1.d0-s2wmz)
      ! Slepton and squark mass matrices at the weak scale
      mass2q(3) = msq**2
      mass2u(3) = msq**2
      mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
      mass2l(3) = mass2d(3)
      mass2e(3) = mass2d(3)
      asofte(3) = 0.d0
      asoftu(3) = atm*msq
      asoftd(3) = abm*msq
      do i=1,2
         mass2q(i) = mass2d(3)
         mass2u(i) = mass2d(3)
         mass2d(i) = mass2d(3)
         mass2l(i) = mass2d(3)
         mass2e(i) = mass2d(3)
         asofte(i) = 0.d0
         asoftu(i) = 0.d0
         asoftd(i) = 0.d0
      enddo

      modeltype=1 ! MSSM-7 model

      return
      end

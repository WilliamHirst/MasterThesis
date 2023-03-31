      subroutine dsgive_model13(amu,am1,am2,ama,atanbe,
     &  amse,amsmu,amstau,amsq,
     &  atm,abm,ataum,aotherm)
c----------------------------------------------------------------------
c
c     To specify the supersymmetric parameters of a model.
c     Inputs:
c        amu - mu parameter (GeV)
c        am1 - M1 parameter (GeV)
c        am2 - M2 parameter (GeV)
c        ama - Mass of the CP-odd Higgs boson A (or H3)
c        atanbe - ratio of Higgs vacuum expecation values, tan(beta)
c        amse - selectron mass scale (GeV)
c        amsmu - smuon mass scale (GeV)
c        amstau - stau mass scale (GeV)
c        amsq - common squark mass scale, M_sq_tilde (GeV)
c        atm - trilinear term in units of amsq, top sector 
c        abm - trilinear term in units of amsq, bottom sector
c        ataum -  trilinear term in units of amstau, stau sector
c        aotherm - trilinear term for remaining squarks and sleptons
c                  in units of amsq, amse and amsmu respectively
c     Outputs:
c        The common blocks are set corresponding to the values above
c     Author: Joakim Edsjo, edsjo@physto.se
c     Date:   2007-12-17
c----------------------------------------------------------------------

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 amu,am1,am2,ama,atanbe,
     &  amse,amsmu,amstau,amsq,
     &  atm,abm,ataum,aotherm,msq
      mu = amu
      m2 = am2
      ! GUT relation on gaugino mass parameters
      m3 =  m2 * (alph3mz*s2wmz)/alphem    ! alph3->alph3mz 020912 (JE)
c      m1 = m2 * 5.0d0/3.0d0 * s2wmz/(1.d0-s2wmz)
      m1=am1
      ma = ama
      tanbe = atanbe
      msq = amsq
c... further assumptions
c      mqtilde=amsq   ! Still needed by FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      ! Slepton and squark mass matrices at the weak scale
      mass2q(3) = msq**2
      mass2u(3) = msq**2
      mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
      asoftu(3) = atm*msq
      asoftd(3) = abm*msq
      do i=1,2
         mass2q(i) = mass2d(3)
         mass2u(i) = mass2d(3)
         mass2d(i) = mass2d(3)
         asoftu(i) = aotherm*msq
         asoftd(i) = aotherm*msq
      enddo

      mass2l(1) = amse**2
      mass2l(2) = amsmu**2
      mass2l(3) = amstau**2
      mass2e(1)=mass2l(1)
      mass2e(2)=mass2l(2)
      mass2e(3)=mass2l(3)

      asofte(1) = aotherm*amse
      asofte(2) = aotherm*amsmu
      asofte(3) = ataum*amstau

      modeltype=5 ! MSSM-13 model

      return
      end

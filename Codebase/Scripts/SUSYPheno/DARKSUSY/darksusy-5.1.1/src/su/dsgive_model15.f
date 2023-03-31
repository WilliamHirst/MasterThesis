      subroutine dsgive_model15(amu,am1,am2,ama,atanbe,
     &  amse,amsmu,amstau,amsq1,amsq2,amsq3,
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
c        amsq1 - squark mass scale, 1st generation (GeV)
c        amsq2 - squark mass scale, 2nd generation (GeV)
c        amsq3 - squark mass scale, 3rd generation (GeV)
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
     &  amse,amsmu,amstau,amsq1,amsq2,amsq3,
     &  atm,abm,ataum,aotherm
      mu = amu
      m2 = am2
      ! GUT relation on gaugino mass parameters
      m3 =  m2 * (alph3mz*s2wmz)/alphem    ! alph3->alph3mz 020912 (JE)
c      m1 = m2 * 5.0d0/3.0d0 * s2wmz/(1.d0-s2wmz)
      m1=am1
      ma = ama
      tanbe = atanbe
c... further assumptions
c      mqtilde=amsq   ! Still needed by FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      ! Slepton and squark mass matrices at the weak scale

      mass2q(1) = amsq1**2
      mass2u(1) = amsq1**2
      mass2d(1) = 2.0d0*mass2q(1)-mass2u(1)
      asoftu(1) = aotherm*amsq1
      asoftd(1) = aotherm*amsq1

      mass2q(2) = amsq2**2
      mass2u(2) = amsq2**2
      mass2d(2) = 2.0d0*mass2q(2)-mass2u(2)
      asoftu(2) = aotherm*amsq2
      asoftd(2) = aotherm*amsq2

      mass2q(3) = amsq3**2
      mass2u(3) = amsq3**2
      mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
      asoftu(3) = atm*amsq3
      asoftd(3) = abm*amsq3

      mass2l(1) = amse**2
      mass2l(2) = amsmu**2
      mass2l(3) = amstau**2
      mass2e(1)=mass2l(1)
      mass2e(2)=mass2l(2)
      mass2e(3)=mass2l(3)

      asofte(1) = aotherm*amse
      asofte(2) = aotherm*amsmu
      asofte(3) = ataum*amstau

      modeltype=8 ! MSSM-15 model

      return
      end

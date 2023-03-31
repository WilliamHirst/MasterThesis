      subroutine dsgive_model25(am1,am2,am3,amu,ama,atanbe,
     &  amsqL1,amsql2,amsql3,amsqRu,amsqRc,amsqRt,amsqRd,
     &  amsqRs,amsqRb,amslL1,amslL2,amslL3,amslRe,
     &  amslRmu,amslRtau,atm,abm,ataum,aemum)
c----------------------------------------------------------------------
c
c     To specify the supersymmetric parameters of a model.
c     Inputs:
c        Gaugino Sector
c        am1 - M1 parameter (GeV)
c        am2 - M2 parameter (GeV)
c        am3 - M3 parameter (GeV)
c
c        Higgs Sector
c        amu - mu parameter (GeV)
c        ama - Mass of the CP-odd Higgs boson A (or H3)
c        atanbe - ratio of Higgs vacuum expecation values, tan(beta)
c
c        Sfermion sector
c        amsqL1 - left type squark mass scale, 1st generation
c        amsqL2 - left type squark mass scale, 2nd generation
c        amsqL3 - left type squark mass scale, 3rd generation
c        amsqRu - right type sup squark mass scale
c        amsqRc - right type scharm squark mass scale
c        amsqRt - right type stop squark mass scale
c        amsqRd - right type sdown squark mass scale
c        amsqRs - right type sstrange squark mass scale
c        amsqRb - right type sbottom squark mass scale
c        amslL1 - left type slepton mass scale, 1st generation
c        amslL2 - left type slepton mass scale, 2nd generation
c        amslL3 - left type slepton mass scale, 3rd generation
c        amslRe - right type selectron mass scale
c        amslRmu - right type smuon mass scale
c        amslRtau - right type stau mass scale
c
c        Trilinear Couplings
c        atm - trilinear term, top sector, (in units of ?) 
c        abm - trilinear term, bottom sector, (in units of ?) 
c        ataum - trilinear term, stau sector
c        aemum - trilinear term, electron and muon sector
c
c     Outputs:
c        The common blocks are set corresponding to the values above
c     
c     Author: Joakim Edsjo, edsjo@physto.se
c     Date:   2007-12-17
c     Modified by Hamish Silverwood, moving from 15 to 25 parameters
c     Date:   2011-06-30
c----------------------------------------------------------------------

      implicit none
      include 'dsmssm.h'
      real*8 am1,am2,am3,amu,ama,atanbe,
     &  amsqL1,amsql2,amsql3,amsqRu,amsqRc,
     &  amsqRt,amsqRd,amsqRs,amsqRb,
     &  amslL1,amslL2,amslL3,amslRe,amslRmu,
     &  amslRtau,atm,abm,ataum,aemum
      

c     Gaugino Sector
      m1 = am1
      m2 = am2
      m3 = am3

c     Higgs Sector
      mu = amu
      ma = ama
      tanbe = atanbe

c     Sfermion sector
      mass2q(1) = amsqL1**2
      mass2q(2) = amsqL2**2
      mass2q(3) = amsqL3**2

      mass2u(1) = amsqRu**2
      mass2u(2) = amsqRc**2
      mass2u(3) = amsqRt**2

      mass2d(1) = amsqRd**2
      mass2d(2) = amsqRs**2
      mass2d(3) = amsqRb
      mass2l(1) = amslL1**2
      mass2l(2) = amslL2**2
      mass2l(3) = amslL3**2

      mass2e(1) = amslRe**2
      mass2e(2) = amslRmu**2
      mass2e(3) = amslRtau**2

c     Trilinear Couplings
      asoftu(1) = 0
      asoftu(2) = 0
      asoftu(3) = atm

      asoftd(1) = 0
      asoftd(2) = 0
      asoftd(3) = abm
      
      asofte(1) = aemum
      asofte(2) = aemum
      asofte(3) = ataum

      modeltype=25 ! MSSM-25 model

c      write(*,*) '******************************************'
c      write(*,*) 'In dsgive_model25 we have'
c      write(*,*) m1,m2,m3 
c      write(*,*) mu,ma,tanbe
c      write(*,*) mass2q(1),mass2q(2),mass2q(3)
c      write(*,*) mass2u(1),mass2u(2),mass2u(3)
c      write(*,*) mass2d(1),mass2d(2),mass2d(3)
c      write(*,*) mass2l(1),mass2l(2),mass2l(3)
c      write(*,*) mass2e(1),mass2e(2),mass2e(3)
c      write(*,*) asoftu(1),asoftu(2),asoftu(3)
c      write(*,*) asoftd(1),asoftd(2),asoftd(3)
c      write(*,*) asofte(1),asofte(2),asofte(3)
c      write(*,*) '******************************************'

      return
      end

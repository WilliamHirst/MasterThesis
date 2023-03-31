      subroutine dsgive_model_isasugra(m0,mhf,a0,sgnmu,tgbeta)
c----------------------------------------------------------------------
c
c     To specify the supersymmetric parameters of a model.
c     Inputs:
c        m0 - m0 parameter (GeV)
c        mhf - m_{1/2} parameter (GeV)
c        a0 - trilinear term (GeV)
c        sgnmu - sign of mu (+1.0d0 or -1.0d0)
c        tgbeta - ratio of Higgs vacuum expecation values, tan(beta)
c     Outputs:
c        The common blocks are set corresponding to the values above
c     Author: Joakim Edsjo, edsjo@physto.se
c        2002-03-12
c----------------------------------------------------------------------

      implicit none
      include 'dsmssm.h'
      integer i
      real*8 m0,mhf,a0,sgnmu,tgbeta

      m0var=m0
      mhfvar = mhf
      a0var = a0
      sgnmuvar = sgnmu
      tgbetavar = tgbeta         

      modeltype=3 ! mSUGRA - Higgs widths routines need to know

      return
      end

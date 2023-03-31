**********************************************************************
*** function dspbtpb gives the antiproton kinetic energy at the helio-
*** sphere as a function of the kinetic energy at the earth.
***   input:
***     tp - antiproton kinetic energy in gev
*** date: 98-02-10
**********************************************************************

      real*8 function dspbtpb(tp)
      implicit none

      include 'dshmcom.h'
      include 'dsidtag.h'
      include 'dsmpconst.h'

c----------------------------------------------------------------------
      real*8 pp,tp,e,kkc,epb,ec,ppb,smod

c      parameter (kkc=1.015d0)  ! critical momentum, gev
c      parameter (smod=.495d0)  ! solar modulation, gev
c      parameter (smod=.365d0)  ! solar modulation, gev

c...parameters for pbar paper with lars, piero and joakim:
      parameter (kkc=0.0d0)    ! critical momentum, gev
      parameter (smod=.500d0)  ! solar modulation, gev

c----------------------------------------------------------------------

c...solar modulation: calculate kinetic energy at heliosphere boundary
      e=m_p+tp
      pp=sqrt(2*m_p*tp+tp**2)
      if (pp.ge.kkc) then
         epb=e+smod
c         tjac=1.0d0
      else
         ec=sqrt(m_p**2+kkc**2)
         epb=kkc*log((pp+e)/(kkc+ec))+ec+smod
c         tjac=kkc/(k+e)
      endif
      ppb=sqrt(max(epb**2-m_p**2,0.0d0))
      dspbtpb=epb-m_p   ! pbar kinetic energy at boundary

      return
      end






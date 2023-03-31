      function dsandwdcoss(costheta)
c_______________________________________________________________________
c  10^15*annihilation differential invariant rate.
c  input:
c    p - initial cm momentum (real) for lsp annihilations via common
c    costheta - cosine of c.m. annihilation angle
c  uses dwdcos
c  used for gaussian integration with gadap.f
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-01-09
c=======================================================================
      implicit none
      real*8 dsandwdcoss,costheta,tmp
      real*8 dsandwdcos
      external dsandwdcos

c-----------------------------------------------------------------------
      real*8 pd
      common /gadint/ pd
c-----------------------------------------------------------------------

      tmp=dsandwdcos(pd,dble(costheta))
      tmp=tmp*1.0d15
      dsandwdcoss=tmp

      end









































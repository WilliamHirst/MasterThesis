      function dsandwdcosy(y)
c_______________________________________________________________________
c  10^15*annihilation differential invariant rate.
c  the integration variable is changed from cos(theta) to
c    y=1/(mx^2+2p^2(1-cos(theta))) for cos(theta)>0 and to
c    y=1/(mx^2+2p^2(1-cos(theta))) for cos(theta)<0.
c  this avoids the poles at cos(theta)=+-1
c  integrate this function from 1/(mx^2+2p^2) to 1/mx^2 to get
c  1d15 times the integral from cos(theta)=0 to 1.
c  input:
c    y - initial cm momentum (real) for lsp annihilations via common
c  uses dwdcos
c  used for gaussian integration with gadap.f
c  author: joakim edsjo (edsjo@physto.se)
c  date: 98-05-03
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 dsandwdcosy,y
      real*8 dsandwdcos,tmp,costheta,jac,mx
      external dsandwdcos

c-----------------------------------------------------------------------
      real*8 pd
      common /gadint/ pd

      real*8 alph,bet
      common /yint/ alph,bet
      save /yint/
c-----------------------------------------------------------------------

      mx=mco(1)
c...cos(theta)=-1..0
c      costheta=-((1.0d0/y-mx**2)/(2.0d0*pd**2)-1.0d0)
c      costheta=-max(((1.0d0/y**2-mx**2)/(2.0d0*pd**2)-1.0d0),-1.0d0)
c      costheta=-max(((y**(1.0d0/alph)-mx**2)/(bet*pd**2)-1.0d0),-1.0d0)
c...cos(theta)=0..1
      costheta=min(1.0d0-(y**(1.0d0/alph)-mx**2)/(bet*pd**2),1.0d0)
      tmp=dsandwdcos(pd,costheta)+dsandwdcos(pd,-costheta)
      tmp=tmp*1.0d15
c...add jaobian from change of variables
      jac=-1.0d0/(alph*bet*pd**2*y**((alph-1.0d0)/alph))
      tmp=tmp*jac
      dsandwdcosy=tmp

      end









































************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98
************************************************************************


************************************************************************
      real*8 function dsepvvcut(e)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 e,u
      ametric=1.0d0/3.0d0**alphexp
      av=ametric/(1.0d0-alphexp)
      u=1.0d0/e
      dsepvvcut=u+av*u**(1.0d0-alphexp)
      return
      end

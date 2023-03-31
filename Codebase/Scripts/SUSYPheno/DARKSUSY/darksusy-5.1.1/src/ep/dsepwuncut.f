************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98
************************************************************************


************************************************************************
      real*8 function dsepwuncut(v)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 v
      dsepwuncut=((1.0d0-alphexp)*v)**
     &  (-((2.0d0-alphexp)/(1.0d0-alphexp)))
      return
      end

************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98
***           jul-06-99 paolo gondolo - calls to dshunt, ee, vv
************************************************************************


************************************************************************
      real*8 function dsepeecut(v,tabindx)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 v
      integer tabindx
      dsepeecut=etable(tabindx)
      return
      end

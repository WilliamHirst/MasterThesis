************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98
************************************************************************


c     this is a power law e^-2
c     this function is [d<sigma velocity>/de](v)
      real*8 function dsepdsigv_de(v,vmin)
      implicit none
      real*8 v,vmin
      if (v.lt.vmin) then
         dsepdsigv_de=0.0d0
         return
      endif
      dsepdsigv_de=1.0e-27*(0.4d0*v)**5
      return
      end

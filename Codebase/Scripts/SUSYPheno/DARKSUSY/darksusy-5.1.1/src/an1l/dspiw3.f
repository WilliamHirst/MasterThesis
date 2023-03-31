c====================================================================
c
c   auxiliary function used in:  
c   dsanglglre.f  dsrepfbox.f  dsrepgh.f  dsrepw.f 
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dspiw3(a,b)
      implicit none
      real*8 a,b,ax,bx
      real*8 sla1,sla2,abserr,alist,blist,
     &     elist,epsabs,epsrel,result,rlist
      integer ier,iord,last,limit,neval
      dimension alist(1000),blist(1000),elist(1000),iord(1000),
     &     rlist(1000)
      external dspiw3i
      common/dspiw3c/ax,bx
      epsabs=1.d-17
      epsrel=1.d-17
      limit=1000
      ax=a
      bx=b
      sla1=0.d0
      sla2=1.d0
      dspiw3=0.d0
      call dqagse(dspiw3i,sla1,sla2,epsabs,epsrel,limit,result,
     &     abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dspiw3=dspiw3+result
      return
      end


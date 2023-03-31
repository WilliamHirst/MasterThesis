**********************************************************************
*** function called in dspbtd15char
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbcharpar1(x)
      implicit none
      include 'dspbprivate.h'
      real*8 x,xx
      real*8 abserr,alist,blist,
     &  elist,epsabs,epsrel,result,rlist,dbesj0
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dspbcharpar2
      common /pbchar/xx
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      limit=5000
      xx=x
      call dqagseb(dspbcharpar2,0.d0,1.d0,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      dspbcharpar1=x*result*dbesj0(zero*x)
      return
      end

c====================================================================
c
c   auxiliary function used in:  
c   dsi_14.f  dsi_24.f  dsi_34.f  dsilp2.f dsti_5.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dslp2(c1,c2)
c this function compute the integral between 0 and 1 of 1/x*log(1+c1*x+c2*x**2)
      implicit none
      include 'dsge.h'
      real*8 c1,c2,cc1,cc2,delta,q,x1i,x2i,dsdilog
      real*8 sla1,sla2,abserr,alist,blist,
     &  elist,epsabs,epsrel,result,rlist
      integer ier,iord,last,limit,neval
      dimension alist(1000),blist(1000),elist(1000),iord(1000),
     & rlist(1000)
      external dsilp2
      common/dslp2c/cc1,cc2
      epsabs=1.d-17
      epsrel=1.d-17
      limit=1000
      sla1=0.d0
      sla2=1.d0
      dslp2=0.d0
      cc1=c1
      cc2=c2
      if(cc1.eq.0.d0.and.cc2.eq.0.d0) then
      dslp2=0.d0
      goto 100
      endif
      delta=cc1*cc1-4.d0*cc2
      if(delta.lt.0) then
c compute the integral with the numerical integration
        call dqagse(dsilp2,sla1,sla2,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c        write(*,*) ier,result,abserr
        dslp2=dslp2+result
      else
c compute the integral in terms of dsdilogarithms
c x1i and x2i are the inverse of the two roots
        if(cc1.eq.0.d0) then
          x1i=2.d0*cc2/dsqrt(delta)
          x2i=-2.d0*cc2/dsqrt(delta)
        else
          q=-0.5d0*(cc1+cc1/dabs(cc1)*dsqrt(delta))
          x1i=cc2/q
          x2i=q
        endif
c first add the term from x1i
        if(x1i.eq.1d0) then
          dslp2=dslp2-pi*pi/6.d0
        elseif(x1i.eq.-1d0) then
          dslp2=dslp2-(-pi*pi/12.d0)
        elseif(x1i.gt.1.d0) then
          dslp2=dslp2-(-dsdilog(1.d0/x1i)+pi*pi/3.d0
     &      -0.5d0*(dlog(x1i))**2)
        elseif(x1i.lt.-1.d0) then
          dslp2=dslp2-(-dsdilog(1.d0/x1i)-pi*pi/6.d0
     &      -0.5d0*(dlog(-x1i))**2)
        else
          dslp2=dslp2-dsdilog(x1i)
        endif
c and then add the term from x2i
        if(x2i.eq.1d0) then
          dslp2=dslp2-pi*pi/6.d0
        elseif(x2i.eq.-1d0) then
          dslp2=dslp2-(-pi*pi/12.d0)
        elseif(x2i.gt.1.d0) then
          dslp2=dslp2-(-dsdilog(1.d0/x2i)+pi*pi/3.d0
     &      -0.5d0*(dlog(x2i))**2)
        elseif(x2i.lt.-1.d0) then
          dslp2=dslp2-(-dsdilog(1.d0/x2i)-pi*pi/6.d0
     &      -0.5d0*(dlog(-x2i))**2)
        else
          dslp2=dslp2-dsdilog(x2i)
        endif
      endif
 100  continue
      return
      end

      subroutine dsrdeqn(wrate,x0,x1,y1,xf,nfcn)
c_______________________________________________________________________
c  solve the relic density evolution equation by means of an implicit
c    trapezoidal method with adaptive stepsize and termination.
c  input:
c    wrate - invariant annihilation rate (real, external)
c    x0 - initial mass/temperature (real)
c    x1 - final mass/temperature (real)
c    y1 - final number/entropy densities (real)
c    nfcn - number of calls to wrate (integer)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdrhs.
c  called by dsrdens.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996
c  modified: joakim edsjo (edsjo@physto.se) 961212
c            Paolo Gondolo, factor added 2003
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wrate,x0,y1,x1,xf,h1,safety,errcon,yo,eps1,xo,eps,factor
      external wrate
      integer maxstp,nstp,kstp,nfcn
      parameter (maxstp=100000,safety=0.9d0,errcon=3.24d-2)
      real*8 x,y,h,lam,ye,hh,lamh,yeh,u,c,rho,cc,yh,xh,y2,err
      data kstp/0/
c      kstp=1 ! PG-DEBUG
c-------------------------------------------------------- initial values
      x=x0
      nfcn=0
      call dsrdrhs(x,wrate,lam,ye,nfcn)
      ! factor=lam
      factor=1.d0   ! results do not depend on factor, but the speed of the
                    ! calculation could be improved by changing factor
      ye=ye*factor
      lam=lam/factor
      if (rderr.ne.0) return
      y=ye
      h1=hstep   ! 0.01 from block data
c      hmin=1.0d-9   ! set in block data
      eps=compeps   ! 0.01 from block data
      x1=xfinal  ! from block data
      yo=y
      eps1=1.0d-2
      xo=1.d5*mco(1)
c----------------------------------------------------------- calculation
      h=sign(h1,x1-x0)
      do 10 nstp=1,maxstp
        if((x+h-x1)*(x+h-x0).gt.0.0d0) h=x1-x
        hh=h
    1   xh=x+hh
        if (xh.eq.x) then
          write (rdluerr,*) 'error in dsrdeqn: stepsize underflow',
     &          ' at x=',x,' y=',y,' ye=',ye,' lam=',lam,' hh=',hh
          write (rdluerr,*) '  for model ',rdtag
          rderr=ibset(rderr,2)
          return
        endif
        call dsrdrhs(xh,wrate,lamh,yeh,nfcn)
        yeh=yeh*factor
        lamh=lamh/factor
        if (rderr.ne.0) return
        u=hh*lamh
        rho=lam/lamh
        c=2.0d0*y+u*((yeh*yeh+rho*ye*ye)-rho*y*y)
        if(1.0d0+u*c.lt.0.0d0) then
          hh=0.5d0*hh
          goto 1
        endif
        yh=c/(1.0d0+sqrt(1.0d0+u*c))
        cc=4.0d0*(y+u*yeh*yeh)
        y2=0.5d0*cc/(1.0d0+sqrt(1.0d0+u*cc))
        err=abs((yh-y2)/(yh*eps))
        if(err.gt.1.0d0) then
          hh=max(safety/sqrt(err),0.1d0)*hh
          goto 1
        else  ! we stepped ahead
          x=xh
          y=yh
          lam=lamh
          ye=yeh
          if (err.eq.0.d0) then
            h=5.0d0*hh
          else
            h=min(safety/sqrt(err),5.0d0)*hh
          endif
          if ((y-ye).le.cfr*ye) xf=x  ! x at freeze out
          if (kstp.ne.0) then
            if (mod(nstp-1,kstp).eq.0.and.rdlulog.gt.0)
     &          write (rdlulog,1000) x,y
1000        format (1x,'dsrdeqn> x=',f10.3:2x,'y=',e12.4)
          endif
        endif
        if((x-x1)*(x1-x0).ge.0.0d0) then
          if(abs((y-yo)/(y*eps)).le.1.0d0.or.x.ge.xo) then
            y1=y/factor
c            write (*,*) 'PG-DEBUG dsrdeqn: y1,y,factor',
c     &           y1,y,factor
            return ! with no error
          endif
          yo=y
          x1=10.0d0*x1
        endif
        if (abs(h).lt.hmin) then
          write (rdluerr,*)
     &          'error in dsrdeqn: stepsize smaller than minimum'
          write (rdluerr,*) '  for model ',rdtag
          call dsrdwres
          rderr=ibset(rderr,3)
          return
        endif
   10 continue
      write (rdluerr,*) 'error in dsrdeqn: too many steps'
      write (rdluerr,*) '  for model ',rdtag
      rderr=ibset(rderr,4)
      return
      end



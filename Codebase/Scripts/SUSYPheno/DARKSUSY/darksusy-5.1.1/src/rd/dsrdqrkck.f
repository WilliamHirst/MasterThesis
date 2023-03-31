      subroutine dsrdqrkck(f,p,wrate,x1,x2,s)
c_______________________________________________________________________
c  numerical integration with runge-kutta method.
c  input:
c    f - integrand (real,external)
c    p - parameter mass/temperature (real)
c    wrate - invariant rate (real,external)
c    x1 - lower limit (real)
c    x2 - upper limit (real)
c  output:
c    s - integral (real)
c  common:
c    'dsrdcom.h' - included common blocks
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 f,p,wrate,x1,x2,s,x
      real*8
     & eps,h,hmi,y,fx,yscal,tiny,safety,pgrow,pshrnk,htemp,xnew,
     & ytemp,yerr,err,
     & k3,k4,k5,k6,a3,a4,a5,a6,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      integer maxstp,nstp,ncalls,nfcn
      external f,wrate
      parameter (maxstp=10000,tiny=1.0d-30,eps=1.0d-6,
     & safety=0.9d0,pgrow=-0.2d0,pshrnk=-0.25d0,
     & a3=0.3d0,a4=0.6d0,a5=1.0d0,a6=0.875d0,c1=37.d0/378.d0,
     & c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,
     & dc1=c1-2825.d0/27648.d0,dc3=c3-18575.d0/48384.d0,
     & dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0,dc6=c6-0.25d0)
c-----------------------------------------------------------------------
c      write(*,*) 'dsrdqrkck called with a=',x1,'  b=',x2,'  x=',p
      hmi=(0.00001d0/maxstp)*(x2-x1)
      x=x1
      h=0.005d0*(x2-x1) ! je change 97-01-07 0.01 -> 0.005
      y=0.0d0
      ncalls=0
      do 10 nstp=1,maxstp
        fx=f(x,p,wrate)
        ncalls=ncalls+1
c        yscal=abs(y)+abs(h*fx)+tiny
        if((x+h-x2)*(x+h-x1).gt.0.0d0) h=x2-x
        htemp=h
    1   k3=f(x+a3*htemp,p,wrate)
        k4=f(x+a4*htemp,p,wrate)
        k5=f(x+a5*htemp,p,wrate)
        k6=f(x+a6*htemp,p,wrate)
        ncalls=ncalls+4
        ytemp=y+htemp*(c1*fx+c3*k3+c4*k4+c6*k6)
        yerr=htemp*(dc1*fx+dc3*k3+dc4*k4+dc5*k5+dc6*k6)
        yscal=abs(h*k5)+tiny
        err=abs(yerr/yscal)/eps
c        write (*,*) 'dsrdqrkck>',x,htemp,'htemp',ytemp,'ytemp',err,'err'
        if (err.gt.1.0d0) then
          htemp=max(0.1d0,safety*(err**pshrnk))*htemp
          xnew=x+htemp
          if (xnew.eq.x) then
            write (rdluerr,*) 'error in dsrdqrkck: stepsize underflow'
            write (rdluerr,*) '  for model ',rdtag
            rderr=ibset(rderr,5)
            return
          endif
          goto 1
        else
          if (err.eq.0.d0) then ! corr mar 95
             h=0.d0
          else
             h=min(5.0d0,safety*(err**pgrow))*htemp
          endif
          x=x+htemp
          y=ytemp
        endif
        if((x-x2)*(x2-x1).ge.0.0d0) then
          s=y
          nfcn=nfcn+ncalls
c         write (*,*) 'dsrdqrkck>',p,nstp,'stp',ncalls,'calls',nfcn,'fcn'
          return
        endif
        if(abs(h).lt.hmi) then
          write (rdluerr,*)
     &	     'error in dsrdqrkck: stepsize smaller than minimum'
          write (rdluerr,*) '  for model ',rdtag
          call dsrdwres
          call dsrdwfunc(p,wrate)
          rderr=ibset(rderr,6)
          return
        endif
   10 continue
      write (rdluerr,*) 'error in dsrdqrkck: too many steps'
      write (rdluerr,*) '  for model ',rdtag
      rderr=ibset(rderr,7)
      end

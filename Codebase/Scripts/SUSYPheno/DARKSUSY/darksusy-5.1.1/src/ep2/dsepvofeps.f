****************************************************************
***                                                          ***
*** functions which give conversions between the variables   ***
*** eps - u - v  for the positron flux calculation           ***
*** the conversion between eps and v and viceversa is done   ***
*** through a tabulation, to reset this tabulation           ***
*** reinitialize the variable vofeset                        ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-02-03                                         ***
****************************************************************

      real*8 function dsepuofeps(eps)
c this is the definition of the function u(eps)
c u = int_eps^epsmax depsp 1/(tau*b(epsp))
c NOTE: this assumes u=1/eps, change this when you implement the general
c formula
      implicit none
      real*8 eps
      dsepuofeps=1.d0/eps
      return
      end

      real*8 function dsepwofu(u)
c this is the definition of the function w(u)
      implicit none
      real*8 u,dsephofu
      dsepwofu=1.d0/u**2/dsephofu(u)
      return
      end


      real*8 function dsephofu(u)
c this is the definition of the function h(u), which is the diffusion
c coefficient divided by k0 and as a function of u
c NOTE: this assumes u=1/eps, change this when you implement the general
c formula
      implicit none
      include 'dsepcom.h'
      real*8 u
      if(dsepdiffmethod.eq.0) then
        dsephofu=u**(-alphexp)
      elseif(dsepdiffmethod.eq.1) then  
        dsephofu=1.d0+u**(-alphexp)/3.d0**alphexp
      elseif(dsepdiffmethod.eq.5) then
        if(u.ge.1.d0/4.d0) then
          dsephofu=1.d0
        else  
          dsephofu=u**(-alphexp)/4.d0**alphexp
        endif  
      else
        write(*,*) 'dsepdiffmethod not properly set' 
        write(*,*) 'dsepdiffmethod = ',dsepdiffmethod 
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end


      real*8 function dsepvofeps(eps)
c v = int_0^u dup h(up) as a functio of eps
      implicit none
c
      real*8 eps
      integer k
      real*8 epsmin,epsmax,logeps,u,dsepuofeps
      integer ndigits,vofeset
      real*8 epsabs
      common/vofesetcom/epsabs,ndigits,vofeset
      real*8 epspoints,xepsp(3001),yepsp(3001),yepsp2(3001)
     & ,xepspb(3001),yepspb(3001),yepspb2(3001)
      common/vofecom/xepsp,yepsp,yepsp2
     & ,xepspb,yepspb,yepspb2,epspoints
c
      real*8 up,low
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsephofu
c
      if(vofeset.ne.123456) then
        write(*,*) 'dsepvofeps tabulation started'
        ndigits=5
        epsabs=1.d-14     !numerical accuracy
        epsrel=10.d0**(-ndigits)
        limit=5000
        epsmin=0.01d0 
        epsmax=1.2d4
        epspoints=3.d3
        do k=0,int(epspoints)-2
          logeps=dlog(epsmin)+(dlog(epsmax)-dlog(epsmin))
     &           /(epspoints-2.d0)*k
          u=dsepuofeps(dexp(logeps)) 
          low=0.d0
          up=u
 20       call dqagse(dsephofu,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(result)/10.d0**ndigits.lt.epsabs) then
            epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 20
          elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
            epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy
          endif
          xepsp(k+1)=logeps
          yepsp(k+1)=dlog(result)
c          write(*,*) k+1,xepsp(k+1),yepsp(k+1)
        enddo
        write(*,*) 'dsepvofeps tabulation is over'
        xepsp(1)=xepsp(2)-(xepsp(3)-xepsp(2))/1.d3
        yepsp(1)=yepsp(2)
        xepsp(int(epspoints))=xepsp(int(epspoints)-1)
     &    +(xepsp(int(epspoints)-1)-xepsp(int(epspoints)-2))/1.d3
        yepsp(int(epspoints))=yepsp(int(epspoints)-1)
c        write(*,*) int(epspoints),xepsp(1),xepsp(int(epspoints))
        call dsspline(xepsp,yepsp,int(epspoints),1.d31,1.d31,yepsp2)
        do k=2,int(epspoints)-1
          xepspb(k)=yepsp(int(epspoints)-k+1)
          yepspb(k)=xepsp(int(epspoints)-k+1)
        enddo
        xepspb(1)=xepspb(2)-(xepspb(3)-xepspb(2))/1.d3
        yepspb(1)=yepspb(2)
        xepspb(int(epspoints))=xepspb(int(epspoints)-1)
     &    +(xepspb(int(epspoints)-1)-xepspb(int(epspoints)-2))/1.d3
        yepspb(int(epspoints))=yepspb(int(epspoints)-1)
        call dsspline(xepspb,yepspb,int(epspoints),1.d31,1.d31,yepspb2)
        do k=1,int(epspoints)
c          write(*,*) k,xepspb(k),yepspb(k),yepspb2(k)
        enddo
        vofeset=123456 
      endif
      if(dlog(eps).ge.xepsp(1).and.dlog(eps).le.xepsp(int(epspoints))) 
     & then
         call dssplint(xepsp,yepsp,yepsp2,int(epspoints),
     &        dlog(eps),result)
        dsepvofeps=dexp(result)
      else
        write(*,*) 'dsepvofeps called out of the defined range, eps =',
     &    eps
        write(*,*) 'epsmin, epsmax = ',dexp(xepsp(1))
     &       ,dexp(xepsp(int(epspoints)))
        write(*,*) 'program stopped'
        stop
      endif
      return
      end



      real*8 function dsepepsofv(v)
      implicit none
c
      real*8 v,result,dummy,dsepvofeps
      real*8 epspoints,xepsp(3001),yepsp(3001),yepsp2(3001)
     & ,xepspb(3001),yepspb(3001),yepspb2(3001)
      common/vofecom/xepsp,yepsp,yepsp2
     & ,xepspb,yepspb,yepspb2,epspoints
      integer ndigits,vofeset
      real*8 epsabs
      common/vofesetcom/epsabs,ndigits,vofeset
c
      if(vofeset.ne.123456) then
        dummy=dsepvofeps(1.d0)
      endif
      if(dlog(v).ge.xepspb(1).and.dlog(v).le.xepspb(int(epspoints))) 
     &  then
        call dssplint(xepspb,yepspb,yepspb2,int(epspoints),
     & dlog(v),result)
        dsepepsofv=dexp(result)
      else
        write(*,*) 'dsepepsofv called out of the defined range, v =',v
        write(*,*) 'vmin, vmax = ',dexp(xepspb(1))
     &       ,dexp(xepspb(int(epspoints)))
        write(*,*) 'program stopped'
        stop
      endif
      return
      end

**********************************************************************
*** function called in dspbtd15x
*** it gives the antiproton diffusion time in units of 10^15 sec
*** it assumes the diffusion model in:
***   bergstrom, edsjo & ullio, ajp 526 (1999) 215
***   but with the DC-like setup as in moskalenko et al.
***      ApJ 565 (2002) 280
*** inputs:
***     tp - antiproton kinetic energy (gev)
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbtd15beum(tp)
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
c      integer i
      real*8 tp,pp,ee,kpc,cc,partial,total,dspbkdiffm,dspbsigmavpbar
      real*8 abserr,alist,blist,elist,epsabs,epsrel,rlist,
     & result,dbesj0,dbesj1,beta,relprec
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dspbbeuparm
      parameter (kpc=3.08567802d0,cc=c_light/1.d5)
      real*8 jzero(nzerojk)
      integer k,s,ii
      integer incr
      parameter(incr=99)
      real*8 dsbessjw
      real*8 sumscheck,seriesvec(100)
      real*8 nusk,Jklocal,Jkplus1squared
      integer js,jscheck
      integer inipbcheck
      common/pbinicom/inipbcheck
      if(inipbcheck.ne.123456) then
        k=0 
        do ii=1,nzerojk
          storage(k,ii,1)=0.d0
        enddo
        relprec=1.d-3
        inipbcheck=123456
      endif  
      jscheck=15
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      epsabs=1.d-6     !numerical accuracy
      epsrel=1.d-3
      limit=5000
      pp=dsqrt(2*m_p*tp+tp**2)
      ee=tp+m_p
      beta=pp/ee
      dh=dspbkdiffm(beta,pp,1)
      dg=dspbkdiffm(beta,pp,2)
      axsec=dspbsigmavpbar(ee)
      partial=0.d0
      total=0.d0
      k=0
      do s=1,nzerojk
 10     nusk=storage(k,s,1)
        if(nusk.lt.1.d-16) then
          CALL DBZEJY(dble(k),s+incr,1,1.d-16,jzero)
          do ii=s,s+incr,1
            storage(k,ii,1)=jzero(ii)
            storage(k,ii,2)=dsbessjw(k,jzero(ii)*pbr0/pbrh)
            storage(k,ii,3)=(dsbessjw(k+1,jzero(ii)))**2
          enddo
          goto 10
        endif
        zero=nusk
        Jklocal=storage(k,s,2)
        Jkplus1squared=storage(k,s,3)
c  lambdag,gamh,lambdah in 10**(-21) cm^-1
        lambdag=dsqrt(axsec*pbng/dg/100.d0+(zero/pbrh/kpc)**2)
        gamh=pbcvel/2.d0/dh/10.d0
        lambdah=dsqrt(gamh**2+axsec*pbnh/dh/100.d0
     &                +(zero/pbrh/kpc)**2)
        partial=0.d0
        call dqagse(dspbbeuparm,0.d0,1.d-3/pbrh,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        partial=partial+result
        write(*,*) 'izero,partial,result ',s,partial,result
        call dqagse(dspbbeuparm,0.d0,1.d0,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        partial=partial+result
        write(*,*) 'izero,partial,result ',s,partial,result
        result=partial
c        call dqagse(dspbbeuparm,0.d0,1.d0,epsabs,epsrel,limit,
c     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c        partial=result*dbesj0(zero*pbr0/pbrh)/(dbesj1(zero))**2
c     &       *2.d0/(pbrh*kpc)**2
        partial=result*Jklocal/Jkplus1squared
     &       *2.d0/(pbrh*kpc)**2
        total=total+partial
        write(*,*) 'izero,totpar,total ',s,total,partial
        if(s.le.jscheck) then
          seriesvec(s)=partial
        else
          do js=2,jscheck 
            seriesvec(js-1)=seriesvec(js)
          enddo  
          seriesvec(jscheck)=partial
          sumscheck=0.d0 
          do js=1,jscheck 
            sumscheck=sumscheck+seriesvec(js)
          enddo 
          write(*,*) 'sumscheck, total = ',sumscheck,total
          if(dabs(sumscheck/total).lt.relprec) then
            goto 20
          endif  
        endif
      enddo
      write(*,*) 'storage matrix has index s too small, k, s = ',k,s
      write(*,*) 'program stopped'
      stop 
c  dspbtd15beum in 10^15 sec 
 20   dspbtd15beum = total
      return
      end









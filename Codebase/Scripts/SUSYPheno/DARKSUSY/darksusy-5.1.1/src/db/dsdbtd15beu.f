**********************************************************************
*** function called in dsdbtd15x
*** it gives the antiproton diffusion time in units of 10^15 sec
*** it assumes the diffusion model in:
***   bergstrom, edsjo & ullio, ajp 526 (1999) 215
*** inputs:
***     td - kinetic energy per nucleon (gev)
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
*** modified: 04-01-22 (pu)
**********************************************************************

      real*8 function dsdbtd15beu(td)
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
      integer i
      real*8 td,nucleon,ee,pp,rig
      real*8 kpc,cc,partial,total,dspbkdiffm,dspbkdiff,
     & dsdbsigmavdbar
      real*8 abserr,alist,blist,elist,epsrel,rlist,
     & result,beta,dsbessjw
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dspbbeuparm
      parameter (kpc=3.08567802d0,cc=c_light/1.d5)
      real*8 jzero(nzerojk)
      integer k,s,ii
      integer incr
      parameter(incr=99)
      real*8 sumscheck,seriesvec(100)
      real*8 nusk,Jklocal,Jkplus1squared,low,up,relprec
      integer js,jscheck
      real*8 epsabs1,epsabs2,epsabs3
      integer beuset,ndigits
      common/dbbeusetcom/ epsabs1,epsabs2,epsabs3,beuset,ndigits
      if(beuset.ne.123456) then
         ndigits=5
         epsabs1=1.d-10          !numerical accuracy
         epsabs2=1.d-10          !numerical accuracy
         epsabs3=1.d-10          !numerical accuracy
         beuset=123456
      endif
      relprec=1.d-4
      epsrel=10.d0**(-ndigits) 
      jscheck=15
      limit=5000
      nucleon=2.d0
      ee=nucleon*td+m_d
      pp=dsqrt(dabs(ee**2-m_d**2))
c      rig=pp/nucleon
      rig=pp                    ! bug fixed PU 111104
      if(pbpropmodel.eq.2) then
        dh=dspbkdiff(rig,1)
        dg=dspbkdiff(rig,2)
      elseif(pbpropmodel.eq.3) then
        beta=pp/ee
        dh=dspbkdiffm(beta,rig,1)
        dg=dspbkdiffm(beta,rig,2)
      else
        write(*,*) 'dsdbtd15beu called with wrong pbpropmodel'
        write(*,*) 'pbpropmodel = ',pbpropmodel
        write(*,*) 'program stopped'
        stop
      endif
      axsec=dsdbsigmavdbar(ee)
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
c split the integration to allow for the option to exclude a cylinder
c corresponding to a central portion of the halo
        if(pbrcy.gt.0.d0.and.pbrcy.lt.pbrh) then
          low=0.d0
          up=pbrcy/pbrh
 21       call dqagse(dspbbeuparm,low,up,epsabs1,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(result)/10.d0**ndigits.lt.epsabs1) then
            epsabs1=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 21
          elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs1) then
            epsabs1=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
          partial=result
          low=up
          up=1.d0
 22       call dqagse(dspbbeuparm,low,up,epsabs2,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(result)/10.d0**ndigits.lt.epsabs2) then
            epsabs2=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 22
          elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs2) then
            epsabs2=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
          partial=partial+result
          result=partial
        else   
 23        call dqagse(dspbbeuparm,0.d0,1.d0,epsabs3,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(result)/10.d0**ndigits.lt.epsabs3) then
            epsabs3=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 23
          elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs3) then
            epsabs3=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
       endif
        partial=result*Jklocal/Jkplus1squared
     &       *2.d0/(pbrh*kpc)**2
        total=total+partial
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
          if(dabs(sumscheck/total).lt.relprec) then
            goto 20
          endif  
        endif
      enddo
      write(*,*) 'storage matrix has index s too small, k, s = ',k,s
      write(*,*) 'program stopped'
      stop 
c  dsdbtd15beu in 10^15 sec 
 20   dsdbtd15beu = total
      return
      end














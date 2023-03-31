**********************************************************************
*** function called in dspbtd15beum
*** it is integrated in the cylindrical coordinate r from 0 to 1 
*** (linear change of variables such that r=1 corresponds to r=pbrh 
*** (radial extent of the diffusion box))
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
*** modified: 04-01-22 (pu)
**********************************************************************

      real*8 function dspbbeuparm(r)
      implicit none
      include 'dspbcom.h'
      include 'dspbprivate.h'
      real*8 zlower,up,low
      real*8 r,rr,rebigh,rebiggs,rebiggc,kpc,dbesj0
      real*8 abserr,alist,blist,elist,epsrel,rlist
     & ,resultgs,resultgc,resulth
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      parameter(kpc=3.08567802d0)
      external dspbbeuparh,dspbbeupargs,dspbbeupargc
      common /czint/rr
      real*8 epsabss,epsabsc,epsabsh
      integer parmset,ndigits
      common/parmsetcom/epsabss,epsabsc,epsabsh,parmset,ndigits
      rr=r
      if(parmset.ne.123456) then
         ndigits=5
         epsabss=1.d-10          !numerical accuracy
         epsabsc=1.d-10          !numerical accuracy
         epsabsh=1.d-10          !numerical accuracy
         parmset=123456
      endif
      epsrel=10.d0**(-ndigits)
      limit=5000
c  rebigh, rebiggs, rebiggc in 10^6 cm sec^-1
      rebigh=1.d0/dsinh(lambdag*pbhg*kpc)
     &  /(dg*lambdag+(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &  /dtanh(lambdag*pbhg*kpc))
      rebiggs=1.d0
     &  /(dg*lambdag+(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &  /dtanh(lambdag*pbhg*kpc)) 
      rebiggs=rebiggs/dg/lambdag*(dh*gamh
     &     +dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
      rebiggc=1.d0
     &  /(dg*lambdag+(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &  /dtanh(lambdag*pbhg*kpc))
c allow for the option to exclude a cylinder corresponding to a central 
c portion of the halo
      if(r*pbrh.lt.pbrcy) then
        zlower=pbzcy
        if(zlower.lt.pbhg.and.zlower.ge.0.d0) then
          up=pbhg/pbhh
          low=zlower/pbhh
 10       call dqagseb(dspbbeupargs,low,up,epsabss,epsrel,limit,
     &      resultgs,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(resultgs)/10.d0**ndigits.lt.epsabss) then
            epsabss=dabs(resultgs)/10.d0**(ndigits+2) !too low accuracy
            goto 10
          elseif(dabs(resultgs)/10.d0**(ndigits+2).gt.epsabss) then
            epsabss=dabs(resultgs)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
 11       call dqagseb(dspbbeupargc,low,up,epsabsc,epsrel,limit,
     &      resultgc,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(resultgc)/10.d0**ndigits.lt.epsabsc) then
            epsabsc=dabs(resultgc)/10.d0**(ndigits+2) !too low accuracy
            goto 11
          elseif(dabs(resultgc)/10.d0**(ndigits+2).gt.epsabsc) then
            epsabsc=dabs(resultgc)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
          low=up
          up=1.d0
 12       call dqagseb(dspbbeuparh,low,up,epsabsh,epsrel,limit,
     &      resulth,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(resulth)/10.d0**ndigits.lt.epsabsh) then
            epsabsh=dabs(resulth)/10.d0**(ndigits+2) !too low accuracy
            goto 12
          elseif(dabs(resulth)/10.d0**(ndigits+2).gt.epsabsh) then
            epsabsh=dabs(resulth)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
        elseif(zlower.le.pbhh.and.zlower.ge.pbhg) then
          resultgs=0.d0
          resultgc=0.d0
          low=zlower/pbhh
          up=1.d0
 13       call dqagseb(dspbbeuparh,low,up,epsabsh,epsrel,limit,
     &      resulth,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(resulth)/10.d0**ndigits.lt.epsabsh) then
            epsabsh=dabs(resulth)/10.d0**(ndigits+2) !too low accuracy
            goto 13
          elseif(dabs(resulth)/10.d0**(ndigits+2).gt.epsabsh) then
            epsabsh=dabs(resulth)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
       else
          write(*,*) 'wrong setup for zlower in dspbbeuparm1'
          write(*,*) 'zlower, pbgh, pbhh = ',zlower,pbhg,pbhh 
          write(*,*) 'program stopped'
          stop
        endif
      else   
 20     call dqagseb(dspbbeupargs,0.d0,pbhg/pbhh,epsabss,epsrel,limit,
     &    resultgs,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
          if(dabs(resultgs)/10.d0**ndigits.lt.epsabss) then
            epsabss=dabs(resultgs)/10.d0**(ndigits+2) !too low accuracy
            goto 20
          elseif(dabs(resultgs)/10.d0**(ndigits+2).gt.epsabsc) then
            epsabss=dabs(resultgs)/10.d0**(ndigits+2)  !too high accuracy 
          endif  
 21     call dqagseb(dspbbeupargc,0.d0,pbhg/pbhh,epsabsc,epsrel,limit,
     &    resultgc,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        if(dabs(resultgc)/10.d0**ndigits.lt.epsabsc) then
          epsabsc=dabs(resultgc)/10.d0**(ndigits+2) !too low accuracy
          goto 21
        elseif(dabs(resultgc)/10.d0**(ndigits+2).gt.epsabsc) then
          epsabsc=dabs(resultgc)/10.d0**(ndigits+2)  !too high accuracy 
        endif  
 22     call dqagseb(dspbbeuparh,pbhg/pbhh,1.d0,epsabsh,epsrel,limit,
     &    resulth,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        if(dabs(resulth)/10.d0**ndigits.lt.epsabsh) then
          epsabsh=dabs(resulth)/10.d0**(ndigits+2) !too low accuracy
          goto 22
        elseif(dabs(resulth)/10.d0**(ndigits+2).gt.epsabsh) then
          epsabsh=dabs(resulth)/10.d0**(ndigits+2)  !too high accuracy 
        endif  
      endif
      dspbbeuparm=(pbrh*kpc)**2*r*dbesj0(zero*r)
     &   *(resulth*rebigh+resultgs*rebiggs+resultgc*rebiggc)
      return
      end

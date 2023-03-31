**********************************************************************
*** function called in dsdbtd15x
*** it gives the antiproton diffusion time in units of 10^15 sec
*** it assumes the diffusion model in:
***   bergstrom, edsjo & ullio, ajp 526 (1999) 215
***   but with the DC-like setup as in moskalenko et al.
***      ApJ 565 (2002) 280
*** inputs:
***     td - kinetic energy per nucleon (gev)
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dsdbtd15beum(td)
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
      integer i
      real*8 td,nucleon,ee,pp,rig
      real*8 kpc,cc,partial,total,dspbkdiffm,dsdbsigmavdbar
      real*8 abserr,alist,blist,elist,epsabs,epsrel,rlist,
     & result,dbesj0,dbesj1,beta
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dspbbeuparm
      parameter (kpc=3.08567802d0,cc=c_light/1.d5)
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      limit=5000
      nucleon=2.d0
      ee=nucleon*td+m_d
      pp=dsqrt(dabs(ee**2-m_d**2))
      beta=pp/ee
      rig=pp/nucleon
      dh=dspbkdiffm(beta,rig,1)
      dg=dspbkdiffm(beta,rig,2)
      axsec=dsdbsigmavdbar(ee)
      partial=0.d0
      total=0.d0
      do i=1,pbnzero
        zero=zerovec(i)
c  lambdag,gamh,lambdah in 10**(-21) cm^-1
        lambdag=dsqrt(axsec*pbng/dg/100.d0+(zero/pbrh/kpc)**2)
        gamh=pbcvel/2.d0/dh/10.d0
        lambdah=dsqrt(gamh**2+axsec*pbnh/dh/100.d0
     &                +(zero/pbrh/kpc)**2)
        call dqagse(dspbbeuparm,0.d0,1.d0,epsabs,epsrel,limit,result,
     &      abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        partial=result*dbesj0(zero*pbr0/pbrh)/(dbesj1(zero))**2
     &       *2.d0/(pbrh*kpc)**2
        total=total+partial
      enddo
c  dsdbtd15beu in 10^15 sec 
      dsdbtd15beum = total
      return
      end








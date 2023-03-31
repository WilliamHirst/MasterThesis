**********************************************************************
*** function called in dspbtd15x
*** it gives the antiproton diffusion time in units of 10^15 sec
*** it assumes the diffusion model in:
***   chardonnet et al., phys. lett. b384 (1996) 161
***   bottino et al.,  phys. rev. d58 (1998) 123503
*** inputs:
***     tp - antiproton kinetic energy (gev)
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbtd15char(tp)
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
      integer i
      real*8 tp,kpc,cc,partial,total,
     &  dspbkdiff,dspbsigmavpbar,coeffa,coeffb,pp,lpbar,dbesj0,dbesj1
      real*8 abserr,alist,blist,elist,epsabs,epsrel,rlist,
     & result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dspbcharpar1
      parameter(kpc=3.08567802d0,cc=c_light/1.d5)
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      limit=5000
      pp=dsqrt(2*m_p*tp+tp**2)
      dh=dspbkdiff(pp,1)
      lpbar=dspbsigmavpbar(tp+m_p)*pbng
      total=0.d0
      partial=0.d0
      do i=1,pbnzero
        zero=zerovec(i)
        call dqagse(dspbcharpar1,0.d0,1.d0,epsabs,epsrel,limit,result,
     &    abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        coeffa=2.d0*(pbhg*kpc)*lpbar/100.d0 ! 100.d0 due to units
        coeffa=coeffa+2.d0*dh*zero/(pbrh*kpc)/
     &           dtanh(pbhh/pbrh*zero)
        coeffb=4.d0*dbesj0(pbr0/pbrh*zero)/(dbesj1(zero))**2*
     &           (pbhh*kpc)
        partial=result*coeffb/coeffa
        total=total+partial
      enddo
c  dspbtd15char in 10^15 sec 
      dspbtd15char=total
      return
      end





************************************************************************
*** positron propagation routines.
*** author: e.a. baltz (eabaltz@astron.berkeley.edu)
*** modified slightly by joakim edsjo (edsjo@physto.se)
*** date: jun-02-98
*** modified: jun-09-98, 2002-11-19
***   2004-01-26 (better r integration)
************************************************************************


************************************************************************
      subroutine dsepmake_tables
***   creates a table of i(delta v) versus delta v.
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 deltav,lndeltav,deltar,dsbessei0
      real*8 rprime,rmin,rmax,integrand,dsepf,dsepimage_sum,
     &  lrmin,lrmax
      real*8 lne,egev,u,dsepideltavint,dsf_int2,imsum
      integer i
      external dsepideltavint

      real*8 abserr,alist,blist,elist,epsrel,rlist,result,epsabs
      integer ier,iord,last,limit,neval,ndigits
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)

      real*8 dv
      common/epint/dv

      k0tau = k27*tau16
      ametric=1.0d0/3.0d0**alphexp
      av=ametric/(1.0d0-alphexp)
      do 30 i=1,22001
         lndeltav=real(i-1)*0.001d0-19.0d0
         deltav=exp(lndeltav)
         dv=4.0d0*k0tau*deltav


c...First calculate the sum of image charges
c...JE change 090109: always sum, don't put to 1 for small dv any more
         imsum=dsepimage_sum(deltav)
         imsum=max(0.0d0,imsum) ! To avoid numerical problems, 021119
      
         if (imsum.gt.0.0d0) then  ! do the radial integration

c...If dv is small, skip the r integration as it is not needed
           if(dv.lt.0.1d0) then
             table(i)=dsepf(r_e)
             goto 50
           endif

c...Integrate over r with dqagse integration
           limit=5000
           ndigits=2
           epsabs=1.d-14     !numerical accuracy
           epsrel=10.d0**(-ndigits)

           rmin=1.0d-6  ! to make sure to catch the cusp at the center
           rmax=r_e+rwid*sqrt(dv)

           lrmin=log(rmin)
           lrmax=log(rmax)

c...full integration in one go
 13        call dqagseb(dsepideltavint,lrmin,lrmax,epsabs,epsrel,
     &       limit,
     &       result,abserr,neval,ier,alist,blist,rlist,elist,iord,
     &       last)
           if(dabs(result)/10.d0**ndigits.lt.epsabs) then
             epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
             goto 13
           endif
           integrand=result/1.d10    ! 1.d10 is from rescaling in dsephalodens2

           table(i)=integrand

         else  ! sum of image charges=0, put table(i) to zero as well
           table(i)=0.0d0
         endif

 50      table(i)=table(i)*imsum
         table(i)=max(table(i),0.0d0)  ! To avoid numerical problems, 021119

 30   continue

      return
      end

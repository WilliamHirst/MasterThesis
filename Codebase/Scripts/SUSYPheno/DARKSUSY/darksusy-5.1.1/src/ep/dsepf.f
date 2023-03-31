************************************************************************
*** This is the average of the halo density squared from
*** z=-l_h to z=+l_h.
*** It is the function f(r) given in Eq. (20) in Baltz & Edsjo, 
*** PRD 59(1999)023511, except that g(r) here is the NORMALIZED
*** halo density. Consequently, n_c=1.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2000-09-03
************************************************************************

************************************************************************
      real*8 function dsepf(r)
************************************************************************
      implicit none
      include 'dsepcom.h'
      real*8 r,dsephalodens2,dsf_int,dseploghalodens2
      external dsephalodens2,dseploghalodens2

      real*8 abserr,alist,blist,elist,epsrel,rlist,result,epsabs
      integer ier,iord,last,limit,neval,ndigits
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)

c----------------------------------------------------------------------

      r_gc = r   ! put in common block that dsephalodens2 can see

c...Integrate over z
      limit=5000
      ndigits=3 
      epsabs=1.d-14     !numerical accuracy
      epsrel=10.d0**(-ndigits)

c 10   call dqagse(dsephalodens2,0.0d0,l_h,epsabs,epsrel,limit,
c     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
 10   call dqagse(dseploghalodens2,log(1.0d-6),log(l_h),
     &      epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 10
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dsepf=result/l_h/1.d10    ! 1.d10 is from rescaling in dsephalodens2
      

c...Old simpler simplex integration
c      dsepf=1.0d0/l_h*
c     &  (dsf_int(dsephalodens2,0.0d0,l_h/100.0d0,0.01d0)
c     &  +dsf_int(dsephalodens2,l_h/100.0d0,l_h,0.01d0))
c     &  /1.d0     ! 1.d10 is from rescaling in dsephalodens2

      return
      end




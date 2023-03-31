      subroutine dsddeta(vmin,t,eta)
c_______________________________________________________________________
c  eta function entering the differential rate: eta = \int {f(v)/v} d^3 v
c
c  Truncated Maxwellian. 
c
c  input:
c    vmin : minimum velocity to deposit energy e, in km/s 
c           vmin=sqrt(M*E/2/mu^2)
c    t : time, in fraction of the year
c  output:
c    eta : in (km/s)^{-1}
c  authors: paolo gondolo (paolo@physics.utah.edu) 2004
c           piero ullio (ullio@sissa.it) 2004
c
c  2010/06/17 [C Savage]: fixed error in iso formula;
c                         added z < y case
c=======================================================================
      implicit none
      include 'dshmcom.h'
      include 'dsge.h'
      real*8 vmin,t,eta,etatmp,v0bar,x,y,z,vearth,erf,Nesc,vmax
ccc
      DOUBLE PRECISION ABSERR,ALIST,BLIST,ELIST,EPSREL,RLIST
     & ,RESULT
      INTEGER IER,IORD,LAST,LIMIT,NEVAL
      DIMENSION ALIST(5000),BLIST(5000),ELIST(5000),IORD(5000),
     & RLIST(5000)
      external dshmuDF
ccc
      real*8 epsabs
      integer dRdQsiset,ndigits
      common/dRdQsicom/epsabs,dRdQsiset,ndigits

      call dsddvearth(t,vearth)

      if (veldf.eq.'gauss'.or.veldf.eq.'iso') then
         v0bar=sqrt(2.d0/3.d0)*vd_3d
         x=vmin/v0bar
         y=vearth/v0bar
         z=vgalesc/v0bar
         Nesc = erf(z)-2.d0*z*exp(-z**2)/sqrt(pi)
         if (x.le.abs(z-y)) then
            if (y.le.z) then
              etatmp = erf(x+y)-erf(x-y)-4/sqrt(pi)*y*exp(-z**2)
            else
              etatmp = 2.d0*Nesc
            end if
         else if (x.lt.z+y) then
            etatmp = erf(z)-erf(x-y)-2/sqrt(pi)*(z+y-x)*exp(-z**2)
         else
            etatmp = 0.d0
         endif
         eta = 0.5d0/vearth/Nesc*etatmp
      else
         if(dRdQsiset.ne.123456) then
            ndigits=3
            epsabs=1.d-10       !numerical accuracy
            dRdQsiset=123456
         endif
         epsrel=10.d0**(-ndigits)
         limit=5000
         vmax=vgalesc+vearth
         if(vmin.gt.vmax) then
            write(*,*) 'vmin larger than vmax: did you rember to set'
            write(*,*) 'vgalesc?   vmin, vmax = ',vmin,vmax
            eta=0.d0
            return
         endif
         v_obs=vearth
 10      call DQAGSE(dshmuDF,vmin,vmax,EPSABS,EPSREL,LIMIT,RESULT,
     &        ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)   
         if(dabs(result)/10.d0**ndigits.lt.epsabs) then
            epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
            goto 10
         endif 
         eta = result

      endif

      return
      end

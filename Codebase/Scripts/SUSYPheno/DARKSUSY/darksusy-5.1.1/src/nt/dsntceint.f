
***********************************************************************
*** l.b. and j.e. 1999-04-06
*** auxiliary function for r-integration
*** input: radius in centimeters
*** output: integrand in cm^-1 s^-1
***********************************************************************

      real*8 function dsntceint(r,foveru)
      implicit none
      real*8 mx,sigsi,mu,muplus,muminus,ma,rx,eps
      real*8 v,res,dsntceint2,r,max,dsntearthvesc,
     &  umin,umax,foveru,sigsd
      integer imass,vtype
      external foveru
      common/ntdkint/mx,ma,sigsi,sigsd,rx,vtype
      common/ntdkint2/imass
      external dsntceint2
      include 'dsntdkcom.h'
      include 'dshmcom.h'

      rx=r
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
c...begin velocity integration
c...
c...determine integration limits
      v=dsntearthvesc(r/100.0d0)   ! escape velocity at r

c...For a general velocity distribution
      umin=0.d0                  ! lower velocity limit
      umax=sqrt(mu/muminus**2)*v ! upper velocity limit
c...Note, we should only allow such velocities that we scatter to
c...velocities lower than the escape velocity. This means that we have to
c...make sure that (mu/muplus^2 > u^2 / (u^2 + v^2 )). This is the same as
c...the condition umax above.

      if (vtype.eq.2) then ! DK velocity distribution
c...velocity limits coming from the actual form of the dk distribution
c...note: change this for a different velocity distribution!!!
        eps=0.18377d0*dklambda
        umin=max(v_earth*sqrt(1.0d0-eps),umin)
        umax=min(v_earth*sqrt(3.0d0-1.0d0/2.6d0),umax)
      endif

c      write(*,*) 'umin=',umin,'  umax=',umax
      if (umin.lt.umax) then
        call dshiprecint2(dsntceint2,foveru,umin,umax,res)
        res=res*1.0d5    ! km/s -> cm/s from u integration
        dsntceint=res*r*r*4.*3.141592  ! assume spherical earth
      else
        dsntceint=0.0d0
      endif

c      if (umax.gt.30.0d0) then
c        write(*,*) 'umax=',umax,'  res=',res
c      endif

      return
      end

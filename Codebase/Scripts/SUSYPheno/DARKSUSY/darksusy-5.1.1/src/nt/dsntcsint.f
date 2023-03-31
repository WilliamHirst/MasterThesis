
***********************************************************************
*** l.b. and j.e. 1999-04-06
*** auxiliary function for r-integration
*** input: radius in centimeters
*** output: integrand in cm^-1 s^-1
*** Adapted for the Sun by J. Edsjo, 2003-11-26
***********************************************************************

      real*8 function dsntcsint(r,foveru)
      implicit none
      real*8 mx,sigsi,mu,muplus,muminus,ma,rx
      real*8 v,res,dsntcsint2,r,max,dsntsunvesc,
     &  umin,umax,foveru,sigsd,vp
      integer iel,vtype
      external foveru
      common/ntdkint/mx,ma,sigsi,sigsd,rx,vtype
      common/ntdkint2/iel
      external dsntcsint2
      include 'dsntdkcom.h'
      include 'dsntcom.h'

      rx=r
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
c...begin velocity integration
c...
c...determine integration limits
      v=dsntsunvesc(r/100.0d0)   ! escape velocity at r

c...Now use a different maximal velocity after the scatter to allow
c...for cutting away WIMPs that reach too far out (to Jupiter).
      vp=sqrt(v**2-veout**2)

c...general velocity limits
      umin=0.d0                  ! lower velocity limit
c      umax=sqrt(mu/muminus**2)*v ! upper velocity limit
c...Correct umax for a differnt choice of maximum velocity
c...veout is typically set to 0 for standard treatment or to
c...the escape veolocity at Jupiter for a more conservative estimate
c...
c...The 10000.d0 maximum is put to avoid numerical problems when mx and
c...ma matches. In principle, that should be set to the maximum velocity
c...for which f(u) is non-zero, but 10000.d0 is a safe choice here.
      if (abs(muminus).gt.1.d-10) then
        umax=sqrt(max(vp**2*mu/muminus**2-veout**2,0.d0))
        umax=min(umax,10000.d0)
      else
        umax=10000.d0
      endif
c...Note, we should only allow such velocities that we scatter to
c...velocities lower than the escape velocity. This means that we have to
c...make sure that (mu/muplus^2 > u^2 / (u^2 + v^2 )). This is the same as
c...the condition umax above.

      if (umin.lt.umax) then
        call dshiprecint2(dsntcsint2,foveru,umin,umax,res)
        res=res*1.0d5    ! km/s -> cm/s from u integration
        dsntcsint=res*r*r*4.*3.141592  ! assume spherical sun
      else
        dsntcsint=0.0d0
      endif

      return
      end

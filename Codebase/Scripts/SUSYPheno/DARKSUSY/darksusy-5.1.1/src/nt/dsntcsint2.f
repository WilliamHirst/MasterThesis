
c...
c...auxiliary function for inner integrand
c...input: velocity relaitve to sun in km/s
c...output: integrand in cm^-4
c...We here follow the analysis in Gould, ApJ 321 (1987) 571 and more
c...specifically the more general expressions in appendix A.
      real*8 function dsntcsint2(u,foveru)
      implicit none
      include 'dssun.h'
      include 'dsmpconst.h'
      include 'dsntcom.h'

      real*8 u,mu,mx,sigsi,sigsd,muplus,muminus,siga,ma
      real*8 r,ni,e0,sla1,sla2,sla3,v,dsntsunvesc,vp,
     &  dsntsundenscomp,foveru
      integer iel,vtype
      external foveru
      common/ntdkint/mx,ma,sigsi,sigsd,r,vtype
      common/ntdkint2/iel
      siga=sigsi*sdaa(iel)**2*(mx*ma/(mx+ma))**2/(mx*m_p/(mx+m_p))**2
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
      v=dsntsunvesc(r/100.0d0)   ! escape velocity at r
c...now correct escape velocity so that we demand a slightly
c...lower escape velocity after to make sure it does not reach Jupiter
      vp=sqrt(v**2-veout**2)

      ni=dsntsundenscomp(r/100.0d0,iel) ! target number density in sun

      if (iel.eq.1) then  ! no form-factor suppression for Hydrogen        

        siga=siga+sigsd  ! add spin-dependent for Hydrogen
c...We below allow for a lower maximal velocity after scatter
c...to choose only scatters that reach out to Jupiter
        sla3=siga*ni*((vp/c_light)**2
     &    -muminus**2/mu*((u/c_light)**2+(veout/c_light)**2)) ! modified gould (2.13)
        dsntcsint2=foveru(u)*sla3 ! gould (2.13) + (2.8)

      else ! heavier eleements
c...Below we assume an exponential form factor, which lets us
c...integrate over the energy loss analytically. For general form factors
c...we would need to generalize the expressions below with an integral
c...over the energy loss (or q^2)
        e0=3./2./ma*0.038938 ! gives units of gev; see gould (a8)
        e0=e0/(0.91*ma**(0.33333)+0.3)**2
c...In the expressions below, we have modified the expressions in Gould
c...(or if you prefer Lundberg & Edsjo, 2004) to allow for scatterings 
c...not only to velocities less than the escape velocity, but
c...instead require a more conservative value from requiring that
c...the WIMP does not reach out further than to a distance where the
c...escape velocity is veout. (veout=0 gives us back the old 
c...Gould expressions).
        sla1=exp(-mx*((u/c_light)**2+(veout/c_light)**2)/2./e0)  ! gould (A6)
        sla2=exp(-mu/muplus**2*mx*((v/c_light)**2+(u/c_light)**2)/2./e0) ! gould (A6)
        sla3=siga*ni*muplus**2*2.*e0/mu/mx ! gould (A6)
        dsntcsint2=foveru(u)*sla3*(sla1-sla2)  ! gould (A6) + (2.8)
c      write(*,*) 'u=',u,'  dsntcsint2=',dsntcsint2,sla3,sla1,sla2
c      write(*,*) 'foveru = ',dsntdkfoveru(u)

      endif


c...we have two factors of c missing, add them to get units cm^-4
      dsntcsint2=dsntcsint2*(3.0d10)**2
c...dsntcsint2 is now equal to f(u)/u * w Omega_v^-(w), i.e.
c...f(u)/u * (A6) in Gould.

c...Note, in principle we should have a theta function here to allow
c...only such scatterings that go to a velocity lower than the escape
c...velocity, theta(mu/muplus^2 - u^2/(u^2+v^2)),  however, this
c...is taken care of in the limits for the u integration and is thus not
c...needed here.
c      if (u**2/(u**2+v**2).gt.mu/muplus**2) then
c        dsntcsint2=0.0d0
c      endif

      return
      end

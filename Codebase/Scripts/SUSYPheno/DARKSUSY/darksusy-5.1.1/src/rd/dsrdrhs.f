      subroutine dsrdrhs(x,wrate,lambda,yeq,nfcn)
c_______________________________________________________________________
c  adimensional annihilation rate lambda in the boltzmann equation
c    y' = -lambda (y**2-yeq**2) and equilibrium dm density in units
c    of the entropy density.
c  input:
c    x - mass/temperature (real)
c    wrate - invariant annihilation rate (real)
c  output:
c    lambda - adimensional parameter in the evolution equation (real)
c    yeq - equilibrium number/entropy densities (real)
c    nfcn - number of calls to wrate (integer)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses qrkck or dgadap.
c  called by dsrdeqn.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996
c  modified: joakim edsjo (edsjo@physto.se) 98-04-28
c    bug fix 98-04-28: y_eq was too small by a factor of 2 (je)
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 dsrdthav
      real*8 x,wrate,lambda,yeq,gentr,
     &  k2,dsbessek2
      integer nfcn,i
      external wrate
      real*8 t,grav,wav,sqrtg,eqcnst,genergy,gpressure
c     grav=mplanck*sqrt(pi/45), eqcnst=45/(4*pi^4)
      parameter (grav=3.2262726d18)
      parameter (eqcnst=0.1154923d0)
      integer k,kk

c--------------------------------------------------------------------dof
      t=mco(1)/x
      call dsrddof(t,sqrtg,gentr,genergy,gpressure)
c      write (*,'(1x,a,1x,e14.8,1x,i3,1x,i3,2(1x,e14.8)))') 
c     & 'dsrdrhs: t,klo,khi,sqrtg,gentr',
c     & t,klo,khi,sqrtg,gentr

c-------------- thermally averaged cross section times relative velocity

      wav=dsrdthav(x,wrate)
      if (rderr.ne.0) return

c---------------------------------------------------------------- lambda

      lambda=grav*sqrtg*wav*mco(1)/(x**2)

c--------------------------------------------------- equilibrium density

      yeq=0.d0

      do i=1,nco
        k2=exp(-x*mco(i)/mco(1))*dsbessek2(x*mco(i)/mco(1))
        yeq=yeq+(mco(i)/mco(1))**2*mdof(i)*k2
      enddo
      yeq=eqcnst*yeq*x**2/gentr
      
      return
      end

      function dsrdfunc(u,x,wrate)
c_______________________________________________________________________
c  invariant annihilation rate times thermal distribution.
c  when integrated over u, the effective thermal average times
c  m_chi^2 is obtained.
c  input:
c    u - integration variable (real)
c    x - mass/temperature (real)
c    wrate - invariant annihilation rate (real, external)
c  common:
c    'dsrdcom.h' - included common blocks
c  called by dsrdrhs, wirate, dsrdwintrp.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996
c  modified: joakim edsjo (edsjo@physto.se) 98-04-29
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 dsrdfunc,x,p,y,u,u2,www,wrate,
     &  denom,ek1,ek2,dsbessek1,dsbessek2
      integer i
      external wrate
      character*160 message

c-----------------------------------------------------------------------

      dsrdfunc=0.d0

c-----------------------------------------------------------------------
      u2=u*u                 ! u = sqrt((sqrt(s)-2m)/T)
      y=1.0d0+0.5d0*u2/x     ! y = sqrt(s)/(2*m)
      p=mco(1)*u*sqrt(0.5d0*(y+1.0d0)/x)
      if (.not.(p.ge.0.0d0)) then   ! je testing
        write(message,*) 'error in dsrdfunc: p=',p,'  mco(1)=',mco(1),
     &        '  u=',u,'  x=',x,'  y=',y
        call dswrite(0,1,message)
      endif
c------------------------------------------- invariant annihilation rate
      www=wrate(p)
      if (rderr.ne.0) return
c-------------------------------------------------------------- besselk1

      ek1=dsbessek1(2.0d0*x*y)  ! exp(2xy)*k_1(2xy) to avoid num. probl.

c-------------------------------------------------------------- besselk2

      denom=0.d0
      do i=1,nco
c...we use exp(x)*k_2(x_i) instead of k_2(x_i) to avoid numerical
c...problems for high x.
         ek2=dsbessek2(x*mco(i)/mco(1))*
     &     exp(x*(1.0d0-mco(i)/mco(1))) ! exp(x)*k2
         denom=denom+ek2*(mco(i)/mco(1))**2*(mdof(i)/mdof(1))
      enddo

c-------------------------------------------------- thermal distribution


      dsrdfunc=p**2*www*ek1/(mco(1)**4*mco(1)/x*denom**2)
c...take away exponentials introduced to avoid num. problems
      dsrdfunc=exp(2.0d0*x*(1.0d0-y))*dsrdfunc
c...change from peff -> u as integration variable
      dsrdfunc=dsrdfunc*u*mco(1)**2*y/(x*p)
c...multiply by m1^2 to make dimensionless and easier to integrate for
c...both low and high masses
      dsrdfunc=dsrdfunc*mco(1)**2

      return
      end

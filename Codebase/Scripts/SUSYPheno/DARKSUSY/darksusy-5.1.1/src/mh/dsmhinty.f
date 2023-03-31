***********************************************************************
*** dsmhinty integrates the ODE f'(x) = derivs(y,x) from xi to xf, 
*** using an adaptive stepsize control for Runge-Kutta 
*** (adapted from code provided by NUMERICAL RECIPES). 
*** further input: eps is the accuracy of the result
***
*** author: torsten bringmann (troms@physto.se), 2008-12-03
***********************************************************************

      subroutine dsmhinty(ystart,xi,xf,eps,hi,hmin,derivs,ier)
      implicit none

      integer MAXSTP,ier
      real*8 eps,hi,hmin,xi,xf,ystart
      external derivs
      parameter (MAXSTP=10000) !1000

      INTEGER nstp
      real*8 h,hdid,hnext,x,dydx,y,yscal

c... needed for the RK step control
      real*8 hh,htr,errmax,htemp,xnew,yerr,ytemp
      real*8 SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
c... The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.

c... needed for a single 5th order step
      real*8 ytempp,ak2,ak3,ak4,ak5,ak6,
     * A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     * B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,
     * DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     * B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     * B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     * B63=575./13824.,B64=44275./110592.,B65=253./4096.,
     * C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,
     * DC1=C1-2825./27648.,DC3=C3-18575./48384.,
     * DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)

      ier=0
      x=xi
      h=hi
      y=ystart
      do 100 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        yscal=abs(y) + abs(h*dydx)   ! to determine general scaling
        if((x+h-xf)*(x+h-xi).gt.0.) h=xf-x ! If stepsize can overshoot, decrease.

c... Take a fifth-order Runge-Kutta step with monitoring of local truncation error
        hh=h ! Set stepsize to the initial trial value.
   10   ytempp=y+B21*hh*dydx            ! First step
        call derivs(x+A2*hh,ytempp,ak2) ! Second step.
        ytempp=y+hh*(B31*dydx+B32*ak2)
        call derivs(x+A3*hh,ytempp,ak3) ! Third step.
        ytempp=y+hh*(B41*dydx+B42*ak2+B43*ak3)
        call derivs(x+A4*hh,ytempp,ak4) ! Fourth step.
        ytempp=y+hh*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs(x+A5*hh,ytempp,ak5) ! Fifth step.
        ytempp=y+hh*(B61*dydx+B62*ak2+B63*ak3+
     *              B64*ak4+B65*ak5)
        call derivs(x+A6*hh,ytempp,ak6) ! Sixth step.
        ytemp=y+hh*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
        yerr=hh*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
c... end single fifth-order RK-step

        errmax=0. ! Evaluate accuracy.
        errmax=max(errmax,abs(yerr/yscal))
        errmax=errmax/eps ! Scale relative to required tolerance.
        if(errmax.gt.1.) then ! Truncation error too large, reduce stepsize.
          htemp=SAFETY*hh*(errmax**PSHRNK)
          hh=sign(max(abs(htemp),0.1*abs(hh)),hh) ! No more than a factor of 10.
          xnew=x+hh
          if(xnew.eq.x) then
             write (*,*) 'stepsize underflow in rkqs'
             stop
          endif
          goto 10 ! For another try.
        else ! Step succeeded. Compute size of next step.
          if(errmax.gt.ERRCON) then
            hnext=SAFETY*hh*(errmax**PGROW)
          else ! No more than a factor of 5 increase.
            hnext=5.*hh
          endif
          hdid=hh
          x=x+hh
          y=ytemp
        endif

        if((x-xf)*(xf-xi).ge.0.) then ! end of interval reached 
          ystart=y
          xf=x
          return 
        endif
        if(abs(hnext).lt.hmin) then
          ier=2
          return
        endif
c        if(abs(hnext).lt.hmin) then
c           write (*,*) 'stepsize smaller than minimum in dsmhinty'
c           stop
c        endif
        h=hnext
 100    continue
       ier=1
c      write(*,*) 'too many steps in dsmhinty:'
c      write(*,*) 'integration only up to x = ',x
c      stop
      return 
      END

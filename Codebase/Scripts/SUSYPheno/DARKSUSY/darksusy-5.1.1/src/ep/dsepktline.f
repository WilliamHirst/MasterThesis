**********************************************************************
*** function dsepktline calculates the flux of e+ from the line
*** annihilation.
*** from kamionkowski & turner, prd 43(1991)1774.
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-02-10
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepktline(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 egev,q,tcont,yr2sec,hc2,
     &  gfunc,tau7,r,brepem
      parameter (yr2sec=3.15576d7)
      parameter (hc2=0.38937966d-27)  ! gev^2 cm^2
      parameter (tau7=1.0d0) ! yr
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsepktline: dshasetup must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      if (egev.le.hamwimp) then 
        brepem=habr(15) ! br to e+ e-   wtot -> sigmav pg990223
c        brepem=0.01d0  ! delete later

c...source strength
        q=(rhox/hamwimp)**2*hasv/(4.0d0*pi)*0.5d0 ! cm^-3 s^-1, JE Corr 03-01-21

c...containment time
        tcont=1d7*tau7*yr2sec*20.0/hamwimp  ! containment time, sec

c...greens function
        r=66.0d0/(tau7*20.0d0)
        gfunc=6.0d26*q/hamwimp**2*(egev/hamwimp)**(r-2.0d0)

c        gfunc=6.0d26*q/egev**2*(egev/hamwimp)**(r-2.0d0) ! just for testing

        dsepktline=gfunc*brepem

      else
        dsepktline=0.0d0
      endif

      return
      end


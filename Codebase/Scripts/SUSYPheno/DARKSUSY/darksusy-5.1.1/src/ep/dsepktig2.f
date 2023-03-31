**********************************************************************
*** function dsepktig2 is the positron spectrum times the greens
*** function from kamionkowski & turner, prd 43(1991)1774.
*** this routine is integrated by dsepktdiff to give the differential
*** positron flux at earth. the independent variable for this
*** routine is x=1/e**2 instead of e as in dsepktig.
*** units: gev^-2 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-02-10
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepktig2(x)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 x,eep,q,tcont,yr2sec,hc2,
     &  gfunc,tau7,r,phiep,dshaloyield
      integer istat
      parameter (yr2sec=3.15576d7)
      parameter (hc2=0.38937966d-27)  ! gev^2 cm^2
      parameter (tau7=1.0d0) ! yr
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsepktig2: dshasetup must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      eep=sqrt(1.0d0/x)
      if (eep.ge.ehere) then
        phiep=dshaloyield(eep,151,istat)

c...source strength
        q=(rhox/hamwimp)**2*hasv/(4.0d0*pi)*0.5d0  ! cm^-3 s^-1, JE Corr 03-01-21

c...containment time
        tcont=1d7*tau7*yr2sec*20.0/eep  ! containment time, sec

c...greens function
        r=66.0d0/(tau7*20.0d0)
        gfunc=6.0d26*q/eep**2*(ehere/eep)**(r-2.0d0)

c        gfunc=6.0d26*q/ehere**2*(ehere/eep)**(r-2.0d0) ! just for testing

        dsepktig2=gfunc*phiep
c...add differential
        dsepktig2=dsepktig2*(eep**3)/2.0d0
      else
        dsepktig2=0.0d0
      endif

      return
      end


**********************************************************************
*** function dsepmsig is the positron spectrum times the greens
*** function from moskalenko and strong prd 60, 063003 (1999)
*** this routine is integrated by dsepmsdiff to give the differential
*** positron flux at earth.
*** units: gev^-2 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se,
***   edward baltz eabaltz@alum.mit.edu
*** date: 2001 10/18
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepmsig(eep)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 eep,q,gfunc,phiep,dshaloyield,logehere,
     +     aa,bb,cc,ww,xx,yy
      integer istat
c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
         write(*,*) 'DS error in dsepmsig: dshasetup must be called',
     &        ' before any rate calculation'
         write(*,*) 'begins for every new model. program stopping.'
         stop
      endif

c...here includes reacceleration within 10 GeV
      if (eep+10.0d0.ge.ehere) then
         phiep=dshaloyield(eep,151,istat)
         logehere=log10(ehere)

c...source strength
         q=(rhox/hamwimp)**2*hasv*0.5d0  ! cm^-3 s^-1, JE Corr 03-01-21

c...green's function
c...first get fit parameters
         call dsepmstable(eep,aa,bb,cc,ww,xx,yy)
         gfunc=1.0d25*q/ehere**2
      
         if (ehere.le.eep) then
            gfunc=gfunc*10.0d0**((aa*logehere+bb)*logehere+cc)
         else
            gfunc=gfunc*10.0d0**((ww*logehere+xx)*logehere+yy)
         endif
      
         dsepmsig=gfunc*phiep
      else
         dsepmsig=0.0d0
      endif
      
      return
      end
      

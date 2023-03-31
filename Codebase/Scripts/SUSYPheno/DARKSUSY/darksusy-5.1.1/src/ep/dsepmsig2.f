**********************************************************************
*** function dsepmsig2 is the positron spectrum times the greens
*** function from moskalenko and strong, prd 60(1999)063003.
*** this routine is integrated by dsepmsdiff to give the differential
*** positron flux at earth. the independent variable for this
*** routine is x=1/e**2 instead of e as in dsepmsig.
*** units: gev^-2 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se,
***   e.a. baltz, eabaltz@alum.mit.edu
*** date: 01-10-18
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepmsig2(x)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 x,eep,q,gfunc,phiep,dshaloyield,logehere,
     +     aa,bb,cc,ww,xx,yy
      integer istat
c----------------------------------------------------------------------
      
      if (.not.dshasetupcalled) then
        write(*,*) 'error in dsepmsig2: dshasetup must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif
      
      eep=sqrt(1.0d0/x)
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
         
         dsepmsig2=gfunc*phiep
      
c...add differential
         dsepmsig2=dsepmsig2*(eep**3)/2.0d0
      else
         dsepmsig2=0.0d0
      endif

      return
      end


**********************************************************************
*** function dsepmsline calculates the flux of e+ from the line
*** annihilation.
*** from moskaleno & strong, prd 60(1999)063003.
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se,
***   e.a. baltz, eabaltz@alum.mit.edu
*** date: 01-10-19
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepmsline(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 egev,q,gfunc,brepem,logehere,
     +     aa,bb,cc,ww,xx,yy

c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
         write(*,*) 'DS error in dsepmsline: dshasetup must be called',
     &        ' before any rate calculation'
         write(*,*) 'begins for every new model. program stopping.'
         stop
      endif

      brepem=habr(15)  ! br to e+ e-   wtot -> sigmav pg990223
c        brepem=0.01d0  ! delete later
      logehere=log10(ehere)

c...source strength
      q=(rhox/hamwimp)**2*hasv*0.5d0  ! cm^-3 s^-1, JE Corr 03-01-21

c...greens function
c...first get fit parameters
      call dsepmstable(hamwimp,aa,bb,cc,ww,xx,yy)
      gfunc=1.0d25*q/ehere**2
      
      if (ehere.le.hamwimp) then
         gfunc=gfunc*10.0d0**((aa*logehere+bb)*logehere+cc)
      else
         gfunc=gfunc*10.0d0**((ww*logehere+xx)*logehere+yy)
      endif
      
      dsepmsline=gfunc*brepem
      
      return
      end

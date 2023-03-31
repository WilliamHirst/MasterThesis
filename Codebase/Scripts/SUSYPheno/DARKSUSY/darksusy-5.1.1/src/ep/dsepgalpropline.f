**********************************************************************
*** function dsepgalpropline calculates the flux of e+ from the line
*** annihilation, from GALPROP
*** units: gev^-1 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se,
***   e.a. baltz, eabaltz@alum.mit.edu
*** date: 4/27/2006
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsepgalpropline(egev)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsgalpropcom.h'

c----------------------------------------------------------------------

      real*8 egev,q,gfunc,brepem,logehere,fidxin,fidxout
      integer idxin,idxout

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
      q=(1.d0/hamwimp)**2*hasv*0.5d0 ! cm^-3 s^-1, JE Corr 03-01-21

c...greens function
c...first get fit parameters
      gfunc=q/ehere**2
      
      fidxin=(log10(hamwimp)+1.d0)*gpnumdecade
      idxin=int(fidxin)+1
      fidxin=fidxin+1.d0-idxin
      if (idxin.lt.1) idxin=1
      if (idxin.gt.gpnumin) idxin=gpnumin
      fidxout=(log10(ehere)+3.d0)*gpnumdecade
      idxout=int(fidxout)+1
      fidxout=fidxout+1.d0-idxout
      if (idxout.lt.1) idxout=1
      if (idxout.gt.gpnumout) idxout=gpnumout
      gfunc=gfunc*(
     +     (1.d0-fidxin)*(1.d0-fidxout)*epgpgf(idxin,idxout)+
     +     (1.d0-fidxin)*fidxout*epgpgf(idxin,idxout+1)+
     +     fidxin*(1.d0-fidxout)*epgpgf(idxin+1,idxout)+
     +     fidxin*fidxout*epgpgf(idxin+1,idxout+1))
      
      dsepgalpropline=gfunc*brepem
      
      return
      end

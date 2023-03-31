**********************************************************************
*** function dsgalpropig2 is the positron spectrum times the greens
*** function from galprop
*** this routine is integrated by dsepgalpropdiff to give the differential
*** positron flux at earth.  the independent variable for this
*** routine is x=1/e**2 instead of e as in dsepgalpropig.
*** units: gev^-2 cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se,
***   edward baltz eabaltz@alum.mit.edu
*** date: 2006 4/27
*** Modified: Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dsgalpropig2(x,pbar)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsgalpropcom.h'

c----------------------------------------------------------------------

      real*8 eep,q,gfunc,phiep,dshaloyield,logehere,x,
     +     fidxin,fidxout
      integer istat,idxin,idxout,pbar
c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
         write(*,*) 'DS error in dsepgalpropig2:',
     &        ' dshasetup must be called',
     &        ' before any rate calculation'
         write(*,*) 'begins for every new model. program stopping.'
         stop
      endif

      eep=sqrt(1.d0/x)
c...here includes reacceleration within 10 GeV
      if (eep+10.0d0.ge.ehere) then
         if (pbar.eq.0) then
            phiep=dshaloyield(eep,151,istat)
         else                   ! antiprotons
            phiep=dshaloyield(eep,154,istat)
         endif
         logehere=log10(ehere)

c...source strength
         q=(1.d0/hamwimp)**2*hasv*0.5d0 ! cm^-3 s^-1, JE Corr 03-01-21

c...green's function
         gfunc=q/ehere**2
      
         fidxin=(log10(eep)+1.d0)*gpnumdecade
         idxin=int(fidxin)+1
         fidxin=fidxin+1.d0-idxin
         if (idxin.lt.1) idxin=1
         if (idxin.gt.gpnumin) idxin=gpnumin
         fidxout=(log10(ehere)+3.d0)*gpnumdecade
         idxout=int(fidxout)+1
         fidxout=fidxout+1.d0-idxout
         if (idxout.lt.1) idxout=1
         if (idxout.gt.gpnumout) idxout=gpnumout
         if (pbar.eq.0) then
            gfunc=gfunc*(
     +           (1.d0-fidxin)*(1.d0-fidxout)*epgpgf(idxin,idxout)+
     +           (1.d0-fidxin)*fidxout*epgpgf(idxin,idxout+1)+
     +           fidxin*(1.d0-fidxout)*epgpgf(idxin+1,idxout)+
     +           fidxin*fidxout*epgpgf(idxin+1,idxout+1))
         else                   ! antiprotons
            gfunc=gfunc*(
     +           (1.d0-fidxin)*(1.d0-fidxout)*pbgpgf(idxin,idxout)+
     +           (1.d0-fidxin)*fidxout*pbgpgf(idxin,idxout+1)+
     +           fidxin*(1.d0-fidxout)*pbgpgf(idxin+1,idxout)+
     +           fidxin*fidxout*pbgpgf(idxin+1,idxout+1))
         endif
         dsgalpropig2=gfunc*phiep

c...  add differential
         dsgalpropig2=dsgalpropig2*(eep**3)/2.0d0
      else
         dsgalpropig2=0.0d0
      endif
      
      return
      end

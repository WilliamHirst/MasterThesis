**********************************************************************
*** function dsepipol interpolates in the table of
*** <sigma v> dn/de to speed up the
*** positron flux routines.
*** input: eep - positron energy
*** output: <sigma v> dn/de in units of cm^3 s^-1 gev^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: jun-02-98
**********************************************************************

      real*8 function dsepipol(eep)
      implicit none

      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dsepcom.h'

c----------------------------------------------------------------------

      integer i
      real*8 eep,ppl,lelow,lehigh

c----------------------------------------------------------------------

      if (eep.lt.emin) then
        dsepipol=0.0d0
        return
      endif

c...find entry in table
      i=int((log(eep)-lemin)*npt/dle)
      lelow=lemin+dle*dble(i)/dble(npt)
      lehigh=lemin+dle*dble(i+1)/dble(npt)
      ppl=(log(eep)-lelow)/(lehigh-lelow)

      dsepipol=(1.0d0-ppl)*dnde(min(i,npt)) + ppl*dnde(min(i+1,npt))

      return
      end







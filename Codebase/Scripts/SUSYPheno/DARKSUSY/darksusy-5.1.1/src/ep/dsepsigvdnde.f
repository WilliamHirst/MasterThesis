**********************************************************************
*** function dsepsigvdnde gives the differential spectrum of positrons
*** when they are created.
*** input: eep in gev
*** output: <sigma v> dn/de in units of cm^3 s^-1 gev^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: jun-01-98
*** modified: 99-07-02 paolo gondolo
**********************************************************************

      real*8 function dsepsigvdnde(eep)
      implicit none

      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      real*8 dshaloyield,eep
      integer istat

c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dsepdiff: dshasetup must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      hasmooth=2  ! 2 to smooth well
      dsepsigvdnde=hasv*dshaloyield(eep,151,istat)
      hasmooth=0

      return
      end






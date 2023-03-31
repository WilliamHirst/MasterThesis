**********************************************************************
*** subroutine dseptab tabulates <sigma v> dn/de to speed up the
*** positron flux routines. the interpolation is done with dsepipol.
*** input: emin - minimal energy for table
***   npts - number of points for the table
*** author: joakim edsjo, edsjo@physto.se
*** date: jun-02-98
**********************************************************************

      subroutine dseptab(em,npts)
      implicit none

      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dsepcom.h'
      include 'dsprep.h'

c----------------------------------------------------------------------

      integer i,npts
      real*8 dsepsigvdnde,e,em
      external dsepsigvdnde

c----------------------------------------------------------------------

      if (.not.dshasetupcalled) then
        write(*,*) 'DS error in dseptab: dshasetup must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
        stop
      endif

      if (npts.gt.npmax) then
        write(*,*) 'warning in dseptab: too many points requested.'
        write(*,*) 'increase npmac in dsepcom.h'
        write(*,*) 'will use npmax=',npmax,' instead.'
        npt=npmax
      else
        npt=npts
      endif

      emin=em
      lemin=log(emin)
      lemax=log(0.9999d0*hamwimp)
      dle=lemax-lemin
      do i=0,npt
        e=exp(lemin+dle*dble(i)/dble(npt))
        dnde(i)=dsepsigvdnde(e)
      enddo

      return
      end

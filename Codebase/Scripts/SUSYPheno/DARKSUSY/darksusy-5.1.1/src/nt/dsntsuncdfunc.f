***********************************************************************
*** dsntsuncdfunc returns the density of protons, neutrons or the total
*** density depending on the common block variable cdt. If
***   cdt='N': the total density is returned
***   cdt='p': the density in protons is returned
***   cdt='n': the density in neutrons is returned
*** the radius should be given in m and the density is returned in
*** g/cm^3.
*** This routine is used by dsntsuncdensint to calculate the column
*** density in the Sun.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-11-24
***********************************************************************

      real*8 function dsntsuncdfunc(r)
      implicit none

      include 'dssun.h'
      real*8 dsntsunmfrac,pfr,dsntsundens,r
      integer i,j

c...Check if data file is loaded
      call dsntsunread

      if (cdt.eq.'N') then
        dsntsuncdfunc=dsntsundens(r)
      elseif (cdt.eq.'p') then
        pfr=dsntsunmfrac(r,1)
        dsntsuncdfunc=dsntsundens(r)*(pfr + 0.5d0*(1.0d0-pfr))
      elseif (cdt.eq.'n') then
        pfr=dsntsunmfrac(r,1)
        dsntsuncdfunc=dsntsundens(r)*(0.5d0*(1.0d0-pfr))
      else
        write(*,*) 'ERROR in dsntsuncdfunc: wrong type: ',cdt
        stop
      endif

      return

      end

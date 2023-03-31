***********************************************************************
*** dsntsunne2x takes an input number density of electrons and
*** converts this to a fractional solar radius, x.
*** Input: n_e [cm^-3]
*** Output: x = r/r_sun [0,1]
*** See dsntsunread for information about which solar model is used.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2006-03-27
***********************************************************************

      real*8 function dsntsunne2x(ne)
      implicit none

      include 'dssun.h'
      include 'dsmpconst.h'

      real*8 ne,nepl,lne
      integer j

c...Check if data file is loaded
      call dsntsunread

      lne=log10(ne/n_avogadro) ! this is what is tabulated
      if (lne.ge.sdne(1)) then
        dsntsunne2x=0.0d0
        return
      endif

      if (lne.le.sdne(sdnne)) then
        dsntsunne2x=sdrne(sdnne)
        return
      endif

      call dshunt(sdne,sdnne,lne,j)
      if (j.lt.sdnne) goto 20
      
      dsntsunne2x=0.0d0
      return

 20   nepl=(lne-sdne(j))/(sdne(j+1)-sdne(j))

      dsntsunne2x=sdrne(j)*(1.0d0-nepl)+sdrne(j+1)*nepl

      return

      end

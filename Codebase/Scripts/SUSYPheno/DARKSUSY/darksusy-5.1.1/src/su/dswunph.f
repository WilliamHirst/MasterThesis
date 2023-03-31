      subroutine dswunph(unit,unphys)
c_______________________________________________________________________
c  write reasons for unphys<>0 to specified unit.
c  input:
c    unit - logical unit to write on (integer)
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer unit,unphys
      if (unit.le.0) return
      if (unphys.ne.0) write (unit,1000) 'unphys: '
      if (unphys.eq.1) write (unit,1002)
     &  'the lsp is not a neutralino. lsp=',pname(lsp)
      if (unphys.eq.-2) write (unit,1000)
     & 'incompatible m_a and tan(beta):',ma,tanbe
      if (unphys.eq.-3) write (unit,1000)
     & 'negative squared h2 mass at tan(beta):', tanbe
      if (unphys.eq.-4) write (unit,1000)
     & 'negative squared h2 mass at tan(beta):', tanbe
      if (unphys.eq.-1) write (unit,1001)
     & 'invalid or unavailable higloop:',higloop
      if (unphys.ne.0) write (unit,*)
1000  format (1x,a:1x,f8.3:1x,f8.3)
1001  format (1x,a,1x,i4)
1002  format (1x,a,a)
      end

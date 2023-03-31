************************************************************************
*** positron propagation routines.
*** author: e.a. baltz eabaltz@alum.mit.edu
*** date: 2001 10/18
************************************************************************


************************************************************************
      subroutine dsepmstable(eep,aa,bb,cc,ww,xx,yy)
************************************************************************
      implicit none

      include 'dsepcom.h'

      integer idxlo,idxhi
      real*8 eep,aa,bb,cc,ww,xx,yy,frac

      call dshunt(mselo,10,eep,idxlo)
      if (idxlo.eq.0) then
         aa=msatab(1)
         bb=msbtab(1)
         cc=msctab(1)
      else if (idxlo.eq.10) then
         aa=msatab(10)
         bb=msbtab(10)
         cc=msctab(10)
      else
         frac=(eep-mselo(idxlo))/(mselo(idxlo+1)-mselo(idxlo))
         aa=msatab(idxlo)+frac*(msatab(idxlo+1)-msatab(idxlo))
         bb=msbtab(idxlo)+frac*(msbtab(idxlo+1)-msbtab(idxlo))
         cc=msctab(idxlo)+frac*(msctab(idxlo+1)-msctab(idxlo))
      endif

      call dshunt(msehi,3,eep,idxhi)
      if (idxhi.eq.0) then
         ww=mswtab(1)
         xx=msxtab(1)
         yy=msytab(1)
      else if (idxhi.eq.3) then
         ww=mswtab(3)
         xx=msxtab(3)
         yy=msytab(3)
      else
         frac=(eep-msehi(idxhi))/(msehi(idxhi+1)-msehi(idxhi))
         ww=mswtab(idxhi)+frac*(mswtab(idxhi+1)-mswtab(idxhi))
         xx=msxtab(idxhi)+frac*(msxtab(idxhi+1)-msxtab(idxhi))
         yy=msytab(idxhi)+frac*(msytab(idxhi+1)-msytab(idxhi))
      endif
      
      end
      

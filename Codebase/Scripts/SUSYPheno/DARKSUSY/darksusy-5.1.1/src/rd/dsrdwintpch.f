      subroutine dsrdwintpch(p,wspline,wlin)
c_______________________________________________________________________
c  check of interpolation of tabulated invariant rate.
c  input:
c    p - initial cm momentum (real)
c  common:
c    'dsrdcom.h' - included common blocks
c  called by dsrdtab.
c  author: joakim edsjo 96-04-10
c  based on wintp.f by p. gondolo but wlin is also calculated
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      integer inc,nrm,ihi,ilo
      real*8 p,y,a,b,h,wspline,wlin
c-----------------------------------------------------------------------
      if (nlo.le.0.or.nlo.gt.nr) then
         nlo=0
         nhi=nr+1
         goto 3
      endif
      inc=1
      if (p.ge.pp(indx(nlo))) then
    1   nhi=nlo+inc
        if (nhi.gt.nr) then
          nhi=nr+1
        else if (p.ge.pp(indx(nhi))) then
          nlo=nhi
          inc=inc+inc
          goto 1
        endif
      else
        nhi=nlo
    2   nlo=nhi-inc
        if (nlo.lt.1) then
          nlo=0
        else if (p.lt.pp(indx(nlo))) then
          nhi=nlo
          inc=inc+inc
          goto 2
        endif
      endif
    3 if (nhi-nlo.ne.1) then
        nrm=(nhi+nlo)/2
        if (p.gt.pp(indx(nrm))) then
          nlo=nrm
        else
          nhi=nrm
        endif
        goto 3
      endif
      if (nlo.eq.0) then
         write (rdluerr,*) 'warning! dsrdwintpch: index undershoot'
         write (rdluerr,*) p,pp(indx(1)),pp(indx(2)),pp(indx(nr))
      else if (nlo.eq.nr) then
         write (rdluerr,*) 'warning! dsrdwintpch: index overshoot'
         write (rdluerr,*) '  p=',p,' max p in table =',pp(indx(nr))
         write (rdluerr,*) '  mco(1) = ',mco(1)
      endif
      ihi=indx(nlo+1)
      ilo=indx(nlo)
      h=pp(ihi)-pp(ilo)
      if (h.eq.0) then
         write (rdluerr,*) 'dsrdwintpch: coincident p''s in wx(p) table'
         write (*,*) 'ilo = ',ilo,'  pp(ilo) = ',pp(ilo)
         write (*,*) 'ihi = ',ihi,'  pp(ihi) = ',pp(ihi)
         stop
      endif
      a=(pp(ihi)-p)/h
      b=(p-pp(ilo))/h
      y=a*yy(ilo)+b*yy(ihi)+
     & ((a**3-a)*yy2(nlo)+(b**3-b)*yy2(nhi))*(h**2)/6.0d0

c      write(*,*)
c      write(*,*) 'p = ',p
c      write(*,*) 'pp(ihi) = ',pp(ihi)
c      write(*,*) 'pp(ilo) = ',pp(ilo)
c      write(*,*) 'h = ',h
c      write(*,*) 'a = ',a
c      write(*,*) 'b = ',b
c      write(*,*) 'yy(ilo) = ',yy(ilo)
c      write(*,*) 'yy(ihi) = ',yy(ihi)
c      write(*,*) 'yy2(nlo) = ',yy2(nlo)
c      write(*,*) 'yy2(nhi) = ',yy2(nhi)
c      write(*,*) 'y = ',y
c      write(*,*) 'exp(y) = ',exp(y)

c...avoid overflow, je dec 11, 1998
      if (y.gt.600.0d0) then
        y=600.0d0
        rdwar=ibset(rdwar,4)
      endif
      wspline=exp(y)

      wlin=exp(a*yy(ilo)+b*yy(ihi))

      end







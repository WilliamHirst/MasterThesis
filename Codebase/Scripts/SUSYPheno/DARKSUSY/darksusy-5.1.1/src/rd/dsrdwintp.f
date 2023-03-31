      function dsrdwintp(p)
c_______________________________________________________________________
c  interpolation of tabulated invariant rate.
c  input:
c    p - initial cm momentum (real)
c  common:
c    'dsrdcom.h' - included common blocks
c  called by dsrdfunc.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 dsrdwintp,ptop
      integer inc,nrm,ihi,ilo
      real*8 p,y,a,b,h,wspline,ptmp
c-----------------------------------------------------------------------

      ptmp=min(p,0.99999d0*pmax)  ! Extrapolate above pmax
      ptmp=max(ptmp,pp(indx(1)))  ! JE corr 05-10-31


c...find index in pp array
      if (nlo.le.0.or.nlo.gt.nr) then
         nlo=0
         nhi=nr+1
         goto 3
      endif
      inc=1
      if (ptmp.ge.pp(indx(nlo))) then
    1   nhi=nlo+inc
        if (nhi.gt.nr) then
          nhi=nr+1
        else if (ptmp.ge.pp(indx(nhi))) then
          nlo=nhi
          inc=inc+inc
          goto 1
        endif
      else
        nhi=nlo
    2   nlo=nhi-inc
        if (nlo.lt.1) then
          nlo=0
        else if (ptmp.lt.pp(indx(nlo))) then
          nhi=nlo
          inc=inc+inc
          goto 2
        endif
      endif
    3 if (nhi-nlo.ne.1) then
        nrm=(nhi+nlo)/2
        if (ptmp.ge.pp(indx(nrm))) then   ! JE corr. 050617: .gt. -> .ge.
          nlo=nrm
        else
          nhi=nrm
        endif
        goto 3
      endif
      if (nlo.eq.0) then
         write (rdluerr,*) 'warning! dsrdwintp: index undershoot'
         write (rdluerr,*) ptmp,pp(indx(1)),pp(indx(2)),pp(indx(nr))
         write (rdluerr,*) '  model: ',rdtag
         dsrdwintp=0.0d0
      else if (nlo.eq.nr) then
         write (rdluerr,*) 'warning! dsrdwintp: index overshoot'
         write (rdluerr,*) '  ptmp = ',ptmp
         write (rdluerr,*) '  model: ',rdtag
         write (rdluerr,*) '  nlo=',nlo,'  nr=',nr
         write (rdluerr,*) '  pp(indx(nr)) = ',pp(indx(nr))
         dsrdwintp=0.0d0
      endif
      ihi=indx(nlo+1)
      ilo=indx(nlo)
      h=pp(ihi)-pp(ilo)
      if (h.eq.0) then
         write (rdluerr,*) 'dsrdwintp: coincident p''s in wx(p) table'
         write (*,*) 'ilo = ',ilo,'  pp(ilo) = ',pp(ilo)
         write (*,*) 'ihi = ',ihi,'  pp(ihi) = ',pp(ihi)
         write (rdluerr,*) '  model: ',rdtag
         stop
      endif

c...interpolate
      a=(pp(ihi)-ptmp)/h
      b=(ptmp-pp(ilo))/h
      y=a*yy(ilo)+b*yy(ihi)+
     & ((a**3-a)*yy2(nlo)+(b**3-b)*yy2(nhi))*(h**2)/6.0d0
      wspline=exp(y)
      dsrdwintp=wspline

c      write (*,*) 'PG-DEBUG dsrdwintp: wspline,y,a,b,yy(ilo),yy(ihi)',
c     &     wspline,y,a,b,yy(ilo),yy(ihi)

 1000 format(a,f6.2,a)

      end

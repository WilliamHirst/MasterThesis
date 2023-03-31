      subroutine dsrdspline
c_______________________________________________________________________
c  set up 2nd derivatives for cubic dsrdspline interpolation.
c  common:
c    'dsrdcom.h' - included common blocks
c  called by dsrdtab.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c  modified by joakim edsjo, edsjo@physto.se, to split spline
c  at thresholds.
c  modified: april 30, 1998.
c  modified: april 24, 2011 by pat scott (patscott@physics.mcgill.ca)
c   to allow code to run with array bounds-checking enabled 
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      integer i,k,ii,ip1,im1,istart,iend
      real*8 p,qn,sig,un,u(nrmax)
      logical dsrdthtest
c-----------------------------------------------------------------------
c...set-up ispl array where the subregions for the spline are stored
      yy2(1)=-0.5d0
c     yy2(1)=0.0d0 if zero 2nd derivative
      u(1)=0.0d0
      istart=1

 10   yy2(istart)=0.0d0
      u(istart)=0.0d0

      do 20 i=istart+1,nr-1
        ii=indx(i)
        ip1=indx(i+1)
        im1=indx(i-1)
        sig=(pp(ii)-pp(im1))/(pp(ip1)-pp(im1))
        p=sig*yy2(i-1)+2.0d0
        yy2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((yy(ip1)-yy(ii))/
     &    (pp(ip1)-pp(ii))-(yy(ii)-yy(im1))/
     &   (pp(ii)-pp(im1)))/(pp(ip1)-pp(im1))-sig*u(i-1))/p
c...check if at end or at threshold
c...PS bug fix 20110424; split evaluation of conditional to prevent indexing pp(0)
        if (i+1.eq.nr) then
          iend=i+1
          goto 30
        endif
        if (dsrdthtest(i+1)) then
          iend=i+1
          goto 30
        endif
   20 continue

c...jump out of loop
 30   qn=0.0d0
      un=0.0d0
      yy2(iend)=(un-qn*u(iend-1))/(qn*yy2(iend-1)+1.0d0)
      do 11 k=iend-1,istart,-1
        yy2(k)=yy2(k)*yy2(k+1)+u(k)
   11 continue

      if (iend.ne.nr) then
        istart=iend+1 ! skip one point
        goto 10
      endif

      end











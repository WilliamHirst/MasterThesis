      real*8 function dsrddpmin(p,dpmin)
c_______________________________________________________________________
c  routine to determine if there is a narrow resonance present which
c  jusifies changing dpmin to some fraction of lambda
c  author: joakim edsjo, edsjo@physto.se
c  date: april 30, 1998
c  modified: april 30, 1998.
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 p,dpmin,lamfr,dplamfr,tmp,p1,p12,p2,p22
      integer i
      parameter(lamfr=8.0d0,dplamfr=0.1d0)
c...if we want to sample well until the resonance has decreased by a
c...factor of x, we should go to lamfr=0.5*sqrt(x)
c-----------------------------------------------------------------------
      tmp=dpmin
      do i=1,nres
        p12=(rgev(i)-lamfr*rwid(i))**2/4.0d0-mco(1)**2
        if (p12.gt.0.0d0) then
          p1=sqrt(p12)
        else
          p1=0.0d0
        endif
        p22=(rgev(i)+lamfr*rwid(i))**2/4.0d0-mco(1)**2
        if (p22.gt.0.0d0) then
          p2=sqrt(p22)
        else
          p2=0.0d0
        endif

        if (p.ge.p1.and.p.le.p2.and.dplamfr*rwid(i).lt.tmp)
     &    tmp=dplamfr*rwid(i)
      enddo

      dsrddpmin=tmp

      return
      end











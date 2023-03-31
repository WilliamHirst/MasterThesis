      subroutine dsantucn(p,ind1,ind2)
c_______________________________________________________________________
c  routine to check when t- or u-resonances occur.
c  called by dwdcosopt.
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-09-17
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      integer h,kh(3),f,k,kp1,kp2,fp,g,ind1,ind2,dsantures
      real*8 p,ptmp

c-----------------------------------------------------------------------

      kp1=kcha(ind1-4)
      kp2=kn(ind2)

      ptmp=p
c      p=max(p,0.0001d0)

      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3

c============================================================ c n -> h+ h
      do h=1,2
        if (dsantures(kp1,kp2,khc,kh(h),p).eq.1) then
          tur(ind1,ind2,h)=.true.
          tur(ind2,ind1,h)=.true.
        endif
c h+ h
      enddo

c============================================================ c n -> h+ a
      if (dsantures(kp1,kp2,khc,kh3,p).eq.1) then
        tur(ind1,ind2,3)=.true.
        tur(ind2,ind1,3)=.true.
      endif
c h+ a

c============================================================ c n -> w+ h
      do h=1,2
        if (dsantures(kp1,kp2,kw,kh(h),p).eq.1) then
          tur(ind1,ind2,3+h)=.true.
          tur(ind2,ind1,3+h)=.true.
        endif
c w+ h
      enddo

c============================================================ c n -> w+ a
      if (dsantures(kp1,kp2,kw,kh3,p).eq.1) then
        tur(ind1,ind2,6)=.true.
        tur(ind2,ind1,6)=.true.
      endif
c h+ a

c============================================================ c n -> z h+
      if (dsantures(kp1,kp2,kz,khc,p).eq.1) then
        tur(ind1,ind2,7)=.true.
        tur(ind2,ind1,7)=.true.
      endif
c z h+

c======================================================== c n -> gamma h+
      if (dsantures(kp1,kp2,kgamma,khc,p).eq.1) then
        tur(ind1,ind2,8)=.true.
        tur(ind2,ind1,8)=.true.
      endif
c gamma h+

c============================================================ c n -> z w+
      if (dsantures(kp1,kp2,kz,kw,p).eq.1) then
        tur(ind1,ind2,9)=.true.
        tur(ind2,ind1,9)=.true.
      endif
c z w+

c======================================================== c n -> gamma w+
      if (dsantures(kp1,kp2,kgamma,kw,p).eq.1) then
        tur(ind1,ind2,10)=.true.
        tur(ind2,ind1,10)=.true.
      endif
c gamma w+

c========================================================= c n -> f' fbar
      do k=1,6
        f = 2*k-1
        do g=1,3
          if (f.eq.1.or.f.eq.3.or.f.eq.5) then          ! nue, numu, nutau
            fp = g*2
          else !if (f.eq.7.or.f.eq.9.or.f.eq.11) then   ! u, c, t
            fp = 2*g+6
          endif
          if (dsantures(kp1,kp2,f,fp,p).eq.1) then
            tur(ind1,ind2,7+k*3+g)=.true.
            tur(ind2,ind1,7+k*3+g)=.true.
          endif
        enddo
      enddo
c f fbar

c=======================================================================

      end












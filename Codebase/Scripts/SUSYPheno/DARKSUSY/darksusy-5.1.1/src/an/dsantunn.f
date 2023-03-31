      subroutine dsantunn(p,ind1,ind2)
c_______________________________________________________________________
c  routine to check if t- or u-resonances occur for neutralino-neutralino
c  annihilation
c  called by dwdcosopt.
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-09-17
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      integer h,i,ii,kh(3),k,kp1,kp2,ind1,ind2,dsantures
      real*8 p,ptmp
c-----------------------------------------------------------------------

      kp1=kn(ind1)
      kp2=kn(ind2)

      ptmp=p
c      p=max(p,0.000001d0)
c      p=max(p,0.0001d0)

      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3

c=================================================== n n -> hi hii / a a
      do k=1,4
        if (k.eq.1) then
          i=1
          ii=1
        elseif (k.eq.2) then
          i=1
          ii=2
        elseif (k.eq.3) then
          i=2
          ii=2
        else
          i=3
          ii=3
        endif
        if (dsantures(kp1,kp2,kh(i),kh(ii),p).eq.1) then
          tur(ind1,ind2,min(i+ii-1,4))=.true.
          tur(ind2,ind1,min(i+ii-1,4))=.true.
        endif

c hi hii / a a
      enddo

c============================================================ n n -> h a
      do h=1,2
        if (dsantures(kp1,kp2,kh(h),kh3,p).eq.1) then
          tur(ind1,ind2,4+h)=.true.
          tur(ind2,ind1,4+h)=.true.
        endif
c h a
      enddo

c========================================================== n n -> h- h+

c h- h+

c============================================================ n n -> h z
      do h=1,2
        if (dsantures(kp1,kp2,kz,kh(h),p).eq.1) then
          tur(ind1,ind2,7+h)=.true.
          tur(ind2,ind1,7+h)=.true.
        endif
c h z
      enddo

c============================================================ n n -> a z
      if (dsantures(kp1,kp2,kz,kh3,p).eq.1) then
        tur(ind1,ind2,10)=.true.
        tur(ind2,ind1,10)=.true.
      endif
c a z


c============================================================ n n -> h w
c     h+ w- calculated and the result is multiplied by 2
      if (dsantures(kp1,kp2,kw,khc,p).eq.1) then
        tur(ind1,ind2,11)=.true.
        tur(ind2,ind1,11)=.true.
      endif
c h w


c============================================================ n n -> z z

c z z


c============================================================ n n -> w w

c w w


c========================================================= n n -> f fbar

c f fbar

      return

      end












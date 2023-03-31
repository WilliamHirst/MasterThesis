      subroutine dsantucc(p,ind1,ind2)
c_______________________________________________________________________
c  routine to check for t- or u-channel resonances.
c  called by dwdcosopt.
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-09-17
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      integer h,i,ii,kh(3),f,k,kp1,kp2,g,fp,ind1,ind2,dsantures
      real*8 p,ptmp

c-----------------------------------------------------------------------
      kp1=kcha(ind1-4)
      kp2=kcha(ind2-4)

      ptmp=p
c      p=max(p,0.0001d0)

      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3

c********************************************** opposite charge channels
c=================================================== c c -> hi hii / a a
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

c============================================================ c c -> h a
      do h=1,2
        if (dsantures(kp1,kp2,kh(h),kh3,p).eq.1) then
          tur(ind1,ind2,4+h)=.true.
          tur(ind2,ind1,4+h)=.true.
        endif
c h a
      enddo

c========================================================== c c -> h+ h-

c h- h+

c============================================================ c c -> h z
      do h=1,2
        if (dsantures(kp1,kp2,kz,kh(h),p).eq.1) then
          tur(ind1,ind2,7+h)=.true.
          tur(ind2,ind1,7+h)=.true.
        endif
c h z
      enddo

c============================================================ c c -> a z
      if (dsantures(kp1,kp2,kz,kh3,p).eq.1) then
        tur(ind1,ind2,10)=.true.
        tur(ind2,ind1,10)=.true.
      endif
c a z


c============================================================ c c -> h w
c     h+ w- calculated and the result is multiplied by 2
      if (dsantures(kp1,kp2,khc,kw,p).eq.1) then
        tur(ind1,ind2,11)=.true.
        tur(ind2,ind1,11)=.true.
      endif
c h w


c============================================================ c c -> z z

c z z

c============================================================ c c -> w w

c w w

c==================================================== c c -> gamma gamma

c gamma gamma

c======================================================= c c -> z0 gamma
      if (dsantures(kp1,kp2,kz,kgamma,p).eq.1) then
        tur(ind1,ind2,15)=.true.
        tur(ind2,ind1,15)=.true.
      endif
c z0 gamma

c======================================================== c c -> f fbar'
      do f=1,12
        do g=1,3
          if (f.eq.1.or.f.eq.3.or.f.eq.5) then          ! nue, numu, nutau
            fp=g*2-1
          else if (f.eq.2.or.f.eq.4.or.f.eq.6) then     ! e, mu, tau
            fp=g*2
          else if (f.eq.7.or.f.eq.9.or.f.eq.11) then    ! u, c, t
            fp=g*2+5
          else !if (f.eq.8.or.f.eq.10.or.f.eq.12) then  ! d, s, b
            fp=g*2+6
          endif
          if (dsantures(kp1,kp2,f,fp,p).eq.1) then
            tur(ind1,ind2,12+f*3+g)=.true.
            tur(ind2,ind1,12+f*3+g)=.true.
          endif
        enddo
      enddo
c f fbar'

c=======================================================================

c************************************************* equal charge channels

c========================================================== c c -> h+ h+

c h+ h+

c========================================================== c c -> h+ w+
      if (dsantures(kp1,kp2,khc,kw,p).eq.1) then
        tur(ind1,ind2,53)=.true.
        tur(ind2,ind1,53)=.true.
      endif
c h+ w+

c========================================================== c c -> w+ w+

c w+ w+

      end












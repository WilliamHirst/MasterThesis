      function dsandwdcoscc(p,costheta,kp1,kp2)
c_______________________________________________________________________
c  annihilation differential invariant rate between particle kp1
c  and kp2 where kp1 and kp2 are charginos
c  input:
c    p - initial cm momentum (real)
c    costheta - cosine of c.m. annihilation angle
c    kp1 - particle code, particle 1
c    kp2 - particle code, particle 2
c  Output:
c  BeginTex
c   \begin{displaymath}
c     \frac{dW_{ij}}{d\cos\theta}
c   \end{displaymath}
c   where
c   \begin{displaymath}
c      W_{ij} = 4 p_{ij} \sqrt{s} \sigma_{ij} 
c      = 4 \sigma_{ij} \sqrt{(p_i \cdot p_j)^2 - m_i^2 m_j^2} 
c      = 4 E_{i} E_{j} \sigma_{ij} v_{ij}
c   \end{displaymath}
c  EndTex
c  The returned dW/dcos(theta) is unitless
c  uses dsanclearaa,dsansumaa
c  called by dsandwdcos.
c  author: joakim edsjo (edsjo@physto.se)
c  date: 96-08-06
c  modified: 01-09-12
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandiacom.h'
      include 'dsandwcom.h'
      real*8 dsandwdcoscc
      integer h,i,ii,j,c,kh(3),f,ksf,k,kp1,kp2,g,fp,ind1,ind2
      real*8 w,p,costheta,sumf,dsansumaa
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c-----------------------------------------------------------------------
      if (p.lt.0.0d0) then
        dsandwdcoscc=0.0d0
        return
      endif

      ind1=kp1-kcha(1)+5
      ind2=kp2-kcha(1)+5

c      p=max(p,0.0001d0)   ! akinvar stable now, ok with p=0

      w=0.0d0
      do i=1,54
        prtial(i)=0.0d0
      enddo
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
        call dsanclearaa
        do c=1,2
          call dsantfffss(p,costheta,kp1,kp2,kcha(c),kh(i),kh(ii))
          call dsantfffssex(p,costheta,kp1,kp2,kcha(c),kh(i),kh(ii))
        enddo
        call dsansffsss(p,costheta,kp1,kp2,kh1,kh(i),kh(ii))
        call dsansffsss(p,costheta,kp1,kp2,kh2,kh(i),kh(ii))
        sumf = dsansumaa()
        if (i.eq.ii) sumf=0.5d0*sumf ! identical in final state
        prtial(min(i+ii-1,4))=(kk/sqrt(s))*sumf/(32.d0*pi)
c hi hii / a a
      enddo

c============================================================ c c -> h a
      do h=1,2
        call dsanclearaa
        do c=1,2
          call dsantfffss(p,costheta,kp1,kp2,kcha(c),kh(h),kh3)
          call dsantfffssex(p,costheta,kp1,kp2,kcha(c),kh(h),kh3)
        enddo
        call dsansffsss(p,costheta,kp1,kp2,kh3,kh(h),kh3)
        call dsansffvss(p,costheta,kp1,kp2,kz,kh(h),kh3)
        prtial(4+h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h a
      enddo

c========================================================== c c -> h+ h-
      call dsanclearaa
      do j=1,4
        call dsantfffss(p,costheta,kp1,kp2,kn(j),khc,khc)
      enddo
      call dsansffsss(p,costheta,kp1,kp2,kh1,khc,khc)
      call dsansffsss(p,costheta,kp1,kp2,kh2,khc,khc)
      call dsansffvss(p,costheta,kp1,kp2,kz,khc,khc)
      call dsansffvss(p,costheta,kp1,kp2,kgamma,khc,khc)
      prtial(7)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h- h+

c============================================================ c c -> h z
      do h=1,2
        call dsanclearaa
        do c=1,2
          call dsantfffvs(p,costheta,kp1,kp2,kcha(c),kz,kh(h))
          call dsantfffvsex(p,costheta,kp1,kp2,kcha(c),kz,kh(h))
        enddo
        call dsansffsvs(p,costheta,kp1,kp2,kh3,kz,kh(h))
        call dsansffvvs(p,costheta,kp1,kp2,kz,kz,kh(h))
        prtial(7+h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h z
      enddo

c============================================================ c c -> a z
      call dsanclearaa
      do c=1,2
        call dsantfffvs(p,costheta,kp1,kp2,kcha(c),kz,kh3)
        call dsantfffvsex(p,costheta,kp1,kp2,kcha(c),kz,kh3)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,kh1,kz,kh3)
      call dsansffsvs(p,costheta,kp1,kp2,kh2,kz,kh3)
      prtial(10)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c a z


c============================================================ c c -> h w
c     h+ w- calculated and the result is multiplied by 2
      call dsanclearaa
      do j=1,4
        call dsantfffsv(p,costheta,kp1,kp2,kn(j),khc,kw)
      enddo
      call dsansffssv(p,costheta,kp1,kp2,kh1,khc,kw)
      call dsansffssv(p,costheta,kp1,kp2,kh2,khc,kw)
      call dsansffssv(p,costheta,kp1,kp2,kh3,khc,kw)
      prtial(11)=2.0d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h w


c============================================================ c c -> z z
      call dsanclearaa
      do c=1,2
        call dsantfffvv(p,costheta,kp1,kp2,kcha(c),kz,kz)
        call dsantfffvvex(p,costheta,kp1,kp2,kcha(c),kz,kz)
      enddo
      call dsansffsvv(p,costheta,kp1,kp2,kh1,kz,kz)
      call dsansffsvv(p,costheta,kp1,kp2,kh2,kz,kz)
      prtial(12)=0.5d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c z z

c============================================================ c c -> w w
      call dsanclearaa
      do j=1,4
        call dsantfffvvex(p,costheta,kp1,kp2,kn(j),kw,kw)
      enddo
      call dsansffsvv(p,costheta,kp1,kp2,kh1,kw,kw)
      call dsansffsvv(p,costheta,kp1,kp2,kh2,kw,kw)
      call dsansffvvv(p,costheta,kp1,kp2,kz,kw,kw)
      call dsansffvvv(p,costheta,kp1,kp2,kgamma,kw,kw)
      prtial(13)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c w w

c==================================================== c c -> gamma gamma
      if (kp1.eq.kp2) then
        call dsanclearaa
        call dsantfffvv(p,costheta,kp1,kp2,kp1,kgamma,kgamma)
        call dsantfffvvex(p,costheta,kp1,kp2,kp1,kgamma,kgamma)
        prtial(14)=0.5d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
      endif
c gamma gamma

c======================================================= c c -> z0 gamma
      call dsanclearaa
      call dsantfffvv(p,costheta,kp1,kp2,kp2,kz,kgamma)
      call dsantfffvvex(p,costheta,kp1,kp2,kp1,kz,kgamma)
      prtial(15)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
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
          call dsanclearaa
c... sfermion exchange
          do i=1,6
            if (f.eq.1.or.f.eq.3.or.f.eq.5) then          ! nue, numu, nutau
              ksf=ksl(i)
              call dsantffsff(p,costheta,kp1,kp2,ksf,f,fp)
            else if (f.eq.2.or.f.eq.4.or.f.eq.6) then     ! e, mu, tau
              if (i.le.3) then
                ksf=ksnu(i)
                call dsanuffsff(p,costheta,kp1,kp2,ksf,f,fp)
              endif
            else if (f.eq.7.or.f.eq.9.or.f.eq.11) then    ! u, c, t
              ksf=ksqd(i)
              call dsantffsff(p,costheta,kp1,kp2,ksf,f,fp)
            else !if (f.eq.8.or.f.eq.10.or.f.eq.12) then  ! d, s, b
              ksf=ksqu(i)
              call dsanuffsff(p,costheta,kp1,kp2,ksf,f,fp)
            endif
          enddo
c... call subroutines to calculate amplitudes
          if (f.eq.1.or.f.eq.3.or.f.eq.5) then
            call dsansffvff(p,costheta,kp1,kp2,kz,f,fp)
          else
            call dsansffsff(p,costheta,kp1,kp2,kh1,f,fp)
            call dsansffsff(p,costheta,kp1,kp2,kh2,f,fp)
            call dsansffsff(p,costheta,kp1,kp2,kh3,f,fp)
            call dsansffvff(p,costheta,kp1,kp2,kz,f,fp)
            call dsansffvff(p,costheta,kp1,kp2,kgamma,f,fp)
          endif
          prtial(12+f*3+g)=ncolor(f)*(kk/sqrt(s))*dsansumaa()/
     &      (32.d0*pi)
        enddo
      enddo
c f fbar'

c=======================================================================

c************************************************* equal charge channels

c========================================================== c c -> h+ h+
      call dsanclearaa
      do j=1,4
        call dsantfffssin(p,costheta,kp1,kp2,kn(j),khc,khc)
        call dsanufffssin(p,costheta,kp1,kp2,kn(j),khc,khc)
      enddo
      prtial(52)=0.5d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h+ h+

c========================================================== c c -> h+ w+
      call dsanclearaa
      do j=1,4
        call dsantfffsvin(p,costheta,kp1,kp2,kn(j),khc,kw)
        call dsanufffsvin(p,costheta,kp1,kp2,kn(j),khc,kw)
      enddo
      prtial(53)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h+ w+

c========================================================== c c -> w+ w+
      call dsanclearaa
      do j=1,4
        call dsantfffvvin(p,costheta,kp1,kp2,kn(j),kw,kw)
        call dsanufffvvin(p,costheta,kp1,kp2,kn(j),kw,kw)
      enddo
      prtial(54)=0.5d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c w+ w+

c------------------------------------------------------- sum up partials
c... note that the invariant rate should be
c... w = 0.5 [sum w(\chi^+ \chi^+ ->...) + sum w(\chi^+ \chi^- ->...)]
c...   = 0.5 [sum w(\chi^- \chi^- ->...) + sum w(\chi^- \chi^+ ->...)]
c... since the chargino is treated as one particle with 4 degrees of
c... freedom.

      w=0.d0
      w=w+prtial(1)             ! h1 h1
      w=w+prtial(2)             ! h1 h2
      w=w+prtial(3)             ! h2 h2
      w=w+prtial(4)             ! a a
      w=w+prtial(5)             ! h1 a
      w=w+prtial(6)             ! h2 a
      w=w+prtial(7)             ! h+ h-
      w=w+prtial(8)             ! h1 z
      w=w+prtial(9)             ! h2 z
      w=w+prtial(10)            ! a z
      w=w+prtial(11)            ! h+ w- & h- w+
      w=w+prtial(12)            ! z z
      w=w+prtial(13)            ! w+ w-
      w=w+prtial(14)            ! gamma gamma
      w=w+prtial(15)            ! z0 gamma
      w=w+prtial(16)            ! nu_e nu_e-bar
      w=w+prtial(17)            ! nu_e nu_mu-bar
      w=w+prtial(18)            ! nu_e nu_tau-bar
      w=w+prtial(19)            ! e e-bar
      w=w+prtial(20)            ! e mu-bar
      w=w+prtial(21)            ! e tau-bar
      w=w+prtial(22)            ! nu_mu nu_e-bar
      w=w+prtial(23)            ! nu_mu nu_mu-bar
      w=w+prtial(24)            ! nu_mu nu_tau-bar
      w=w+prtial(25)            ! mu e-bar
      w=w+prtial(26)            ! mu mu-bar
      w=w+prtial(27)            ! mu tau-bar
      w=w+prtial(28)            ! nu_tau nu_e-bar
      w=w+prtial(29)            ! nu_tau nu_mu-bar
      w=w+prtial(30)            ! nu_tau nu_tau-bar
      w=w+prtial(31)            ! tau e-bar
      w=w+prtial(32)            ! tau mu-bar
      w=w+prtial(33)            ! tau tau-bar
      w=w+prtial(34)            ! u u-bar
      w=w+prtial(35)            ! u c-bar
      w=w+prtial(36)            ! u t-bar
      w=w+prtial(37)            ! d d-bar
      w=w+prtial(38)            ! d s-bar
      w=w+prtial(39)            ! d b-bar
      w=w+prtial(40)            ! c u-bar
      w=w+prtial(41)            ! c c-bar
      w=w+prtial(42)            ! c t-bar
      w=w+prtial(43)            ! s d-bar
      w=w+prtial(44)            ! s s-bar
      w=w+prtial(45)            ! s b-bar
      w=w+prtial(46)            ! t u-bar
      w=w+prtial(47)            ! t c-bar
      w=w+prtial(48)            ! t t-bar
      w=w+prtial(49)            ! b d-bar
      w=w+prtial(50)            ! b s-bar
      w=w+prtial(51)            ! b b-bar
      w=w+prtial(52)            ! h+ h+
      w=w+prtial(53)            ! h+ w+
      w=w+prtial(54)            ! w+ w+

      dsandwdcoscc = 0.5d0*w        ! see above for the 0.5d0

c      write(*,*) 'cc: w = ',w

      end












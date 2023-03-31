      function dsandwdcoscn(p,costheta,kp1,kp2)
c_______________________________________________________________________
c  annihilation differential invariant rate between particle kp1
c  and kp2 where kp1 is a chargino and kp2 is a neutralino.
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
      real*8 dsandwdcoscn
      integer h,i,j,c,kh(3),f,ksf1,ksf2,k,kp1,kp2,fp,g,ind1,ind2
      real*8 w,p,costheta,dsansumaa
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c-----------------------------------------------------------------------
      if (p.lt.0.0d0) then
        dsandwdcoscn=0.0d0
        return
      endif

      ind1=kp1-kcha(1)+5
      ind2=kp2-kn(1)+1

c      p=max(p,0.001d0)     ! akinvar stable now, ok with p=0

      w=0.0d0
      do i=1,54
        prtial(i)=0.0d0
      enddo
      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3

c============================================================ c n -> h+ h
      do h=1,2
        call dsanclearaa
        do c=1,2
          call dsantfffssex(p,costheta,kp1,kp2,kcha(c),khc,kh(h))
        enddo
        do j=1,4
          call dsantfffss(p,costheta,kp1,kp2,kn(j),khc,kh(h))
        enddo
        call dsansffsss(p,costheta,kp1,kp2,khc,khc,kh(h))
        call dsansffvss(p,costheta,kp1,kp2,kw,khc,kh(h))
        prtial(h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h+ h
      enddo

c============================================================ c n -> h+ a
      call dsanclearaa
      do c=1,2
        call dsantfffssex(p,costheta,kp1,kp2,kcha(c),khc,kh3)
      enddo
      do j=1,4
        call dsantfffss(p,costheta,kp1,kp2,kn(j),khc,kh3)
      enddo
      call dsansffvss(p,costheta,kp1,kp2,kw,khc,kh3)
      prtial(3)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h+ a

c============================================================ c n -> w+ h
      do h=1,2
        call dsanclearaa
        do c=1,2
          call dsantfffvsex(p,costheta,kp1,kp2,kcha(c),kw,kh(h))
        enddo
        do j=1,4
          call dsantfffvs(p,costheta,kp1,kp2,kn(j),kw,kh(h))
        enddo
        call dsansffsvs(p,costheta,kp1,kp2,khc,kw,kh(h))
        call dsansffvvs(p,costheta,kp1,kp2,kw,kw,kh(h))
        prtial(3+h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c w+ h
      enddo

c============================================================ c n -> w+ a
      call dsanclearaa
      do c=1,2
        call dsantfffvsex(p,costheta,kp1,kp2,kcha(c),kw,kh3)
      enddo
      do j=1,4
        call dsantfffvs(p,costheta,kp1,kp2,kn(j),kw,kh3)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,khc,kw,kh3)
      prtial(6)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h+ a

c============================================================ c n -> z h+
      call dsanclearaa
      do j=1,4
        call dsantfffvsex(p,costheta,kp1,kp2,kn(j),kz,khc)
      enddo
      do c=1,2
        call dsantfffvs(p,costheta,kp1,kp2,kcha(c),kz,khc)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,khc,kz,khc)
      prtial(7)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c z h+

c======================================================== c n -> gamma h+
      call dsanclearaa
      do c=1,2
        call dsantfffvs(p,costheta,kp1,kp2,kcha(c),kgamma,khc)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,khc,kgamma,khc)
      prtial(8)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c gamma h+

c============================================================ c n -> z w+
      call dsanclearaa
      do j=1,4
        call dsantfffvvex(p,costheta,kp1,kp2,kn(j),kz,kw)
      enddo
      do c=1,2
        call dsantfffvv(p,costheta,kp1,kp2,kcha(c),kz,kw)
      enddo
      call dsansffvvv(p,costheta,kp1,kp2,kw,kz,kw)
      prtial(9)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c z w+

c======================================================== c n -> gamma w+
      call dsanclearaa
      do c=1,2
        call dsantfffvv(p,costheta,kp1,kp2,kcha(c),kgamma,kw)
      enddo
      call dsansffvvv(p,costheta,kp1,kp2,kw,kgamma,kw)
      prtial(10)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
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
          call dsanclearaa
c... sfermion exchange
          do i=1,6
            if (f.eq.1.or.f.eq.3.or.f.eq.5) then        ! nue, numu, nutau
              ksf1=ksl(i)
              ksf2=ksnu(i)
            else !if (f.eq.7.or.f.eq.9.or.f.eq.11) then ! u, c, t
              ksf1=ksqd(i)
              ksf2=ksqu(i)
            endif
            call dsantffsff(p,costheta,kp1,kp2,ksf1,f,fp)
            if (ksf2.gt.0) call dsanuffsff(p,costheta,kp1,kp2,ksf2,f,fp)
          enddo
          call dsansffsff(p,costheta,kp1,kp2,khc,f,fp)
          call dsansffvff(p,costheta,kp1,kp2,kw,f,fp)
          prtial(7+k*3+g)=ncolor(f)*(kk/sqrt(s))*dsansumaa()/
     &      (32.d0*pi)
        enddo
      enddo
c f fbar

c=======================================================================
c------------------------------------------------------- sum up partials
c... note that the summed invariant rate should be
c... w = sum w(\chi^0 + \chi^+ -> ... ) = sum(\chi^0 + \chi^0 -> ...)
c... i.e., either \chi^0 \chi^+ or \chi^0 \chi^- and not the sum of the
c... two since the chargino is treated as one particle with 4 degrees of
c... freedom.

      w=0.d0
      w=w+prtial(1)             ! h+ h1
      w=w+prtial(2)             ! h+ h2
      w=w+prtial(3)             ! h+ a
      w=w+prtial(4)             ! w+ h1
      w=w+prtial(5)             ! w+ h2
      w=w+prtial(6)             ! w+ h3
      w=w+prtial(7)             ! z h+
      w=w+prtial(8)             ! gamma h+
      w=w+prtial(9)             ! z w+
      w=w+prtial(10)            ! gamma w+
      w=w+prtial(11)            ! nu_e e-bar
      w=w+prtial(12)            ! nu_e mu-bar
      w=w+prtial(13)            ! nu_e tau-bar
      w=w+prtial(14)            ! nu_mu e-bar
      w=w+prtial(15)            ! nu_mu mu-bar
      w=w+prtial(16)            ! nu_mu tau-bar
      w=w+prtial(17)            ! nu_tau e-bar
      w=w+prtial(18)            ! nu_tau mu-bar
      w=w+prtial(19)            ! nu_tau tau-bar
      w=w+prtial(20)            ! u d-bar
      w=w+prtial(21)            ! u s-bar
      w=w+prtial(22)            ! u b-bar
      w=w+prtial(23)            ! c d-bar
      w=w+prtial(24)            ! c s-bar
      w=w+prtial(25)            ! c b-bar
      w=w+prtial(26)            ! t d-bar
      w=w+prtial(27)            ! t s-bar
      w=w+prtial(28)            ! t b-bar

      dsandwdcoscn = w

c      write(*,*)
c      write(*,*) 'dsandwdcoscn called with:'
c      write(*,*) '  kp1=',kp1,' kp2=',kp2,' p=',p,' cth=',costheta
c      write(*,*) 'cn:  w = ',w
c      do i=7,10
c        write(*,*) '  prtial(',i,') = ',prtial(i)
c      enddo

      end

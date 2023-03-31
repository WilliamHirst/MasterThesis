      function dsandwdcosnn(p,costheta,kp1,kp2)
c_______________________________________________________________________
c  Annihilation differential invariant rate between particle kp1
c  and kp2 where kp1 and kp2 are neutralinos
c  Input:
c    p - initial cm momentum (real)
c    costheta - cosine of c.m. annihilation angle
c    kp1 - particle code for particle 1
c    kp2 - particle code for particle 2
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
c  The returned dW/dcos(theta) is unitless.
c     
c  uses dsanclearaa,dsansumaa
c  called by dsandwdcos.
c  note: the 32pi in the partial cross sections is 8pi g_1^2, g_1=2
c  author: joakim edsjo (edsjo@physto.se)
c  date: 96-02-21
c  modified: 01-09-12
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsprep.h'
      include 'dsidtag.h'
      include 'dsandiacom.h'
      include 'dsandwcom.h'
      real*8 dsandwdcosnn
      integer h,i,ii,j,c,kh(3),f,ksf,k,kp1,kp2,calc,ind1,ind2
      real*8 w,p,costheta,
     & sumf,dsansumaa,imres,reres,gluons,gammas,zgam
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      save gluons,gammas,zgam

c-----------------------------------------------------------------------

      calc=0 ! just for printing when gluons calculated
c...initialize akinvar variables
      call dsankinvar(0.0d0,0.0d0,0,0,0,0,0)

      ind1=kp1-kn(1)+1
      ind2=kp2-kn(1)+1

      if (p.lt.0.0d0) then
        dsandwdcosnn=0.0d0
        return
      endif

c      p=max(p,0.001d0)  ! akinvar stable now, ok with p=0

      w=0.0d0
      do i=1,54
        prtial(i)=0.0d0
      enddo
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
        call dsanclearaa
        do j=1,4
          call dsantfffss(p,costheta,kp1,kp2,kn(j),kh(i),kh(ii))
          call dsanufffss(p,costheta,kp1,kp2,kn(j),kh(i),kh(ii))
        enddo
        call dsansffsss(p,costheta,kp1,kp2,kh1,kh(i),kh(ii))
        call dsansffsss(p,costheta,kp1,kp2,kh2,kh(i),kh(ii))
        sumf = dsansumaa()
        if (i.eq.ii) sumf=0.5d0*sumf ! identical in final state
        prtial(min(i+ii-1,4))=(kk/sqrt(s))*sumf/(32.d0*pi)
c hi hii / a a
      enddo

c============================================================ n n -> h a
      do h=1,2
        call dsanclearaa
        do j=1,4
          call dsantfffss(p,costheta,kp1,kp2,kn(j),kh(h),kh3)
          call dsanufffss(p,costheta,kp1,kp2,kn(j),kh(h),kh3)
        enddo
        call dsansffsss(p,costheta,kp1,kp2,kh3,kh(h),kh3)
        call dsansffvss(p,costheta,kp1,kp2,kz,kh(h),kh3)
        prtial(4+h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h a
      enddo

c========================================================== n n -> h- h+
      call dsanclearaa
      do c=1,2
        call dsantfffss(p,costheta,kp1,kp2,kcha(c),khc,khc)
        call dsanufffss(p,costheta,kp1,kp2,kcha(c),khc,khc)
      enddo
      call dsansffsss(p,costheta,kp1,kp2,kh1,khc,khc)
      call dsansffsss(p,costheta,kp1,kp2,kh2,khc,khc)
      call dsansffvss(p,costheta,kp1,kp2,kz,khc,khc)
      prtial(7)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h- h+

c============================================================ n n -> h z
      do h=1,2
        call dsanclearaa
        do j=1,4
          call dsantfffvs(p,costheta,kp1,kp2,kn(j),kz,kh(h))
          call dsanufffvs(p,costheta,kp1,kp2,kn(j),kz,kh(h))
        enddo
        call dsansffsvs(p,costheta,kp1,kp2,kh3,kz,kh(h))
        call dsansffvvs(p,costheta,kp1,kp2,kz,kz,kh(h))
        prtial(7+h)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h z
      enddo

c============================================================ n n -> a z
      call dsanclearaa
      do j=1,4
        call dsantfffvs(p,costheta,kp1,kp2,kn(j),kz,kh3)
        call dsanufffvs(p,costheta,kp1,kp2,kn(j),kz,kh3)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,kh1,kz,kh3)
      call dsansffsvs(p,costheta,kp1,kp2,kh2,kz,kh3)
      prtial(10)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c a z


c============================================================ n n -> h w
c     h+ w- calculated and the result is multiplied by 2
      call dsanclearaa
      do c=1,2
        call dsantfffvs(p,costheta,kp1,kp2,kcha(c),kw,khc)
        call dsanufffvs(p,costheta,kp1,kp2,kcha(c),kw,khc)
      enddo
      call dsansffsvs(p,costheta,kp1,kp2,kh1,kw,khc)
      call dsansffsvs(p,costheta,kp1,kp2,kh2,kw,khc)
      call dsansffsvs(p,costheta,kp1,kp2,kh3,kw,khc)
      prtial(11)=2.0d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c h w


c============================================================ n n -> z z
      call dsanclearaa
      do j=1,4
        call dsantfffvv(p,costheta,kp1,kp2,kn(j),kz,kz)
        call dsanufffvv(p,costheta,kp1,kp2,kn(j),kz,kz)
      enddo
      call dsansffsvv(p,costheta,kp1,kp2,kh1,kz,kz)
      call dsansffsvv(p,costheta,kp1,kp2,kh2,kz,kz)
      prtial(12)=0.5d0*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c z z


c============================================================ n n -> w w
      call dsanclearaa
      do c=1,2
        call dsantfffvv(p,costheta,kp1,kp2,kcha(c),kw,kw)
        call dsanufffvv(p,costheta,kp1,kp2,kcha(c),kw,kw)
      enddo
      call dsansffsvv(p,costheta,kp1,kp2,kh1,kw,kw)
      call dsansffsvv(p,costheta,kp1,kp2,kh2,kw,kw)
      call dsansffvvv(p,costheta,kp1,kp2,kz,kw,kw)
      prtial(13)=(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c w w


c========================================================= n n -> f fbar
      do f=1,12
        call dsanclearaa
c... sfermion exchange
        do i=1,6
          if (f.eq.1.or.f.eq.3.or.f.eq.5) then          ! nue, numu, nutau
            ksf=ksnu(i)
          else if (f.eq.2.or.f.eq.4.or.f.eq.6) then     ! e, mu, tau
            ksf=ksl(i)
          else if (f.eq.7.or.f.eq.9.or.f.eq.11) then    ! u, c, t
            ksf=ksqu(i)
          else !if (f.eq.8.or.f.eq.10.or.f.eq.12) then  ! d, s, b
            ksf=ksqd(i)
          endif
          if (ksf.gt.0) then
            call dsantffsff(p,costheta,kp1,kp2,ksf,f,f)
            call dsanuffsff(p,costheta,kp1,kp2,ksf,f,f)
          endif
        enddo
c... call subroutines to calculate amplitudes
        call dsansffsff(p,costheta,kp1,kp2,kh1,f,f)
        call dsansffsff(p,costheta,kp1,kp2,kh2,f,f)
        call dsansffsff(p,costheta,kp1,kp2,kh3,f,f)
        call dsansffvff(p,costheta,kp1,kp2,kz,f,f)

        prtial(13+f)=ncolor(f)*(kk/sqrt(s))*dsansumaa()/(32.d0*pi)
c f fbar
      enddo

      if (kp1.eq.kn(1).and.kp2.eq.kn(1)) then

c===================================================== n n -> glue glue
c... routines by l. bergstrom 97-02-17
c... modified: 97-04-22 by l. bergstrom
c... modified: 97-05-14 by j. edsjo to use new w-convention
        if (incglue) then
          if (newmodelandwdcosnn) then
            call dsanglglim(imres)
            call dsanglglre(reres)
c right normalization:
c        prtial(26)=(imres**2+reres**2)*(s/4.)*alph3**2/8./pi**3*
c     &     mass(kn(kln))**2
            prtial(26)=(imres**2+reres**2)*alph3**2/8./pi**3*
     &       mass(kn(kln))**2 ! only p=0 computed; guess 1/s behaviour
                                     ! maybe it should be constant?
c division by 2 to compensate integration over dtheta
            prtial(26)=prtial(26)/2.
c multiplication by 4 mx^2 to get inv. rate
            prtial(26)=prtial(26)*4.0d0*mass(kn(kln))**2
            gluons=prtial(26)
          else
            prtial(26)=gluons
          endif
        else
          prtial(26)=0.0d0
        endif

c===================================================== n n -> qq glue
        prtial(27)=0.d0 ! no qq glue yet


c=======================================================n n -> 2 gamma
c... modified: 97-06-05 by j. edsjo to use new w-convention
        if (incgaga) then
          if (newmodelandwdcosnn) then
            imres=0.0d0
            reres=0.0d0
            call dsanggim(imres)
            call dsanggre(reres)
            prtial(28)=(imres**2+reres**2)*(mass(kn(kln))**2)*
     &        alphem**2/16./pi**3
c division by 2 to compensate integration over dtheta
            prtial(28) =prtial(28)/2.
c multiplication by 4 mx^2 to get inv. rate
            prtial(28)=prtial(28)*4.0d0*mass(kn(1))**2
c correct for w propagator, je addition 97-06-12
c          prtial(28)=prtial(28)*
c     &      (mass(kcha(2))**2-mass(kn(1))**2-mass(kw)**2)**2/
c     &      ((mass(kcha(2))**2-mass(kn(1))**2-mass(kw)**2)**2
c     &      +(width(kw)*mass(kw))**2)
            gammas=prtial(28)
          else
            prtial(28)=gammas
          endif
        else
          prtial(28)=0.0d0
        endif
c   gamma gamma

c=======================================================n n -> z gamma
c... modified: 97-06-05 by j. edsjo to use new w-convention
        if (incgaz) then
          if (newmodelandwdcosnn) then
            imres=0.0d0
            reres=0.0d0
            if (mass(kn(kln)).gt.mass(kz)/2.d0) then
            call dsanzg(reres,imres)
            endif
            prtial(29)=(imres**2+reres**2)*(mass(kn(kln))**2-
     &                mass(kz)**2/4.d0)/mass(kn(kln))**4*
     &        alphem/32./pi**4  !  k**3 replaced by k  - lbe 970627
c division by 2 to compensate integration over dtheta
            prtial(29) =prtial(29)/2.
c multiplication by 4 mx^2 to get inv. rate
            prtial(29)=prtial(29)*4.0d0*mass(kn(1))**2
c correct for w propagator, je addition 97-06-12
c          prtial(29)=prtial(29)*
c     &      (mass(kcha(2))**2-mass(kn(1))**2-mass(kw)**2)**2/
c     &      ((mass(kcha(2))**2-mass(kn(1))**2-mass(kw)**2)**2
c     &      +(width(kw)*mass(kw))**2)
            zgam=prtial(29)
          else
            prtial(29)=zgam
          endif
        else
          prtial(29)=0.0d0
        endif
c   z gamma

        if (newmodelandwdcosnn) then
          newmodelandwdcosnn=.false.
          calc=1
        endif
c=======================================================================

      else
        prtial(26)=0.0d0
        prtial(27)=0.0d0
        prtial(28)=0.0d0
        prtial(29)=0.0d0
      endif
c=======================================================================
c------------------------------------------------------- sum up partials
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
      w=w+prtial(14)            ! nue nuebar
      w=w+prtial(15)            ! e+ e-
      w=w+prtial(16)            ! numu numubar
      w=w+prtial(17)            ! mu+ mu-
      w=w+prtial(18)            ! nutau nutaubar
      w=w+prtial(19)            ! tau+ tau-
      w=w+prtial(20)            ! u ubar
      w=w+prtial(21)            ! d dbar
      w=w+prtial(22)            ! c cbar
      w=w+prtial(23)            ! s sbar
      w=w+prtial(24)            ! t tbar
      w=w+prtial(25)            ! b bbar
c... loop diagrams
      w=w+prtial(26)            ! glue glue
      w=w+prtial(27)            ! qq glue
      w=w+prtial(28)            ! gamma gamma
      w=w+prtial(29)            ! z gamma

      dsandwdcosnn = w

      end












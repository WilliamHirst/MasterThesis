      subroutine dsvertx1
c_______________________________________________________________________
c  some couplings used in neutralino-neutralino, neutralino-chargino
c  and chargino-chargino annihilation.
c  common:
c    'dsmssm.h' - file with susy common blocks
c  called by susyin.
c  needs neusct, chasct, sfesct, higsct.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1999
c  history:
c    951110 complex vertex constants
c    970213 joakim edsjo
c    990724 paolo gondolo trilinear higgs and goldstone couplings
c    020710 Joakim Edsjo, chargino-(up-squark)-(down-quark) sign-change
c    020903 Mia Schelke, Higgs-squark-squark, A-terms sign-change
c=======================================================================
c
c  vertices included:
c     z-w-w, z-h-h, w-h-a, w-h-h, h-w-w, h-z-z, z-a-h, h-h-h, h-a-a,
c     h-h-h, a-f-f, h-f-f, z-f-f, a-n-n, h-n-n, z-n-n, w-n-c, h-n-c,
c     squark-gluino-quark, sf-n-f, h-c-c, a-c-c, squark-squark-higgs,
c     w-f-f', h-f-f', gamma-w-w, gamma-h-h, z-c-c, gamma-c-c
c     gamma-f-f
c     gld-h-h, gld-gld-h, gld-n-c
c     Z-f~-f~, gamma-f~-f~, gluon-f~-f~
c
      implicit none
      include 'dsmssm.h'
      real*8 g2o2cw,aux,la3t,ca,sa,s2a,c2a,ca2,sa2,ca3,sa3,beta,cb,sb,
     &     c2b,s2b,cb2,sb2,cb3,sb3,cbpa,sbpa,cbma,sbma,e_charge,sw2
      real*8 ghaa(2),gh2h(2),ghhh(2,2,2),ghff,gaff,af,vf,gzah,ghzz,
     & gwhh,gwha,gzhh,ghww,gzww,gwff,glhff,grhff,ggww,gghh
      real*8 rmass(12) ! temporary array for lepton and quark running masses
                       ! currently extracted from yukawa 
                       ! one day we will have rmass(k,q)
      integer i,j,k,l,c,g,d
      complex*16 qij,sij,gjll,gjlr,gjrl,gjrr,
     & glscf,grscf,glhnc,grhnc,glwnc,grwnc,gzij,ghij,gaij,
     & qcd,scd,qdc,sdc,ocdl,ocdr
      real*8 rgl(50*50*50*2),rgr(50*50*50*2)
      equivalence (gl,rgl)
      equivalence (gr,rgr)

c...... temporary fix for running masses (PG 2002-06-16)
      aux=dsqrt(2.d0)*mass(kw)/g2weak
      do i=1,3
         rmass(knu(i)) = aux*yukawa(knu(i))*sinbe
         rmass(kl(i)) = aux*yukawa(kl(i))*cosbe
         rmass(kqu(i)) = aux*yukawa(kqu(i))*sinbe
         rmass(kqd(i)) = aux*yukawa(kqd(i))*cosbe
      enddo
c...... end of temporary fix for running masses
      g2o2cw=g2weak/(2.0d0*costhw)
      e_charge=g2weak*sinthw
      sw2=sinthw**2
      ca=cos(alpha)
      sa=sin(alpha)
      ca2=ca*ca
      sa2=sa*sa
      s2a=2.d0*sa*ca
      c2a=ca2-sa2
      ca3=ca2*ca
      sa3=sa2*sa
      beta=atan(tanbe)
      cb=cos(beta)
      sb=sin(beta)
      cb2=cb*cb
      sb2=sb*sb
      s2b=2.d0*sb*cb
      c2b=cb2-sb2
      cb3=cb2*cb
      sb3=sb2*sb
      cbpa=cos(beta+alpha) ! corr june 94
      sbpa=sin(beta+alpha) ! corr june 94
      cbma=cos(beta-alpha) ! corr june 94
      sbma=sin(beta-alpha) ! corr june 94

c--------------------------------------------------- clear vertex matrix
c...Note: This is needed on some compilers. Don't comment out
      do i=1,50*50*50*2
         rgl(i)=0.d0
         rgr(i)=0.d0
      enddo
c---------------------------------------------------------- gauge bosons

      gzww=g2weak*costhw
      gl(kz,kw,kw)=dcmplx(gzww,0.d0)
      gl(kw,kz,kw)=-gl(kz,kw,kw)          ! corr. je 96-02-27
      gl(kw,kw,kz)=gl(kz,kw,kw)
      gr(kz,kw,kw)=gl(kz,kw,kw)
      gr(kw,kz,kw)=gl(kw,kz,kw)           ! corr. je 96-02-27
      gr(kw,kw,kz)=gl(kw,kw,kz)

      ggww=g2weak*sinthw                  ! added by je 96-02-21
      gl(kgamma,kw,kw)=dcmplx(ggww,0.d0)
      gl(kw,kgamma,kw)=-gl(kgamma,kw,kw)  ! corr. je 96-02-27
      gl(kw,kw,kgamma)=gl(kgamma,kw,kw)
      gr(kgamma,kw,kw)=gl(kgamma,kw,kw)
      gr(kw,kgamma,kw)=gl(kw,kgamma,kw)   ! corr. je 96-02-27
      gr(kw,kw,kgamma)=gl(kw,kw,kgamma)

c   gluon-gluon-gluon. addition by pg 02-03-21
c     Feynman rule = -i g3stro f^{a_1 a_2 a_3} *
c                    [ (k_1-k_2)_{mu_3} g_{mu_1 mu_2} + cyclic perm's ]
c   made explicitly complex by je 02-07-11
      gl(kgluon,kgluon,kgluon)=dcmplx(0.0d0,g3stro)
      gr(kgluon,kgluon,kgluon)=gl(kgluon,kgluon,kgluon)

c------------------------------------------------ gauge and higgs bosons

      gzhh=-g2o2cw*(costhw**2-sinthw**2)
      gl(kz,khc,khc)=dcmplx(gzhh,0.d0)
      gl(khc,kz,khc)=gl(kz,khc,khc)
      gl(khc,khc,kz)=gl(kz,khc,khc)
      gr(kz,khc,khc)=gl(kz,khc,khc)
      gr(khc,kz,khc)=gl(khc,kz,khc)
      gr(khc,khc,kz)=gl(khc,khc,kz)

      gghh=-g2weak*sinthw                 ! added by je 96-02-21
      gl(kgamma,khc,khc)=dcmplx(gghh,0.d0)
      gl(khc,kgamma,khc)=gl(kgamma,khc,khc)
      gl(khc,khc,kgamma)=gl(kgamma,khc,khc)
      gr(kgamma,khc,khc)=gl(kgamma,khc,khc)
      gr(khc,kgamma,khc)=gl(khc,kgamma,khc)
      gr(khc,khc,kgamma)=gl(khc,khc,kgamma)

      gwha=-0.5d0*g2weak
      gl(kw,khc,kh3)=dcmplx(0.d0,gwha)     ! corr. je 960212 w h+ a
      gl(khc,kw,kh3)=gl(kw,khc,kh3)
      gl(khc,kh3,kw)=gl(kw,khc,kh3)
      gl(kw,kh3,khc)=dcmplx(0.d0,-gwha)    ! corr. je 960212 w a h+
      gl(kh3,kw,khc)=gl(kw,kh3,khc)
      gl(kh3,khc,kw)=gl(kw,kh3,khc)
      gr(kw,khc,kh3)=gl(kw,khc,kh3)
      gr(khc,kw,kh3)=gl(khc,kw,kh3)
      gr(khc,kh3,kw)=gl(khc,kh3,kw)
      gr(kw,kh3,khc)=gl(kw,kh3,khc)
      gr(kh3,kw,khc)=gl(kh3,kw,khc)
      gr(kh3,khc,kw)=gl(kh3,khc,kw)

      gwhh=+0.5d0*g2weak*sbma
      gl(kw,khc,kh1)=dcmplx(gwhh,0.d0)     ! w h+ h1
      gl(khc,kw,kh1)=gl(kw,khc,kh1)
      gl(khc,kh1,kw)=gl(kw,khc,kh1)
      gl(kw,kh1,khc)=dcmplx(gwhh,0.d0)     ! w h1 h+
      gl(kh1,kw,khc)=gl(kw,kh1,khc)
      gl(kh1,khc,kw)=gl(kw,kh1,khc)
      gr(kw,khc,kh1)=gl(kw,khc,kh1)
      gr(khc,kw,kh1)=gl(khc,kw,kh1)
      gr(khc,kh1,kw)=gl(khc,kh1,kw)
      gr(kw,kh1,khc)=gl(kw,kh1,khc)
      gr(kh1,kw,khc)=gl(kh1,kw,khc)
      gr(kh1,khc,kw)=gl(kh1,khc,kw)

      gwhh=-0.5d0*g2weak*cbma
      gl(kw,khc,kh2)=dcmplx(gwhh,0.d0)     ! w h+ h2
      gl(khc,kw,kh2)=gl(kw,khc,kh2)
      gl(khc,kh2,kw)=gl(kw,khc,kh2)
      gl(kw,kh2,khc)=dcmplx(gwhh,0.d0)     ! w h2 h+
      gl(kh2,kw,khc)=gl(kw,kh2,khc)
      gl(kh2,khc,kw)=gl(kw,kh2,khc)
      gr(kw,khc,kh2)=gl(kw,khc,kh2)
      gr(khc,kw,kh2)=gl(khc,kw,kh2)
      gr(khc,kh2,kw)=gl(khc,kh2,kw)
      gr(kw,kh2,khc)=gl(kw,kh2,khc)
      gr(kh2,kw,khc)=gl(kh2,kw,khc)
      gr(kh2,khc,kw)=gl(kh2,khc,kw)

      ghww=+g2weak*cbma
      gl(kh1,kw,kw)=dcmplx(ghww,0.d0)
      gl(kw,kh1,kw)=gl(kh1,kw,kw)
      gl(kw,kw,kh1)=gl(kh1,kw,kw)
      gr(kh1,kw,kw)=gl(kh1,kw,kw)
      gr(kw,kh1,kw)=gl(kw,kh1,kw)
      gr(kw,kw,kh1)=gl(kw,kw,kh1)

      ghww=+g2weak*sbma
      gl(kh2,kw,kw)=dcmplx(ghww,0.d0)
      gl(kw,kh2,kw)=gl(kh2,kw,kw)
      gl(kw,kw,kh2)=gl(kh2,kw,kw)
      gr(kh2,kw,kw)=gl(kh2,kw,kw)
      gr(kw,kh2,kw)=gl(kw,kh2,kw)
      gr(kw,kw,kh2)=gl(kw,kw,kh2)

      ghzz=+2.0d0*g2o2cw*cbma/costhw     ! corr. je 960208
      gl(kh1,kz,kz)=dcmplx(ghzz,0.d0)
      gl(kz,kh1,kz)=gl(kh1,kz,kz)
      gl(kz,kz,kh1)=gl(kh1,kz,kz)
      gr(kh1,kz,kz)=gl(kh1,kz,kz)
      gr(kz,kh1,kz)=gl(kz,kh1,kz)
      gr(kz,kz,kh1)=gl(kz,kz,kh1)

      ghzz=+2.0d0*g2o2cw*sbma/costhw     ! corr. je 960208
      gl(kh2,kz,kz)=dcmplx(ghzz,0.d0)
      gl(kz,kh2,kz)=gl(kh2,kz,kz)
      gl(kz,kz,kh2)=gl(kh2,kz,kz)
      gr(kh2,kz,kz)=gl(kh2,kz,kz)
      gr(kz,kh2,kz)=gl(kz,kh2,kz)
      gr(kz,kz,kh2)=gl(kz,kz,kh2)

      gzah=+g2o2cw*sbma
      gl(kz,kh1,kh3)=dcmplx(0.d0,-gzah)   ! corr. je 960208
      gl(kh1,kz,kh3)=gl(kz,kh1,kh3)
      gl(kh1,kh3,kz)=gl(kz,kh1,kh3)
      gl(kz,kh3,kh1)=dcmplx(0.d0,gzah)    ! corr. je 960208
      gl(kh3,kz,kh1)=gl(kz,kh3,kh1)
      gl(kh3,kh1,kz)=gl(kz,kh3,kh1)
      gr(kz,kh1,kh3)=gl(kz,kh1,kh3)
      gr(kh1,kz,kh3)=gl(kh1,kz,kh3)
      gr(kh1,kh3,kz)=gl(kh1,kh3,kz)
      gr(kz,kh3,kh1)=gl(kz,kh3,kh1)
      gr(kh3,kz,kh1)=gl(kh3,kz,kh1)
      gr(kh3,kh1,kz)=gl(kh3,kh1,kz)

      gzah=-g2o2cw*cbma
      gl(kz,kh2,kh3)=dcmplx(0.d0,-gzah)   ! corr. je 960208
      gl(kh2,kz,kh3)=gl(kz,kh2,kh3)
      gl(kh2,kh3,kz)=gl(kz,kh2,kh3)
      gl(kz,kh3,kh2)=dcmplx(0.d0,gzah)    ! corr. je 960208
      gl(kh3,kz,kh2)=gl(kz,kh3,kh2)
      gl(kh3,kh2,kz)=gl(kz,kh3,kh2)
      gr(kz,kh2,kh3)=gl(kz,kh2,kh3)
      gr(kh2,kz,kh3)=gl(kh2,kz,kh3)
      gr(kh2,kh3,kz)=gl(kh2,kh3,kz)
      gr(kz,kh3,kh2)=gl(kz,kh3,kh2)
      gr(kh3,kz,kh2)=gl(kh3,kz,kh2)
      gr(kh3,kh2,kz)=gl(kh3,kh2,kz)

c---------------------------------------------- higgs bosons (trilinear)
      la3t=lam3+lam4+lam5

      ghhh(1,1,1)=-6.d0/g2weak*(lam1*ca3*cb+lam2*sa3*sb+
     &     la3t*sa*ca*sbpa+lam6*ca2*(3.d0*sa*cb+ca*sb)+
     &     lam7*sa2*(3.d0*ca*sb+sa*cb))
      gl(kh1,kh1,kh1)=dcmplx(ghhh(1,1,1),0.d0)
      gr(kh1,kh1,kh1)=gl(kh1,kh1,kh1)

      ghhh(1,1,2)=6.d0/g2weak*(lam1*ca2*cb*sa-lam2*sa2*sb*ca+
     &     la3t*(sa3*cb-ca3*sb+2.d0/3.d0*sbma)-
     &     lam6*ca*(cb*c2a-sa*sbpa)-lam7*ca*(sb*c2a+ca*sbpa))
      gl(kh1,kh1,kh2)=dcmplx(ghhh(1,1,2),0.d0)
      gl(kh1,kh2,kh1)=gl(kh1,kh1,kh2)
      gl(kh2,kh1,kh1)=gl(kh1,kh1,kh2)
      gr(kh1,kh1,kh2)=gl(kh1,kh1,kh2)
      gr(kh1,kh2,kh1)=gl(kh1,kh2,kh1)
      gr(kh2,kh1,kh1)=gl(kh2,kh1,kh1)

      ghhh(1,2,2)=-6.d0/g2weak*(lam1*sa2*cb*ca+lam2*ca2*sb*sa+
     &     la3t*(sa3*sb+ca3*cb-2.d0/3.d0*cbma)-
     &     lam6*sa*(cb*c2a+ca*cbpa)+lam7*ca*(sb*c2a+sa*cbpa))
      gl(kh1,kh2,kh2)=dcmplx(ghhh(1,2,2),0.d0)
      gl(kh2,kh1,kh2)=gl(kh1,kh2,kh2)
      gl(kh2,kh2,kh1)=gl(kh1,kh2,kh2)
      gr(kh1,kh2,kh2)=gl(kh1,kh2,kh2)
      gr(kh2,kh1,kh2)=gl(kh2,kh1,kh2)
      gr(kh2,kh2,kh1)=gl(kh2,kh2,kh1)

      ghhh(2,2,2)=6.d0/g2weak*(lam1*sa3*cb-lam2*ca3*sb+
     &     la3t*sa*ca*cbpa-lam6*sa2*(3.d0*ca*cb-sa*sb)+
     &     lam7*ca2*(3.d0*sa*sb-ca*cb))
      gl(kh2,kh2,kh2)=dcmplx(ghhh(2,2,2),0.d0)
      gr(kh2,kh2,kh2)=gl(kh2,kh2,kh2)

      ghaa(1)=-2.d0/g2weak*(lam1*sb2*cb*ca+lam2*cb2*sb*sa+
     &     la3t*(sb3*sa+cb3*ca)-2.d0*lam5*cbma-
     &     lam6*sb*(cb*cbpa+ca*c2b)-lam7*cb*(sb*cbpa+sa*c2b))
      gl(kh1,kh3,kh3)=dcmplx(ghaa(1),0.d0)
      gl(kh3,kh1,kh3)=gl(kh1,kh3,kh3)
      gl(kh3,kh3,kh1)=gl(kh1,kh3,kh3)
      gr(kh1,kh3,kh3)=gl(kh1,kh3,kh3)
      gr(kh3,kh1,kh3)=gl(kh3,kh1,kh3)
      gr(kh3,kh3,kh1)=gl(kh3,kh3,kh1)

      ghaa(2)=2.d0/g2weak*(lam1*sb2*cb*sa-lam2*cb2*sb*ca-
     &     la3t*(sb3*ca-cb3*sa)+2.d0*lam5*sbma-
     &     lam6*sb*(cb*sbpa+sa*c2b)-lam7*cb*(ca*c2b-sb*sbpa))
      gl(kh2,kh3,kh3)=dcmplx(ghaa(2),0.d0)
      gl(kh3,kh2,kh3)=gl(kh2,kh3,kh3)
      gl(kh3,kh3,kh2)=gl(kh2,kh3,kh3)
      gr(kh2,kh3,kh3)=gl(kh2,kh3,kh3)
      gr(kh3,kh2,kh3)=gl(kh3,kh2,kh3)
      gr(kh3,kh3,kh2)=gl(kh3,kh3,kh2)

      gh2h(1)=ghaa(1)-2.d0/g2weak*(lam5-lam4)*cbma
      gl(kh1,khc,khc)=dcmplx(gh2h(1),0.d0)
      gl(khc,kh1,khc)=gl(kh1,khc,khc)
      gl(khc,khc,kh1)=gl(kh1,khc,khc)
      gr(kh1,khc,khc)=gl(kh1,khc,khc)
      gr(khc,kh1,khc)=gl(khc,kh1,khc)
      gr(khc,khc,kh1)=gl(khc,khc,kh1)

      gh2h(2)=ghaa(2)-2.d0/g2weak*(lam5-lam4)*sbma
      gl(kh2,khc,khc)=dcmplx(gh2h(2),0.d0)
      gl(khc,kh2,khc)=gl(kh2,khc,khc)
      gl(khc,khc,kh2)=gl(kh2,khc,khc)
      gr(kh2,khc,khc)=gl(kh2,khc,khc)
      gr(khc,kh2,khc)=gl(khc,kh2,khc)
      gr(khc,khc,kh2)=gl(khc,khc,kh2)

      gl(kh2,kgold0,kgold0)=dcmplx(
     &     -g2weak*(mass(kh2)/mass(kw))**2*sbma/2.d0,0.d0)
      gl(kgold0,kh2,kgold0)=gl(kh2,kgold0,kgold0)
      gl(kgold0,kgold0,kh2)=gl(kh2,kgold0,kgold0)
      gr(kh2,kgold0,kgold0)=gl(kh2,kgold0,kgold0)
      gr(kgold0,kh2,kgold0)=gl(kh2,kgold0,kgold0)
      gr(kgold0,kgold0,kh2)=gl(kh2,kgold0,kgold0)

      gl(kh1,kgold0,kgold0)=dcmplx(
     &    -g2weak*(mass(kh1)/mass(kw))**2*cbma/2.d0,0.d0)
      gl(kgold0,kh1,kgold0)=gl(kh1,kgold0,kgold0)
      gl(kgold0,kgold0,kh1)=gl(kh1,kgold0,kgold0)
      gr(kh1,kgold0,kgold0)=gl(kh1,kgold0,kgold0)
      gr(kgold0,kh1,kgold0)=gl(kh1,kgold0,kgold0)
      gr(kgold0,kgold0,kh1)=gl(kh1,kgold0,kgold0)

      gl(kh2,kgoldc,kgoldc)=gl(kh2,kgold0,kgold0)
      gl(kgoldc,kh2,kgoldc)=gl(kh2,kgoldc,kgoldc)
      gl(kgoldc,kgoldc,kh2)=gl(kh2,kgoldc,kgoldc)
      gr(kh2,kgoldc,kgoldc)=gl(kh2,kgoldc,kgoldc)
      gr(kgoldc,kh2,kgoldc)=gl(kh2,kgoldc,kgoldc)
      gr(kgoldc,kgoldc,kh2)=gl(kh2,kgoldc,kgoldc)

      gl(kh1,kgoldc,kgoldc)=gl(kh1,kgold0,kgold0)
      gl(kgoldc,kh1,kgoldc)=gl(kh1,kgoldc,kgoldc)
      gl(kgoldc,kgoldc,kh1)=gl(kh1,kgoldc,kgoldc)
      gr(kh1,kgoldc,kgoldc)=gl(kh1,kgoldc,kgoldc)
      gr(kgoldc,kh1,kgoldc)=gl(kh1,kgoldc,kgoldc)
      gr(kgoldc,kgoldc,kh1)=gl(kh1,kgoldc,kgoldc)

      gl(kh2,kh3,kgold0)=dcmplx(-g2weak/2.d0*
     &     (mass(kh2)**2-mass(kh3)**2)/mass(kw)**2*cbma,0.d0)
      gl(kh3,kh2,kgold0)=gl(kh2,kh3,kgold0)
      gl(kh3,kgold0,kh2)=gl(kh2,kh3,kgold0)
      gl(kh2,kgold0,kh3)=gl(kh2,kh3,kgold0)
      gl(kgold0,kh2,kh3)=gl(kh2,kh3,kgold0)
      gl(kgold0,kh3,kh2)=gl(kh2,kh3,kgold0)
      gr(kh2,kh3,kgold0)=gl(kh2,kh3,kgold0)
      gr(kh3,kh2,kgold0)=gl(kh2,kh3,kgold0)
      gr(kh3,kgold0,kh2)=gl(kh2,kh3,kgold0)
      gr(kh2,kgold0,kh3)=gl(kh2,kh3,kgold0)
      gr(kgold0,kh2,kh3)=gl(kh2,kh3,kgold0)
      gr(kgold0,kh3,kh2)=gl(kh2,kh3,kgold0)

      gl(kh1,kh3,kgold0)=dcmplx(-g2weak/2.d0*
     &     (mass(kh1)**2-mass(kh3)**2)/mass(kw)**2*sbma,0.d0)
      gl(kh3,kh1,kgold0)=gl(kh1,kh3,kgold0)
      gl(kh3,kgold0,kh1)=gl(kh1,kh3,kgold0)
      gl(kh1,kgold0,kh3)=gl(kh1,kh3,kgold0)
      gl(kgold0,kh1,kh3)=gl(kh1,kh3,kgold0)
      gl(kgold0,kh3,kh1)=gl(kh1,kh3,kgold0)
      gr(kh1,kh3,kgold0)=gl(kh1,kh3,kgold0)
      gr(kh3,kh1,kgold0)=gl(kh1,kh3,kgold0)
      gr(kh3,kgold0,kh1)=gl(kh1,kh3,kgold0)
      gr(kh1,kgold0,kh3)=gl(kh1,kh3,kgold0)
      gr(kgold0,kh1,kh3)=gl(kh1,kh3,kgold0)
      gr(kgold0,kh3,kh1)=gl(kh1,kh3,kgold0)

      gl(kh2,khc,kgoldc)=dcmplx(g2weak/2.d0*
     &    (mass(khc)**2-mass(kh2)**2)/mass(kw)**2*cbma,0.d0)
      gl(khc,kh2,kgoldc)=gl(kh2,khc,kgoldc)
      gl(khc,kgoldc,kh2)=gl(kh2,khc,kgoldc)
      gl(kh2,kgoldc,khc)=gl(kh2,khc,kgoldc)
      gl(kgoldc,kh2,khc)=gl(kh2,kgoldc,khc)
      gl(kgoldc,khc,kh2)=gl(kh2,kgoldc,khc)
      gr(kh2,khc,kgoldc)=gl(kh2,khc,kgoldc)
      gr(khc,kh2,kgoldc)=gl(khc,kh2,kgoldc)
      gr(khc,kgoldc,kh2)=gl(khc,kgoldc,kh2)
      gr(kh2,kgoldc,khc)=gl(kh2,kgoldc,khc)
      gr(kgoldc,kh2,khc)=gl(kgoldc,kh2,khc)
      gr(kgoldc,khc,kh2)=gl(kgoldc,khc,kh2)

      gl(kh1,khc,kgoldc)=dcmplx(-g2weak/2.d0*
     &    (mass(khc)**2-mass(kh1)**2)/mass(kw)**2*sbma,0.d0)
      gl(khc,kh1,kgoldc)=gl(kh1,khc,kgoldc)
      gl(khc,kgoldc,kh1)=gl(kh1,khc,kgoldc)
      gl(kh1,kgoldc,khc)=gl(kh1,khc,kgoldc)
      gl(kgoldc,kh1,khc)=gl(kh1,kgoldc,khc)
      gl(kgoldc,khc,kh1)=gl(kh1,kgoldc,khc)
      gr(kh1,khc,kgoldc)=gl(kh1,khc,kgoldc)
      gr(khc,kh1,kgoldc)=gl(khc,kh1,kgoldc)
      gr(khc,kgoldc,kh1)=gl(khc,kgoldc,kh1)
      gr(kh1,kgoldc,khc)=gl(kh1,kgoldc,khc)
      gr(kgoldc,kh1,khc)=gl(kgoldc,kh1,khc)
      gr(kgoldc,khc,kh1)=gl(kgoldc,khc,kh1)

      gl(kh3,khc,kgoldc)=dcmplx(-g2weak/2.d0*
     &     (mass(khc)**2-mass(kh3)**2)/mass(kw)**2,0.d0)
      gl(khc,kh3,kgoldc)=gl(kh3,khc,kgoldc)
      gl(khc,kgoldc,kh3)=gl(kh3,khc,kgoldc)
      gl(kh3,kgoldc,khc)=-gl(kh3,khc,kgoldc)
      gl(kgoldc,kh3,khc)=gl(kh3,kgoldc,khc)
      gl(kgoldc,khc,kh3)=gl(kh3,kgoldc,khc)
      gr(kh3,khc,kgoldc)=gl(kh3,khc,kgoldc)
      gr(khc,kh3,kgoldc)=gl(khc,kh3,kgoldc)
      gr(khc,kgoldc,kh3)=gl(khc,kgoldc,kh3)
      gr(kh3,kgoldc,khc)=gl(kh3,kgoldc,khc)
      gr(kgoldc,kh3,khc)=gl(kgoldc,kh3,khc)
      gr(kgoldc,khc,kh3)=gl(kgoldc,khc,kh3)

c------------------------------------ fermions and gauge or higgs bosons

      do g=1,3

c     neutrinos
        vf=-g2o2cw*(0.5d0)
        af=-g2o2cw*(0.5d0)
        gl(kz,knu(g),knu(g))=dcmplx(vf+af,0.d0)
        gl(knu(g),kz,knu(g))=gl(kz,knu(g),knu(g))
        gl(knu(g),knu(g),kz)=gl(kz,knu(g),knu(g))

c     charged leptons
        aux=g2weak*rmass(kl(g))/(2.d0*mass(kw)*cosbe)

        gaff=+aux*sinbe ! corr 951108 cosbe->sinbe
        gl(kh3,kl(g),kl(g))=dcmplx(0.d0,-gaff)
        gl(kl(g),kh3,kl(g))=gl(kh3,kl(g),kl(g))
        gl(kl(g),kl(g),kh3)=gl(kh3,kl(g),kl(g))
        gr(kh3,kl(g),kl(g))=dcmplx(0.d0,gaff)
        gr(kl(g),kh3,kl(g))=gr(kh3,kl(g),kl(g))
        gr(kl(g),kl(g),kh3)=gr(kh3,kl(g),kl(g))

        ghff=-aux*ca ! corr 951108 sa->ca
        gl(kh1,kl(g),kl(g))=dcmplx(ghff,0.d0)
        gl(kl(g),kh1,kl(g))=gl(kh1,kl(g),kl(g))
        gl(kl(g),kl(g),kh1)=gl(kh1,kl(g),kl(g))
        gr(kh1,kl(g),kl(g))=dcmplx(ghff,0.d0)
        gr(kl(g),kh1,kl(g))=gr(kh1,kl(g),kl(g))
        gr(kl(g),kl(g),kh1)=gr(kh1,kl(g),kl(g))

        ghff=+aux*sa ! corr 951108 ca->-sa
        gl(kh2,kl(g),kl(g))=dcmplx(ghff,0.d0)
        gl(kl(g),kh2,kl(g))=gl(kh2,kl(g),kl(g))
        gl(kl(g),kl(g),kh2)=gl(kh2,kl(g),kl(g))
        gr(kh2,kl(g),kl(g))=dcmplx(ghff,0.d0)
        gr(kl(g),kh2,kl(g))=gr(kh2,kl(g),kl(g))
        gr(kl(g),kl(g),kh2)=gr(kh2,kl(g),kl(g))

        vf=-g2o2cw*(-0.5d0-2.0d0*echarg(ke)*sinthw**2)
        af=-g2o2cw*(-0.5d0)
        gl(kz,kl(g),kl(g))=vf+af
        gl(kl(g),kz,kl(g))=gl(kz,kl(g),kl(g))
        gl(kl(g),kl(g),kz)=gl(kz,kl(g),kl(g))
        gr(kz,kl(g),kl(g))=vf-af
        gr(kl(g),kz,kl(g))=gr(kz,kl(g),kl(g))
        gr(kl(g),kl(g),kz)=gr(kz,kl(g),kl(g))

        gl(kgamma,kl(g),kl(g))=-g2weak*sinthw*echarg(kl(g)) ! added je 960306
        gl(kl(g),kgamma,kl(g))=gl(kgamma,kl(g),kl(g))
        gl(kl(g),kl(g),kgamma)=gl(kgamma,kl(g),kl(g))
        gr(kgamma,kl(g),kl(g))=gl(kgamma,kl(g),kl(g))
        gr(kl(g),kgamma,kl(g))=gl(kl(g),kgamma,kl(g))
        gr(kl(g),kl(g),kgamma)=gl(kl(g),kl(g),kgamma)

c     up-quarks
        aux=g2weak*rmass(kqu(g))/(2.d0*mass(kw)*sinbe)

        gaff=+aux*cosbe ! corr 951108 sinbe->cosbe
        gl(kh3,kqu(g),kqu(g))=dcmplx(0.d0,-gaff)
        gl(kqu(g),kh3,kqu(g))=gl(kh3,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kh3)=gl(kh3,kqu(g),kqu(g))
        gr(kh3,kqu(g),kqu(g))=dcmplx(0.d0,gaff)
        gr(kqu(g),kh3,kqu(g))=gr(kh3,kqu(g),kqu(g))
        gr(kqu(g),kqu(g),kh3)=gr(kh3,kqu(g),kqu(g))

        ghff=-aux*sa ! corr 951108 ca->sa
        gl(kh1,kqu(g),kqu(g))=dcmplx(ghff,0.d0)
        gl(kqu(g),kh1,kqu(g))=gl(kh1,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kh1)=gl(kh1,kqu(g),kqu(g))
        gr(kh1,kqu(g),kqu(g))=dcmplx(ghff,0.d0)
        gr(kqu(g),kh1,kqu(g))=gr(kh1,kqu(g),kqu(g))
        gr(kqu(g),kqu(g),kh1)=gr(kh1,kqu(g),kqu(g))

        ghff=-aux*ca ! corr 951108 sa->-ca
        gl(kh2,kqu(g),kqu(g))=dcmplx(ghff,0.d0)
        gl(kqu(g),kh2,kqu(g))=gl(kh2,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kh2)=gl(kh2,kqu(g),kqu(g))
        gr(kh2,kqu(g),kqu(g))=dcmplx(ghff,0.d0)
        gr(kqu(g),kh2,kqu(g))=gr(kh2,kqu(g),kqu(g))
        gr(kqu(g),kqu(g),kh2)=gr(kh2,kqu(g),kqu(g))

        vf=-g2o2cw*(0.5d0-2.0d0*echarg(ku)*sinthw**2)
        af=-g2o2cw*(0.5d0)
        gl(kz,kqu(g),kqu(g))=dcmplx(vf+af,0.d0)
        gl(kqu(g),kz,kqu(g))=gl(kz,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kz)=gl(kz,kqu(g),kqu(g))
        gr(kz,kqu(g),kqu(g))=dcmplx(vf-af,0.d0)
        gr(kqu(g),kz,kqu(g))=gr(kz,kqu(g),kqu(g))
        gr(kqu(g),kqu(g),kz)=gr(kz,kqu(g),kqu(g))

        gl(kgamma,kqu(g),kqu(g))=-g2weak*sinthw*echarg(kqu(g)) ! je 960306
        gl(kqu(g),kgamma,kqu(g))=gl(kgamma,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kgamma)=gl(kgamma,kqu(g),kqu(g))
        gr(kgamma,kqu(g),kqu(g))=gl(kgamma,kqu(g),kqu(g))
        gr(kqu(g),kgamma,kqu(g))=gl(kqu(g),kgamma,kqu(g))
        gr(kqu(g),kqu(g),kgamma)=gl(kqu(g),kqu(g),kgamma)

c     gluon f f vertices. addition by pg 02-03-21
c       L = - g_s A^a_mu \bar{q}_j \gamma^mu T^a_{jk} q_k
        gl(kgluon,kqu(g),kqu(g))=-g3stro
        gl(kqu(g),kgluon,kqu(g))=gl(kgluon,kqu(g),kqu(g))
        gl(kqu(g),kqu(g),kgluon)=gl(kgluon,kqu(g),kqu(g))
        gr(kgluon,kqu(g),kqu(g))=gl(kgluon,kqu(g),kqu(g))
        gr(kqu(g),kgluon,kqu(g))=gl(kqu(g),kgluon,kqu(g))
        gr(kqu(g),kqu(g),kgluon)=gl(kqu(g),kqu(g),kgluon)

c     down-quarks
        aux=g2weak*rmass(kqd(g))/(2.d0*mass(kw)*cosbe)


        gaff=+aux*sinbe
        gl(kh3,kqd(g),kqd(g))=dcmplx(0.d0,-gaff)
        gl(kqd(g),kh3,kqd(g))=gl(kh3,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kh3)=gl(kh3,kqd(g),kqd(g))
        gr(kh3,kqd(g),kqd(g))=dcmplx(0.d0,gaff)    ! corr. je 960212 kqu->kqd
        gr(kqd(g),kh3,kqd(g))=gr(kh3,kqd(g),kqd(g))
        gr(kqd(g),kqd(g),kh3)=gr(kh3,kqd(g),kqd(g))

        ghff=-aux*ca
        gl(kh1,kqd(g),kqd(g))=dcmplx(ghff,0.d0)
        gl(kqd(g),kh1,kqd(g))=gl(kh1,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kh1)=gl(kh1,kqd(g),kqd(g))
        gr(kh1,kqd(g),kqd(g))=dcmplx(ghff,0.d0)
        gr(kqd(g),kh1,kqd(g))=gr(kh1,kqd(g),kqd(g))
        gr(kqd(g),kqd(g),kh1)=gr(kh1,kqd(g),kqd(g))

        ghff=+aux*sa
        gl(kh2,kqd(g),kqd(g))=dcmplx(ghff,0.d0)
        gl(kqd(g),kh2,kqd(g))=gl(kh2,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kh2)=gl(kh2,kqd(g),kqd(g))
        gr(kh2,kqd(g),kqd(g))=dcmplx(ghff,0.d0)
        gr(kqd(g),kh2,kqd(g))=gr(kh2,kqd(g),kqd(g))
        gr(kqd(g),kqd(g),kh2)=gr(kh2,kqd(g),kqd(g))

        vf=-g2o2cw*(-0.5d0-2.0d0*echarg(kd)*sinthw**2)
        af=-g2o2cw*(-0.5d0)
        gl(kz,kqd(g),kqd(g))=dcmplx(vf+af,0.d0)
        gl(kqd(g),kz,kqd(g))=gl(kz,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kz)=gl(kz,kqd(g),kqd(g))
        gr(kz,kqd(g),kqd(g))=dcmplx(vf-af,0.d0)
        gr(kqd(g),kz,kqd(g))=gr(kz,kqd(g),kqd(g))
        gr(kqd(g),kqd(g),kz)=gr(kz,kqd(g),kqd(g))

        gl(kgamma,kqd(g),kqd(g))=-g2weak*sinthw*echarg(kqd(g)) ! je 960306
        gl(kqd(g),kgamma,kqd(g))=gl(kgamma,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kgamma)=gl(kgamma,kqd(g),kqd(g))
        gr(kgamma,kqd(g),kqd(g))=gl(kgamma,kqd(g),kqd(g))
        gr(kqd(g),kgamma,kqd(g))=gl(kqd(g),kgamma,kqd(g))
        gr(kqd(g),kqd(g),kgamma)=gl(kqd(g),kqd(g),kgamma)

c     gluon f f vertices. addition by pg 02-03-21
c       L = - g_s A^a_mu \bar{q}_j \gamma^mu T^a_{jk} q_k
        gl(kgluon,kqd(g),kqd(g))=-g3stro
        gl(kqd(g),kgluon,kqd(g))=gl(kgluon,kqd(g),kqd(g))
        gl(kqd(g),kqd(g),kgluon)=gl(kgluon,kqd(g),kqd(g))
        gr(kgluon,kqd(g),kqd(g))=gl(kgluon,kqd(g),kqd(g))
        gr(kqd(g),kgluon,kqd(g))=gl(kqd(g),kgluon,kqd(g))
        gr(kqd(g),kqd(g),kgluon)=gl(kqd(g),kqd(g),kgluon)

c     w f f' vertices. addition by je 96-02-21
        gwff=-g2weak/sqrt(2.0d0)
c       neutrinos and charged leptons
        gl(kw,kl(g),knu(g))=dcmplx(gwff,0.0d0)
        gl(kl(g),kw,knu(g))=gl(kw,kl(g),knu(g))
        gl(kl(g),knu(g),kw)=gl(kw,kl(g),knu(g))
        gl(kw,knu(g),kl(g))=conjg(gl(kw,kl(g),knu(g)))
        gl(knu(g),kw,kl(g))=gl(kw,knu(g),kl(g))
        gl(knu(g),kl(g),kw)=gl(kw,knu(g),kl(g))

c       up- and down quarks
        do k=1,3
          gl(kw,kqd(g),kqu(k))=dcmplx(gwff,0.0d0)*ckm(k,g)
          gl(kqd(g),kw,kqu(k))=gl(kw,kqd(g),kqu(k))
          gl(kqd(g),kqu(k),kw)=gl(kw,kqd(g),kqu(k))
          gl(kw,kqu(k),kqd(g))=conjg(gl(kw,kqd(g),kqu(k)))
          gl(kqu(k),kw,kqd(g))=gl(kw,kqu(k),kqd(g))
          gl(kqu(k),kqd(g),kw)=gl(kw,kqu(k),kqd(g))
        enddo

c     h+ f f' verticies
        glhff=g2weak*rmass(kl(g))*tanbe/(sqrt(2.0d0)*mass(kw))
        grhff=0.0d0
c       neutrinos and charged leptons
        gl(khc,kl(g),knu(g))=dcmplx(glhff,0.0d0)
        gl(kl(g),khc,knu(g))=gl(khc,kl(g),knu(g))
        gl(kl(g),knu(g),khc)=gl(khc,kl(g),knu(g))
        gr(khc,knu(g),kl(g))=conjg(gl(khc,kl(g),knu(g)))
        gr(knu(g),khc,kl(g))=gr(khc,knu(g),kl(g))
        gr(knu(g),kl(g),khc)=gr(khc,knu(g),kl(g))

c       up- and down quarks
        do k=1,3
          glhff=g2weak*rmass(kqd(g))*tanbe/(sqrt(2.0d0)*mass(kw))
          grhff=g2weak*rmass(kqu(k))/(sqrt(2.0d0)*mass(kw)*tanbe)
          gl(khc,kqd(g),kqu(k))=dcmplx(glhff,0.0d0)*ckm(k,g)
          gl(kqd(g),khc,kqu(k))=gl(khc,kqd(g),kqu(k))
          gl(kqd(g),kqu(k),khc)=gl(khc,kqd(g),kqu(k))
          gr(khc,kqd(g),kqu(k))=dcmplx(grhff,0.0d0)*ckm(k,g)
          gr(kqd(g),khc,kqu(k))=gr(khc,kqd(g),kqu(k))
          gr(kqd(g),kqu(k),khc)=gr(khc,kqd(g),kqu(k))
          gl(khc,kqu(k),kqd(g))=conjg(gr(khc,kqd(g),kqu(k)))
          gl(kqu(k),khc,kqd(g))=gl(khc,kqu(k),kqd(g))
          gl(kqu(k),kqd(g),khc)=gl(khc,kqu(k),kqd(g))
          gr(khc,kqu(k),kqd(g))=conjg(gl(khc,kqd(g),kqu(k)))
          gr(kqu(k),khc,kqd(g))=gr(khc,kqu(k),kqd(g))
          gr(kqu(k),kqd(g),khc)=gr(khc,kqu(k),kqd(g))
        enddo

      enddo

c----------------------------------- neutralinos, fermions and sfermions
      do j=1,4

        do g=1,3

c       neutrinos
          gjll=-(g2weak*neunmx(j,2)-gyweak*neunmx(j,1))/sqrt(2.0d0)
          do k=1,3
            gl(ksnu(k),kn(j),knu(g))=conjg(gjll)*slulmx(k,g)
            gl(kn(j),ksnu(g),knu(g))=gl(ksnu(g),kn(j),knu(g))
            gl(kn(j),knu(g),ksnu(g))=gl(ksnu(g),kn(j),knu(g))
            gr(ksnu(k),knu(g),kn(j))=conjg(gl(ksnu(k),kn(j),knu(g)))
            gr(knu(g),ksnu(g),kn(j))=gr(ksnu(g),knu(g),kn(j))
            gr(knu(g),kn(j),ksnu(g))=gr(ksnu(g),knu(g),kn(j))
          enddo

c       charged leptons
          gjll=(g2weak*neunmx(j,2)+gyweak*neunmx(j,1))/sqrt(2.0d0)
          gjlr=-yukawa(kl(g))*neunmx(j,3)
          gjrl=-yukawa(kl(g))*neunmx(j,3)
          gjrr=sqrt(2.0d0)*echarg(ke)*gyweak*neunmx(j,1)
          do k=1,6
            gl(ksl(k),kn(j),kl(g))=
     &            conjg(gjll)*sldlmx(k,g)+conjg(gjrl)*sldrmx(k,g)
            gl(kn(j),ksl(k),kl(g))=gl(ksl(k),kn(j),kl(g))
            gl(kn(j),kl(g),ksl(k))=gl(ksl(k),kn(j),kl(g))
            gr(ksl(k),kn(j),kl(g))=
     &           gjlr*sldlmx(k,g)+gjrr*sldrmx(k,g)
            gr(kn(j),ksl(k),kl(g))=gr(ksl(k),kn(j),kl(g))
            gr(kn(j),kl(g),ksl(k))=gr(ksl(k),kn(j),kl(g))
            gl(ksl(k),kl(g),kn(j))=conjg(gr(ksl(k),kn(j),kl(g)))
            gl(kl(g),ksl(k),kn(j))=gl(ksl(k),kl(g),kn(j))
            gl(kl(g),kn(j),ksl(k))=gl(ksl(k),kl(g),kn(j))
            gr(ksl(k),kl(g),kn(j))=conjg(gl(ksl(k),kn(j),kl(g)))
            gr(kl(g),ksl(k),kn(j))=gr(ksl(k),kl(g),kn(j))
            gr(kl(g),kn(j),ksl(k))=gr(ksl(k),kl(g),kn(j))
          enddo

c       up-type quarks
          gjll=-(g2weak*neunmx(j,2)+
     &     gyweak*neunmx(j,1)/3.0d0)/sqrt(2.0d0)
          gjlr=-yukawa(kqu(g))*neunmx(j,4)
          gjrl=-yukawa(kqu(g))*neunmx(j,4)
          gjrr=sqrt(2.0d0)*echarg(ku)*gyweak*neunmx(j,1)
          do k=1,6
            gl(ksqu(k),kn(j),kqu(g))=
     &           conjg(gjll)*squlmx(k,g)+conjg(gjrl)*squrmx(k,g)
            gl(kn(j),ksqu(k),kqu(g))=gl(ksqu(k),kn(j),kqu(g))
            gl(kn(j),kqu(g),ksqu(k))=gl(ksqu(k),kn(j),kqu(g))
            gr(ksqu(k),kn(j),kqu(g))=
     &           gjlr*squlmx(k,g)+gjrr*squrmx(k,g)
            gr(kn(j),ksqu(k),kqu(g))=gr(ksqu(k),kn(j),kqu(g))
            gr(kn(j),kqu(g),ksqu(k))=gr(ksqu(k),kn(j),kqu(g))
            gl(ksqu(k),kqu(g),kn(j))=conjg(gr(ksqu(k),kn(j),kqu(g)))
            gl(kqu(g),ksqu(k),kn(j))=gl(ksqu(k),kqu(g),kn(j))
            gl(kqu(g),kn(j),ksqu(k))=gl(ksqu(k),kqu(g),kn(j))
            gr(ksqu(k),kqu(g),kn(j))=conjg(gl(ksqu(k),kn(j),kqu(g)))
            gr(kqu(g),ksqu(k),kn(j))=gr(ksqu(k),kqu(g),kn(j))
            gr(kqu(g),kn(j),ksqu(k))=gr(ksqu(k),kqu(g),kn(j))
          enddo

c       down-type quarks
          gjll=-(-g2weak*neunmx(j,2)+
     &     gyweak*neunmx(j,1)/3.0d0)/sqrt(2.0d0)
          gjlr=-yukawa(kqd(g))*neunmx(j,3)
          gjrl=-yukawa(kqd(g))*neunmx(j,3)
          gjrr=sqrt(2.0d0)*echarg(kd)*gyweak*neunmx(j,1)
          do k=1,6
            gl(ksqd(k),kn(j),kqd(g))=
     &           conjg(gjll)*sqdlmx(k,g)+conjg(gjrl)*sqdrmx(k,g)
            gl(kn(j),ksqd(k),kqd(g))=gl(ksqd(k),kn(j),kqd(g))
            gl(kn(j),kqd(g),ksqd(k))=gl(ksqd(k),kn(j),kqd(g))
            gr(ksqd(k),kn(j),kqd(g))=
     &           gjlr*sqdlmx(k,g)+gjrr*sqdrmx(k,g)
            gr(kn(j),ksqd(k),kqd(g))=gr(ksqd(k),kn(j),kqd(g))
            gr(kn(j),kqd(g),ksqd(k))=gr(ksqd(k),kn(j),kqd(g))
            gl(ksqd(k),kqd(g),kn(j))=conjg(gr(ksqd(k),kn(j),kqd(g)))
            gl(kqd(g),ksqd(k),kn(j))=gl(ksqd(k),kqd(g),kn(j))
            gl(kqd(g),kn(j),ksqd(k))=gl(ksqd(k),kqd(g),kn(j))
            gr(ksqd(k),kqd(g),kn(j))=conjg(gl(ksqd(k),kn(j),kqd(g)))
            gr(kqd(g),ksqd(k),kn(j))=gr(ksqd(k),kqd(g),kn(j))
            gr(kqd(g),kn(j),ksqd(k))=gr(ksqd(k),kqd(g),kn(j))
          enddo

        enddo
      enddo

c--------------------------------- neutralinos and gauge or higgs bosons

      do j=1,4
        do i=1,4
          qij=0.5d0*(
     &         neunmx(i,3)*(g2weak*neunmx(j,2)-gyweak*neunmx(j,1))+
     &         neunmx(j,3)*(g2weak*neunmx(i,2)-gyweak*neunmx(i,1)))
          sij=0.5d0*(
     &         neunmx(i,4)*(g2weak*neunmx(j,2)-gyweak*neunmx(j,1))+
     &         neunmx(j,4)*(g2weak*neunmx(i,2)-gyweak*neunmx(i,1)))

          gaij=-qij*sinbe+sij*cosbe
          gl(kh3,kn(i),kn(j))=-dcmplx(0.d0,1.d0)*conjg(gaij) ! corr. je 960208
          gl(kn(i),kh3,kn(j))=gl(kh3,kn(i),kn(j))
          gl(kn(i),kn(j),kh3)=gl(kh3,kn(i),kn(j))
          gr(kh3,kn(i),kn(j))=+dcmplx(0.d0,1.d0)*gaij
          gr(kn(i),kh3,kn(j))=gr(kh3,kn(i),kn(j))
          gr(kn(i),kn(j),kh3)=gr(kh3,kn(i),kn(j))

          ghij=-qij*ca+sij*sa
          gl(kh1,kn(i),kn(j))=conjg(ghij)
          gl(kn(i),kh1,kn(j))=gl(kh1,kn(i),kn(j))
          gl(kn(i),kn(j),kh1)=gl(kh1,kn(i),kn(j))
          gr(kh1,kn(i),kn(j))=ghij
          gr(kn(i),kh1,kn(j))=gr(kh1,kn(i),kn(j))
          gr(kn(i),kn(j),kh1)=gr(kh1,kn(i),kn(j))

          ghij=qij*sa+sij*ca
          gl(kh2,kn(i),kn(j))=conjg(ghij)
          gl(kn(i),kh2,kn(j))=gl(kh2,kn(i),kn(j))
          gl(kn(i),kn(j),kh2)=gl(kh2,kn(i),kn(j))
          gr(kh2,kn(i),kn(j))=ghij
          gr(kn(i),kh2,kn(j))=gr(kh2,kn(i),kn(j))
          gr(kn(i),kn(j),kh2)=gr(kh2,kn(i),kn(j))

          gzij=g2o2cw*(conjg(neunmx(i,3))*neunmx(j,3)-
     &                 conjg(neunmx(i,4))*neunmx(j,4))
          gl(kz,kn(i),kn(j))=-conjg(gzij)
          gl(kn(i),kz,kn(j))=gl(kz,kn(i),kn(j))
          gl(kn(i),kn(j),kz)=gl(kz,kn(i),kn(j))
          gr(kz,kn(i),kn(j))=gzij
          gr(kn(i),kz,kn(j))=gr(kz,kn(i),kn(j))
          gr(kn(i),kn(j),kz)=gr(kz,kn(i),kn(j))

        enddo
      enddo

c---------------------- neutralinos, charginos and gauge or higgs bosons
      do j=1,4
        do c=1,2

          glwnc=g2weak*(-neunmx(j,4)*conjg(chavmx(c,2))/sqrt(2.0d0)+
     &     neunmx(j,2)*conjg(chavmx(c,1)))
          grwnc=g2weak*(conjg(neunmx(j,3))*chaumx(c,2)/sqrt(2.0d0)+
     &     conjg(neunmx(j,2))*chaumx(c,1))
          gl(kw,kn(j),kcha(c))=glwnc
          gl(kn(j),kw,kcha(c))=gl(kw,kn(j),kcha(c))
          gl(kn(j),kcha(c),kw)=gl(kw,kn(j),kcha(c))
          gr(kw,kn(j),kcha(c))=grwnc
          gr(kn(j),kw,kcha(c))=gr(kw,kn(j),kcha(c))
          gr(kn(j),kcha(c),kw)=gr(kw,kn(j),kcha(c))
          gl(kw,kcha(c),kn(j))=conjg(glwnc)
          gl(kcha(c),kw,kn(j))=gl(kw,kcha(c),kn(j))
          gl(kcha(c),kn(j),kw)=gl(kw,kcha(c),kn(j))
          gr(kw,kcha(c),kn(j))=conjg(grwnc)
          gr(kcha(c),kw,kn(j))=gr(kw,kcha(c),kn(j))
          gr(kcha(c),kn(j),kw)=gr(kw,kcha(c),kn(j))

          glhnc=-cosbe*conjg(g2weak*neunmx(j,4)*chavmx(c,1)+
     &     (g2weak*neunmx(j,2)+gyweak*neunmx(j,1))*
     &         chavmx(c,2)/sqrt(2.0d0))
          grhnc=-sinbe*(g2weak*neunmx(j,3)*chaumx(c,1)-
     &     (g2weak*neunmx(j,2)+gyweak*neunmx(j,1))*
     &         chaumx(c,2)/sqrt(2.0d0))
          gl(khc,kn(j),kcha(c))=glhnc
          gl(kn(j),khc,kcha(c))=gl(khc,kn(j),kcha(c))
          gl(kn(j),kcha(c),khc)=gl(khc,kn(j),kcha(c))
          gr(khc,kn(j),kcha(c))=grhnc
          gr(kn(j),khc,kcha(c))=gr(khc,kn(j),kcha(c))
          gr(kn(j),kcha(c),khc)=gr(khc,kn(j),kcha(c))
          gl(khc,kcha(c),kn(j))=conjg(grhnc)
          gl(kcha(c),khc,kn(j))=gl(khc,kcha(c),kn(j))
          gl(kcha(c),kn(j),khc)=gl(khc,kcha(c),kn(j))
          gr(khc,kcha(c),kn(j))=conjg(glhnc)
          gr(kcha(c),khc,kn(j))=gr(khc,kcha(c),kn(j))
          gr(kcha(c),kn(j),khc)=gr(khc,kcha(c),kn(j))

        enddo

      enddo

c---------------------------------------- sfermions, charginos, fermions

      do k=1,6
        do c=1,2
          do g=1,3

c     charged sleptons (ie neutrinos)  ! corr. je 960313 no sum!!
            glscf =
     &        conjg(chaumx(c,2))*sldrmx(k,g)*yukawa(kl(g))-
     &        g2weak*conjg(chaumx(c,1))*sldlmx(k,g)
            grscf = dcmplx(0.d0,0.d0)
c             corr. je 960312 ksqd->ksl, kqu->knu
            gl(ksl(k),kcha(c),knu(g))=glscf
            gl(kcha(c),ksl(k),knu(g))=gl(ksl(k),kcha(c),knu(g))
            gl(kcha(c),knu(g),ksl(k))=gl(ksl(k),kcha(c),knu(g))
            gr(ksl(k),kcha(c),knu(g))=grscf
            gr(kcha(c),ksl(k),knu(g))=gr(ksl(k),kcha(c),knu(g))
            gr(kcha(c),kqd(g),ksl(k))=gr(ksl(k),kcha(c),knu(g))
            gl(ksl(k),knu(g),kcha(c))=conjg(gr(ksl(k),kcha(c),knu(g)))
            gl(knu(g),ksl(k),kcha(c))=gl(ksl(k),knu(g),kcha(c))
            gl(knu(g),kcha(c),ksl(k))=gl(ksl(k),knu(g),kcha(c))
            gr(ksl(k),knu(g),kcha(c))=conjg(gl(ksl(k),kcha(c),knu(g)))
            gr(knu(g),ksl(k),kcha(c))=gr(ksl(k),knu(g),kcha(c))
            gr(knu(g),kcha(c),ksl(k))=gr(ksl(k),knu(g),kcha(c))

            if (k.le.3) then
c...The sign convention for the following changed by je 020710.
c     sneutrinos (ie charged leptons)  ! corr. je 960313  no sum!!
c     attention: this is special
c        \tilde{\nu}^* \bar{\chi^{c.c.}} (glscf p_l + grscf p_r) e
c        \tilde{\nu} \bar{e} (glsfc p_l + grsfc p_r) \chi^{c.c}
c     the first case corresponds to both fermions incoming and is fig. 23
c     in gunion and haber with a minus sign.
c     the second case corresponds to both fermions outgoing ans is fig.
c     22 in gunion and haber.
c     The corresponding Feynman rules are
c     A. Both fermions outgoing:
c        (g^L_{sneutrino,lepton,chargino) P_L 
c         + g^R_{sneutrino,lepton,chargino} P_R) C
c 
c     B. Both fermions ingoing:
c        -C^{-1}(g^L_{sneutrino,chargino,lepton) P_L 
c                + g^R_{sneutrino,chargino,lepton} P_R)

              glscf =
     &          -g2weak*conjg(chavmx(c,1))*slulmx(k,g) ! corr. je 020710
              grscf =
     &             +chaumx(c,2)*slulmx(k,g)*yukawa(kl(g)) ! corr. je 020710
c               corr. je 960312 ksqu->ksnu, kqd->kl
              gl(ksnu(k),kcha(c),kl(g))=glscf
              gl(kcha(c),ksnu(k),kl(g))=gl(ksnu(k),kcha(c),kl(g))
              gl(kcha(c),kl(g),ksnu(k))=gl(ksnu(k),kcha(c),kl(g))
              gr(ksnu(k),kcha(c),kl(g))=grscf
              gr(kcha(c),ksnu(k),kl(g))=gr(ksnu(k),kcha(c),kl(g))
              gr(kcha(c),kl(g),ksnu(k))=gr(ksnu(k),kcha(c),kl(g))
              gl(ksnu(k),kl(g),kcha(c))=
     &          conjg(gr(ksnu(k),kcha(c),kl(g)))
              gl(kl(g),ksnu(k),kcha(c))=gl(ksnu(k),kl(g),kcha(c))
              gl(kl(g),kcha(c),ksnu(k))=gl(ksnu(k),kl(g),kcha(c))
              gr(ksnu(k),kl(g),kcha(c))=
     &          conjg(gl(ksnu(k),kcha(c),kl(g)))
              gr(kl(g),ksnu(k),kcha(c))=gr(ksnu(k),kl(g),kcha(c))
              gr(kl(g),kcha(c),ksnu(k))=gr(ksnu(k),kl(g),kcha(c))
            endif

c     down-type squarks (ie up-type quarks)
            glscf = dcmplx(0.d0,0.d0)
            grscf = dcmplx(0.d0,0.d0)
            do j=1,3
              glscf = glscf +
     &          (conjg(chaumx(c,2))*sqdrmx(k,j)*yukawa(kqd(j))-
     &           g2weak*conjg(chaumx(c,1))*sqdlmx(k,j))*conjg(ckm(g,j))
              grscf = grscf +
     &          chavmx(c,2)*sqdlmx(k,j)*conjg(ckm(g,j))*yukawa(kqu(g))
            enddo
            gl(ksqd(k),kcha(c),kqu(g))=glscf
            gl(kcha(c),ksqd(k),kqu(g))=gl(ksqd(k),kcha(c),kqu(g))
            gl(kcha(c),kqu(g),ksqd(k))=gl(ksqd(k),kcha(c),kqu(g))
            gr(ksqd(k),kcha(c),kqu(g))=grscf
            gr(kcha(c),ksqd(k),kqu(g))=gr(ksqd(k),kcha(c),kqu(g))
            gr(kcha(c),kqd(g),ksqd(k))=gr(ksqd(k),kcha(c),kqu(g))
            gl(ksqd(k),kqu(g),kcha(c))=conjg(gr(ksqd(k),kcha(c),kqu(g)))
            gl(kqu(g),ksqd(k),kcha(c))=gl(ksqd(k),kqu(g),kcha(c))
            gl(kqu(g),kcha(c),ksqd(k))=gl(ksqd(k),kqu(g),kcha(c))
            gr(ksqd(k),kqu(g),kcha(c))=conjg(gl(ksqd(k),kcha(c),kqu(g)))
            gr(kqu(g),ksqd(k),kcha(c))=gr(ksqd(k),kqu(g),kcha(c))
            gr(kqu(g),kcha(c),ksqd(k))=gr(ksqd(k),kqu(g),kcha(c))

c...The sign convention for the following changed by je 020710.
c     up-type squarks (ie down-type quarks) ! corr. je 960312
c     attention: this is special
c        \tilde{u}^* \bar{\chi^{c.c.}} (glscf p_l + grscf p_r) d
c        \tilde{u} \bar{d} (glsfc p_l + grsfc p_r) \chi^{c.c}
c     the first case corresponds to both fermions incoming and is fig. 23
c     in gunion and haber with a minus sign.
c     the second case corresponds to both fermions outgoing ans is fig.
c     22 in gunion and haber.
c     The corresponding Feynman rules are
c     A. Both fermions outgoing:
c        (g^L_{squark,quark,chargino) P_L 
c         + g^R_{squark,quark,chargino} P_R) C
c 
c     B. Both fermions ingoing:
c        -C^{-1}(g^L_{squark,chargino,quark) P_L 
c                + g^R_{squark,chargino,quark} P_R)

            glscf = dcmplx(0.d0,0.d0)
            grscf = dcmplx(0.d0,0.d0)
            do j=1,3
              glscf = glscf +
     &          (conjg(chavmx(c,2))*squrmx(k,j)*yukawa(kqu(j))-
     &           g2weak*conjg(chavmx(c,1))*squlmx(k,j))*ckm(j,g)
              grscf = grscf +
     &          chaumx(c,2)*squlmx(k,j)*ckm(j,g)*yukawa(kqd(g))
            enddo
            gl(ksqu(k),kcha(c),kqd(g))=glscf
            gl(kcha(c),ksqu(k),kqd(g))=gl(ksqu(k),kcha(c),kqd(g))
            gl(kcha(c),kqd(g),ksqu(k))=gl(ksqu(k),kcha(c),kqd(g))
            gr(ksqu(k),kcha(c),kqd(g))=grscf
            gr(kcha(c),ksqu(k),kqd(g))=gr(ksqu(k),kcha(c),kqd(g))
            gr(kcha(c),kqd(g),ksqu(k))=gr(ksqu(k),kcha(c),kqd(g))
            gl(ksqu(k),kqd(g),kcha(c))=
     &        conjg(gr(ksqu(k),kcha(c),kqd(g)))
            gl(kqd(g),ksqu(k),kcha(c))=gl(ksqu(k),kqd(g),kcha(c))
            gl(kqd(g),kcha(c),ksqu(k))=gl(ksqu(k),kqd(g),kcha(c))
            gr(ksqu(k),kqd(g),kcha(c))=
     &        conjg(gl(ksqu(k),kcha(c),kqd(g)))
            gr(kqd(g),ksqu(k),kcha(c))=gr(ksqu(k),kqd(g),kcha(c))
            gr(kqd(g),kcha(c),ksqu(k))=gr(ksqu(k),kqd(g),kcha(c))

          enddo
        enddo
      enddo


c------------------------------------------ sfermions, gluinos, fermions

      do k=1,6
        do g=1,3

c     up-type (s)quarks
          gl(ksqu(k),kgluin,kqu(g)) = -sqrt(2.d0)*g3stro*squlmx(k,g)
          gl(kgluin,ksqu(k),kqu(g)) = gl(ksqu(k),kgluin,kqu(g))
          gl(kgluin,kqu(g),ksqu(k)) = gl(ksqu(k),kgluin,kqu(g))
          gr(ksqu(k),kgluin,kqu(g)) = +sqrt(2.d0)*g3stro*squrmx(k,g)
          gr(kgluin,ksqu(k),kqu(g)) = gr(ksqu(k),kgluin,kqu(g))
          gr(kgluin,kqu(g),ksqu(k)) = gr(ksqu(k),kgluin,kqu(g))
          gl(ksqu(k),kqu(g),kgluin) = conjg(gr(ksqu(k),kgluin,kqu(g)))
          gl(kqu(g),ksqu(k),kgluin) = gl(ksqu(k),kqu(g),kgluin)
          gl(kqu(g),kgluin,ksqu(k)) = gl(ksqu(k),kqu(g),kgluin)
          gr(ksqu(k),kqu(g),kgluin) = conjg(gl(ksqu(k),kgluin,kqu(g)))
          gr(kqu(g),ksqu(k),kgluin) = gr(ksqu(k),kqu(g),kgluin)
          gr(kqu(g),kgluin,ksqu(k)) = gr(ksqu(k),kqu(g),kgluin)

c     down-type (s)quarks
          gl(ksqd(k),kgluin,kqd(g)) = -sqrt(2.d0)*g3stro*sqdlmx(k,g)
          gl(kgluin,ksqd(k),kqd(g)) = gl(ksqd(k),kgluin,kqd(g))
          gl(kgluin,kqd(g),ksqd(k)) = gl(ksqd(k),kgluin,kqd(g))
          gr(ksqd(k),kgluin,kqd(g)) = +sqrt(2.d0)*g3stro*sqdrmx(k,g)
          gr(kgluin,ksqd(k),kqd(g)) = gr(ksqd(k),kgluin,kqd(g))
          gr(kgluin,kqd(g),ksqd(k)) = gr(ksqd(k),kgluin,kqd(g))
          gl(ksqd(k),kqd(g),kgluin) = conjg(gr(ksqd(k),kgluin,kqd(g)))
          gl(kqd(g),ksqd(k),kgluin) = gl(ksqd(k),kqd(g),kgluin)
          gl(kqd(g),kgluin,ksqd(k)) = gl(ksqd(k),kqd(g),kgluin)
          gr(ksqd(k),kqd(g),kgluin) = conjg(gl(ksqd(k),kgluin,kqd(g)))
          gr(kqd(g),ksqd(k),kgluin) = gr(ksqd(k),kqd(g),kgluin)
          gr(kqd(g),kgluin,ksqd(k)) = gr(ksqd(k),kqd(g),kgluin)
        enddo
      enddo

c----------------------------------- charginos and gauge or higgs bosons

      do c=1,2                   ! added by je 96-02-23
        do d=1,2

          qcd=g2weak*conjg(chaumx(c,2)*chavmx(d,1))/sqrt(2.0d0)
          scd=g2weak*conjg(chaumx(c,1)*chavmx(d,2))/sqrt(2.0d0)
          qdc=g2weak*chaumx(d,2)*chavmx(c,1)/sqrt(2.0d0)
          sdc=g2weak*chaumx(d,1)*chavmx(c,2)/sqrt(2.0d0)

c     charginos and higgs bosons
          gl(kh1,kcha(c),kcha(d))=-qcd*ca-scd*sa
          gl(kcha(c),kh1,kcha(d))=gl(kh1,kcha(c),kcha(d))
          gl(kcha(c),kcha(d),kh1)=gl(kh1,kcha(c),kcha(d))
          gr(kh1,kcha(c),kcha(d))=-qdc*ca-sdc*sa
          gr(kcha(c),kh1,kcha(d))=gr(kh1,kcha(c),kcha(d))
          gr(kcha(c),kcha(d),kh1)=gr(kh1,kcha(c),kcha(d))

          gl(kh2,kcha(c),kcha(d))=qcd*sa-scd*ca
          gl(kcha(c),kh2,kcha(d))=gl(kh2,kcha(c),kcha(d))
          gl(kcha(c),kcha(d),kh2)=gl(kh2,kcha(c),kcha(d))
          gr(kh2,kcha(c),kcha(d))=qdc*sa-sdc*ca
          gr(kcha(c),kh2,kcha(d))=gr(kh2,kcha(c),kcha(d))
          gr(kcha(c),kcha(d),kh2)=gr(kh2,kcha(c),kcha(d))

          gl(kh3,kcha(c),kcha(d))=dcmplx(0.d0,1.d0)*
     &      (qcd*sinbe+scd*cosbe)
          gl(kcha(c),kh3,kcha(d))=gl(kh3,kcha(c),kcha(d))
          gl(kcha(c),kcha(d),kh3)=gl(kh3,kcha(c),kcha(d))
          gr(kh3,kcha(c),kcha(d))=-dcmplx(0.d0,1.d0)*
     &      (qdc*sinbe+sdc*cosbe)
          gr(kcha(c),kh3,kcha(d))=gr(kh3,kcha(c),kcha(d))
          gr(kcha(c),kcha(d),kh3)=gr(kh3,kcha(c),kcha(d))

c     charginos and gauge bosons (z)
          ocdl=-g2weak/costhw*(chavmx(c,1)*conjg(chavmx(d,1))
     &      +0.5d0*chavmx(c,2)*conjg(chavmx(d,2)))
          if (c.eq.d) ocdl=ocdl+g2weak/costhw*sinthw**2
          ocdr=-g2weak/costhw*(conjg(chaumx(c,1))*chaumx(d,1)
     &      +0.5d0*conjg(chaumx(c,2))*chaumx(d,2))
          if (c.eq.d) ocdr=ocdr+g2weak/costhw*sinthw**2

          gl(kz,kcha(c),kcha(d))=ocdl
          gl(kcha(c),kz,kcha(d))=gl(kz,kcha(c),kcha(d))
          gl(kcha(c),kcha(d),kz)=gl(kz,kcha(c),kcha(d))
          gr(kz,kcha(c),kcha(d))=ocdr
          gr(kcha(c),kz,kcha(d))=gr(kz,kcha(c),kcha(d))
          gr(kcha(c),kcha(d),kz)=gr(kz,kcha(c),kcha(d))

        enddo
      enddo

c     charginos and gamma
      do c=1,2
        gl(kgamma,kcha(c),kcha(c))=-dcmplx(g2weak*sinthw,0.0d0)
        gl(kcha(c),kgamma,kcha(c))=gl(kgamma,kcha(c),kcha(c))
        gl(kcha(c),kcha(c),kgamma)=gl(kgamma,kcha(c),kcha(c))
        gr(kgamma,kcha(c),kcha(c))=-dcmplx(g2weak*sinthw,0.0d0)
        gr(kcha(c),kgamma,kcha(c))=gr(kgamma,kcha(c),kcha(c))
        gr(kcha(c),kcha(c),kgamma)=gr(kgamma,kcha(c),kcha(c))
      enddo

c----------------------------------------- goldstone neutralino chargino

      do j=1,4
        do c=1,2

          glhnc=-sinbe*conjg(g2weak*neunmx(j,4)*chavmx(c,1)+
     &     (g2weak*neunmx(j,2)+gyweak*neunmx(j,1))*
     &          chavmx(c,2)/sqrt(2.0d0))
          grhnc= cosbe*(g2weak*neunmx(j,3)*chaumx(c,1)-
     &     (g2weak*neunmx(j,2)+gyweak*neunmx(j,1))*
     &         chaumx(c,2)/sqrt(2.0d0))

          gl(kgoldc,kn(j),kcha(c))=glhnc
          gl(kn(j),kgoldc,kcha(c))=gl(kgoldc,kn(j),kcha(c))
          gl(kn(j),kcha(c),kgoldc)=gl(kgoldc,kn(j),kcha(c))
          gr(kgoldc,kn(j),kcha(c))=grhnc
          gr(kn(j),kgoldc,kcha(c))=gr(kgoldc,kn(j),kcha(c))
          gr(kn(j),kcha(c),kgoldc)=gr(kgoldc,kn(j),kcha(c))
          gl(kgoldc,kcha(c),kn(j))=conjg(grhnc)
          gl(kcha(c),kgoldc,kn(j))=gl(kgoldc,kcha(c),kn(j))
          gl(kcha(c),kn(j),kgoldc)=gl(kgoldc,kcha(c),kn(j))
          gr(kgoldc,kcha(c),kn(j))=conjg(glhnc)
          gr(kcha(c),kgoldc,kn(j))=gr(kgoldc,kcha(c),kn(j))
          gr(kcha(c),kn(j),kgoldc)=gr(kgoldc,kcha(c),kn(j))

        enddo
      enddo

c-----------------------------------------------------squark-squark-higgs

      do j=1,6
         do k=1,6
c h1-u~-u~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh1,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh1,ksqu(j),ksqu(k)) = gl(kh1,ksqu(j),ksqu(k))
     &              -g2o2cw*mass(kz)*cbpa*(1.d0-2.d0*echarg(ku)*sw2)*
     &              squlmx(j,l)*conjg(squlmx(k,l))
     &              -g2o2cw*mass(kz)*cbpa*(2.d0*echarg(ku)*sw2)*
     &              squrmx(j,l)*conjg(squrmx(k,l))
     &              -g2weak*sa/sb/mass(kw)*(
     &              squrmx(j,l)*rmass(kqu(l))**2*conjg(squrmx(k,l))+
     &              squlmx(j,l)*rmass(kqu(l))**2*conjg(squlmx(k,l)))
     &              -g2weak/sb/2.d0/mass(kw)*(
     &              squrmx(j,l)*rmass(kqu(l))*(-ca*mu+sa*asoftu(l))*
     &              conjg(squlmx(k,l))+
     &              conjg(squrmx(k,l)*rmass(kqu(l))*
     &              (-ca*mu+sa*asoftu(l))*conjg(squlmx(j,l))))
            enddo
            gl(ksqu(j),kh1,ksqu(k)) = gl(kh1,ksqu(j),ksqu(k))
            gl(ksqu(j),ksqu(k),kh1) = gl(kh1,ksqu(j),ksqu(k))
            gr(kh1,ksqu(j),ksqu(k)) = gl(kh1,ksqu(j),ksqu(k))
            gr(ksqu(j),kh1,ksqu(k)) = gl(kh1,ksqu(j),ksqu(k))
            gr(ksqu(j),ksqu(k),kh1) = gl(kh1,ksqu(j),ksqu(k))
c h1-d~-d~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh1,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh1,ksqd(j),ksqd(k)) = gl(kh1,ksqd(j),ksqd(k))
     &              -g2o2cw*mass(kz)*cbpa*(-1.d0-2.d0*echarg(kd)*sw2)*
     &              sqdlmx(j,l)*conjg(sqdlmx(k,l))
     &              -g2o2cw*mass(kz)*cbpa*(2.d0*echarg(kd)*sw2)*
     &              sqdrmx(j,l)*conjg(sqdrmx(k,l))
     &              -g2weak*ca/cb/mass(kw)*(
     &              sqdrmx(j,l)*rmass(kqd(l))**2*conjg(sqdrmx(k,l))+
     &              sqdlmx(j,l)*rmass(kqd(l))**2*conjg(sqdlmx(k,l)))
     &              -g2weak/cb/2.d0/mass(kw)*(
     &              sqdrmx(j,l)*rmass(kqd(l))*(-sa*mu+ca*asoftd(l))*
     &              conjg(sqdlmx(k,l))+
     &              conjg(sqdrmx(k,l)*rmass(kqd(l))*
     &              (-sa*mu+ca*asoftd(l))*conjg(sqdlmx(j,l))))
            enddo
            gl(ksqd(j),kh1,ksqd(k)) = gl(kh1,ksqd(j),ksqd(k))
            gl(ksqd(j),ksqd(k),kh1) = gl(kh1,ksqd(j),ksqd(k))
            gr(kh1,ksqd(j),ksqd(k)) = gl(kh1,ksqd(j),ksqd(k))
            gr(ksqd(j),kh1,ksqd(k)) = gl(kh1,ksqd(j),ksqd(k))
            gr(ksqd(j),ksqd(k),kh1) = gl(kh1,ksqd(j),ksqd(k))
c h2-u~-u~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh2,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh2,ksqu(j),ksqu(k)) = gl(kh2,ksqu(j),ksqu(k))
     &              +g2o2cw*mass(kz)*sbpa*(1.d0-2.d0*echarg(ku)*sw2)*
     &              squlmx(j,l)*conjg(squlmx(k,l))
     &              +g2o2cw*mass(kz)*sbpa*(2.d0*echarg(ku)*sw2)*
     &              squrmx(j,l)*conjg(squrmx(k,l))
     &              -g2weak*ca/sb/mass(kw)*(
     &              squrmx(j,l)*rmass(kqu(l))**2*conjg(squrmx(k,l))+
     &              squlmx(j,l)*rmass(kqu(l))**2*conjg(squlmx(k,l)))
     &              +g2weak/sb/2.d0/mass(kw)*(
     &              squrmx(j,l)*rmass(kqu(l))*(-sa*mu-ca*asoftu(l))*
     &              conjg(squlmx(k,l))+
     &              conjg(squrmx(k,l)*rmass(kqu(l))*
     &              (-sa*mu-ca*asoftu(l))*conjg(squlmx(j,l))))
            enddo
            gl(ksqu(j),kh2,ksqu(k)) = gl(kh2,ksqu(j),ksqu(k))
            gl(ksqu(j),ksqu(k),kh2) = gl(kh2,ksqu(j),ksqu(k))
            gr(kh2,ksqu(j),ksqu(k)) = gl(kh2,ksqu(j),ksqu(k))
            gr(ksqu(j),kh2,ksqu(k)) = gl(kh2,ksqu(j),ksqu(k))
            gr(ksqu(j),ksqu(k),kh2) = gl(kh2,ksqu(j),ksqu(k))
c h2-d~-d~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh2,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh2,ksqd(j),ksqd(k)) = gl(kh2,ksqd(j),ksqd(k))
     &              +g2o2cw*mass(kz)*sbpa*(-1.d0-2.d0*echarg(kd)*sw2)*
     &              sqdlmx(j,l)*conjg(sqdlmx(k,l))
     &              +g2o2cw*mass(kz)*sbpa*(2.d0*echarg(kd)*sw2)*
     &              sqdrmx(j,l)*conjg(sqdrmx(k,l))
     &              +g2weak*sa/cb/mass(kw)*(
     &              sqdrmx(j,l)*rmass(kqd(l))**2*conjg(sqdrmx(k,l))+
     &              sqdlmx(j,l)*rmass(kqd(l))**2*conjg(sqdlmx(k,l)))
     &              -g2weak/cb/2.d0/mass(kw)*(
     &              sqdrmx(j,l)*rmass(kqd(l))*(-ca*mu-sa*asoftd(l))*
     &              conjg(sqdlmx(k,l))+
     &              conjg(sqdrmx(k,l)*rmass(kqd(l))*
     &              (-ca*mu-sa*asoftd(l))*conjg(sqdlmx(j,l))))
            enddo
            gl(ksqd(j),kh2,ksqd(k)) = gl(kh2,ksqd(j),ksqd(k))
            gl(ksqd(j),ksqd(k),kh2) = gl(kh2,ksqd(j),ksqd(k))
            gr(kh2,ksqd(j),ksqd(k)) = gl(kh2,ksqd(j),ksqd(k))
            gr(ksqd(j),kh2,ksqd(k)) = gl(kh2,ksqd(j),ksqd(k))
            gr(ksqd(j),ksqd(k),kh2) = gl(kh2,ksqd(j),ksqd(k))
         enddo
      enddo

c--------------------------------------------------- the end of dsvertx1

      return
      end

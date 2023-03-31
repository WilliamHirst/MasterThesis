      subroutine dsvertx3
c_______________________________________________________________________
c  some couplings used in sfermion coannihilations
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c  history:
c    0110XX-020618 paolo gondolo
c    020903 Mia Schelke, Higgs-sfermion-sfermion, A-terms sign-change
c=======================================================================
c
c  vertices included:
c     Z-f~-f~, gamma-f~-f~, gluon-f~-f~, W-f~-F~, h-f~-f~
c  quartic vertices included:
c     Z-Z-h-h, W-W-h-h, Z-W-h-h, W-W-f~-f~, 
c     gamma-gamma-f~-f~, Z-Z-f~-f~, gamma-Z-f~-f~, gamma-W-f~-F~,
c     gluon-gluon-f~-f~, gluon-W-f~-F~, gluon-gamma-f~-f~, gluon-Z-f~-f~,
c     h-h-f~-f~, h-goldstone-f~-f~, goldstone-goldstone-f~-f~,
c     h-h-h-h, h-h-h-goldstone, h-h-goldstone-goldstone, 
c     h-goldstone-goldstone-goldstone, 4 x goldstone
c
      implicit none
      include 'dsmssm.h'
      real*8 g2o2cw,ca,sa,c2a,s2a,ca2,sa2,ca3,sa3,beta,cb,sb,c2b,s2b,
     & cb2,sb2,cb3,sb3,cbpa,sbpa,cbma,sbma,e_charge,aux,sw2,cw2,tw2,gsq,
     & aux2
      real*8 rmass(12) ! temporary array for lepton and quark running masses
                       ! currently extracted from yukawa 
                       ! one day we will have rmass(k,q)
      integer i,j,k,l,m,n,a,b,c
      complex*16 
     & gmell(6,6),gmull(6,6),gmdll(6,6),
     & gmerr(6,6),gmurr(6,6),gmdrr(6,6),
     & gmnll(3,3),gmudll(6,6),gmnell(3,6),
     & gm2ell(6,6),gm2ull(6,6),gm2dll(6,6),
     & gm2err(6,6),gm2urr(6,6),gm2drr(6,6),
     & gme2nll(3,3),gmd2ull(6,6),gmu2dll(6,6)
      integer hj(4)
      real*8 auxl1,auxl2,auxr1,auxr2,cc(4),dd(4,2)
      complex*16 auxc,auxc1,auxc2,
     &     rr(4),ss(4),tt(4),vv(4),rp(4),sp(4),tp(4),vp(4)
c...... temporary fix for Asoft (PG 2002-02-16)
      complex*16 au(3,3),ad(3,3),ae(3,3)
      do i=1,3
         do j=1,3
            au(i,j)=dcmplx(0.d0,0.d0)
            ad(i,j)=dcmplx(0.d0,0.d0)
            ae(i,j)=dcmplx(0.d0,0.d0)
         enddo
      enddo
      au(1,1)=asoftu(1)
      au(2,2)=asoftu(2)
      au(3,3)=asoftu(3)
      ad(1,1)=asoftd(1)
      ad(2,2)=asoftd(2)
      ad(3,3)=asoftd(3)
      ae(1,1)=asofte(1)
      ae(2,2)=asofte(2)
      ae(3,3)=asofte(3)
c...... end of temporary fix for Asoft
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
      gsq=g2weak**2
      e_charge=g2weak*sinthw
      sw2=sinthw**2
      cw2=costhw**2
      tw2=sw2/cw2
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
      cbpa=cos(beta+alpha)
      sbpa=sin(beta+alpha)
      cbma=cos(beta-alpha)
      sbma=sin(beta-alpha)

c------------------------------- matrices useful for sfermion couplings
      
      do i=1,6
         do j=1,6
            gmell(i,j)=dcmplx(0.d0,0.d0)
            gmull(i,j)=dcmplx(0.d0,0.d0)
            gmdll(i,j)=dcmplx(0.d0,0.d0)
            gmerr(i,j)=dcmplx(0.d0,0.d0)
            gmurr(i,j)=dcmplx(0.d0,0.d0)
            gmdrr(i,j)=dcmplx(0.d0,0.d0)
            gm2ell(i,j)=dcmplx(0.d0,0.d0)
            gm2ull(i,j)=dcmplx(0.d0,0.d0)
            gm2dll(i,j)=dcmplx(0.d0,0.d0)
            gm2err(i,j)=dcmplx(0.d0,0.d0)
            gm2urr(i,j)=dcmplx(0.d0,0.d0)
            gm2drr(i,j)=dcmplx(0.d0,0.d0)
            gmd2ull(i,j)=dcmplx(0.d0,0.d0)
            gmu2dll(i,j)=dcmplx(0.d0,0.d0)
            do k=1,3
               gmell(i,j)=gmell(i,j)+sldlmx(i,k)*conjg(sldlmx(j,k))
               gmull(i,j)=gmull(i,j)+squlmx(i,k)*conjg(squlmx(j,k))
               gmdll(i,j)=gmdll(i,j)+sqdlmx(i,k)*conjg(sqdlmx(j,k))
               gmerr(i,j)=gmerr(i,j)+sldrmx(i,k)*conjg(sldrmx(j,k))
               gmurr(i,j)=gmurr(i,j)+squrmx(i,k)*conjg(squrmx(j,k))
               gmdrr(i,j)=gmdrr(i,j)+sqdrmx(i,k)*conjg(sqdrmx(j,k))
               gm2ell(i,j)=gm2ell(i,j)+
     &              sldlmx(i,k)*rmass(kl(k))**2*conjg(sldlmx(j,k))
               gm2ull(i,j)=gm2ull(i,j)+
     &              squlmx(i,k)*rmass(kqu(k))**2*conjg(squlmx(j,k))
               gm2dll(i,j)=gm2dll(i,j)+
     &              sqdlmx(i,k)*rmass(kqd(k))**2*conjg(sqdlmx(j,k))
               gm2err(i,j)=gm2err(i,j)+
     &              sldrmx(i,k)*rmass(kl(k))**2*conjg(sldrmx(j,k))
               gm2urr(i,j)=gm2urr(i,j)+
     &              squrmx(i,k)*rmass(kqu(k))**2*conjg(squrmx(j,k))
               gm2drr(i,j)=gm2drr(i,j)+
     &              sqdrmx(i,k)*rmass(kqd(k))**2*conjg(sqdrmx(j,k))
               do a=1,3
                  do b=1,3
                     gmu2dll(i,j)=gmu2dll(i,j)+
     &                    sqdlmx(i,k)*conjg(ckm(a,k))*
     &                    rmass(kqu(a))**2*ckm(a,b)*conjg(sqdlmx(j,b))
                     gmd2ull(i,j)=gmd2ull(i,j)+
     &                    squlmx(i,k)*ckm(k,a)*rmass(kqd(a))**2*
     &                    conjg(ckm(b,a))*conjg(squlmx(j,b))
                  enddo
               enddo
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            gmnll(i,j)=dcmplx(0.d0,0.d0)
            gme2nll(i,j)=dcmplx(0.d0,0.d0)
            do k=1,3
               gmnll(i,j)=gmnll(i,j)+slulmx(i,k)*conjg(slulmx(j,k))
               gme2nll(i,j)=gme2nll(i,j)+
     &              slulmx(i,k)*rmass(kl(k))**2*conjg(slulmx(j,k))
            enddo
         enddo
      enddo
      do i=1,6
         do j=1,6
            gmudll(i,j)=dcmplx(0.d0,0.d0)
            do k=1,3
               do l=1,3
                  gmudll(i,j)=gmudll(i,j)+
     &                 squlmx(i,k)*ckm(k,l)*conjg(sqdlmx(j,l))
               enddo
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,6
            gmnell(i,j)=dcmplx(0.d0,0.d0)
            do k=1,3
               gmnell(i,j)=gmnell(i,j)+slulmx(i,k)*conjg(sldlmx(j,k))
            enddo
         enddo
      enddo

c      do k=1,6
c         write (*,*) 'gmull',(gmull(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gmurr',(gmurr(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2ull',(gm2ull(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2urr',(gm2urr(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gmdll',(gmdll(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gmdrr',(gmdrr(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2dll',(gm2dll(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2drr',(gm2drr(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gmell',(gmell(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gmerr',(gmerr(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2ell',(gm2ell(k,l),l=1,6)
c      enddo
c      do k=1,6
c         write (*,*) 'gm2err',(gm2err(k,l),l=1,6)
c      enddo
c      do k=1,3
c         write (*,*) 'gmnll',(gmnll(k,l),l=1,3)
c      enddo
c      do k=1,6
c         write (*,*) 'gmudll',(gmudll(k,l),l=1,6)
c      enddo
c      do k=1,3
c         write (*,*) 'gmnell',(gmnell(k,l),l=1,6)
c      enddo

c-----------------------------------------------sfermion-sfermion-gauge
c     gl=gr ; L = gl(V,fj,fk) V fj* T-color i double-arrow-partial fk

      do j=1,6
         do k=1,6
c Z-u~-u~
            gl(kz,ksqu(j),ksqu(k)) = 
     &           -g2o2cw*(1.d0-2.d0*echarg(ku)*sw2)*gmull(j,k) 
     &           -g2o2cw*(-2.d0*echarg(ku)*sw2)*gmurr(j,k) 
            gl(ksqu(j),kz,ksqu(k)) = gl(kz,ksqu(j),ksqu(k))
            gl(ksqu(j),ksqu(k),kz) = gl(kz,ksqu(j),ksqu(k))
            gr(kz,ksqu(j),ksqu(k)) = gl(kz,ksqu(j),ksqu(k))
            gr(ksqu(j),kz,ksqu(k)) = gl(kz,ksqu(j),ksqu(k))
            gr(ksqu(j),ksqu(k),kz) = gl(kz,ksqu(j),ksqu(k))
c Z-d~-d~
            gl(kz,ksqd(j),ksqd(k)) = 
     &           -g2o2cw*(-1.d0-2.d0*echarg(kd)*sw2)*gmdll(j,k) 
     &           -g2o2cw*(-2.d0*echarg(kd)*sw2)*gmdrr(j,k) 
            gl(ksqd(j),kz,ksqd(k)) = gl(kz,ksqd(j),ksqd(k))
            gl(ksqd(j),ksqd(k),kz) = gl(kz,ksqd(j),ksqd(k))
            gr(kz,ksqd(j),ksqd(k)) = gl(kz,ksqd(j),ksqd(k))
            gr(ksqd(j),kz,ksqd(k)) = gl(kz,ksqd(j),ksqd(k))
            gr(ksqd(j),ksqd(k),kz) = gl(kz,ksqd(j),ksqd(k))
c Z-e~-e~
            gl(kz,ksl(j),ksl(k)) = 
     &           -g2o2cw*(-1.d0-2.d0*echarg(ke)*sw2)*gmell(j,k) 
     &           -g2o2cw*(-2.d0*echarg(ke)*sw2)*gmerr(j,k) 
            gl(ksl(j),kz,ksl(k)) = gl(kz,ksl(j),ksl(k))
            gl(ksl(j),ksl(k),kz) = gl(kz,ksl(j),ksl(k))
            gr(kz,ksl(j),ksl(k)) = gl(kz,ksl(j),ksl(k))
            gr(ksl(j),kz,ksl(k)) = gl(kz,ksl(j),ksl(k))
            gr(ksl(j),ksl(k),kz) = gl(kz,ksl(j),ksl(k))
c gamma-u~-u~
            if (j.eq.k) then
               gl(kgamma,ksqu(j),ksqu(k)) = -e_charge*echarg(ku)
               gl(ksqu(j),kgamma,ksqu(k)) = gl(kgamma,ksqu(j),ksqu(k))
               gl(ksqu(j),ksqu(k),kgamma) = gl(kgamma,ksqu(j),ksqu(k))
               gr(kgamma,ksqu(j),ksqu(k)) = gl(kgamma,ksqu(j),ksqu(k))
               gr(ksqu(j),kgamma,ksqu(k)) = gl(kgamma,ksqu(j),ksqu(k))
               gr(ksqu(j),ksqu(k),kgamma) = gl(kgamma,ksqu(j),ksqu(k))
            else
               gl(kgamma,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
               gl(ksqu(j),kgamma,ksqu(k)) = dcmplx(0.d0,0.d0)
               gl(ksqu(j),ksqu(k),kgamma) = dcmplx(0.d0,0.d0)
               gr(kgamma,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
               gr(ksqu(j),kgamma,ksqu(k)) = dcmplx(0.d0,0.d0)
               gr(ksqu(j),ksqu(k),kgamma) = dcmplx(0.d0,0.d0)
            endif
c gamma-d~-d~
            if (j.eq.k) then
               gl(kgamma,ksqd(j),ksqd(k)) = -e_charge*echarg(kd)
               gl(ksqd(j),kgamma,ksqd(k)) = gl(kgamma,ksqd(j),ksqd(k))
               gl(ksqd(j),ksqd(k),kgamma) = gl(kgamma,ksqd(j),ksqd(k))
               gr(kgamma,ksqd(j),ksqd(k)) = gl(kgamma,ksqd(j),ksqd(k))
               gr(ksqd(j),kgamma,ksqd(k)) = gl(kgamma,ksqd(j),ksqd(k))
               gr(ksqd(j),ksqd(k),kgamma) = gl(kgamma,ksqd(j),ksqd(k))
            else
               gl(kgamma,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
               gl(ksqd(j),kgamma,ksqd(k)) = dcmplx(0.d0,0.d0)
               gl(ksqd(j),ksqd(k),kgamma) = dcmplx(0.d0,0.d0)
               gr(kgamma,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
               gr(ksqd(j),kgamma,ksqd(k)) = dcmplx(0.d0,0.d0)
               gr(ksqd(j),ksqd(k),kgamma) = dcmplx(0.d0,0.d0)
            endif
c gamma-e~-e~
            if (j.eq.k) then
               gl(kgamma,ksl(j),ksl(k)) = -e_charge*echarg(ke)
               gl(ksl(j),kgamma,ksl(k)) = gl(kgamma,ksl(j),ksl(k))
               gl(ksl(j),ksl(k),kgamma) = gl(kgamma,ksl(j),ksl(k))
               gr(kgamma,ksl(j),ksl(k)) = gl(kgamma,ksl(j),ksl(k))
               gr(ksl(j),kgamma,ksl(k)) = gl(kgamma,ksl(j),ksl(k))
               gr(ksl(j),ksl(k),kgamma) = gl(kgamma,ksl(j),ksl(k))
            else
               gl(kgamma,ksl(j),ksl(k)) = dcmplx(0.d0,0.d0)
               gl(ksl(j),kgamma,ksl(k)) = dcmplx(0.d0,0.d0)
               gl(ksl(j),ksl(k),kgamma) = dcmplx(0.d0,0.d0)
               gr(kgamma,ksl(j),ksl(k)) = dcmplx(0.d0,0.d0)
               gr(ksl(j),kgamma,ksl(k)) = dcmplx(0.d0,0.d0)
               gr(ksl(j),ksl(k),kgamma) = dcmplx(0.d0,0.d0)
            endif
c gluon-u~-u~
            if (j.eq.k) then
               gl(kgluon,ksqu(j),ksqu(k)) = -g3stro
               gl(ksqu(j),kgluon,ksqu(k)) = gl(kgluon,ksqu(j),ksqu(k))
               gl(ksqu(j),ksqu(k),kgluon) = gl(kgluon,ksqu(j),ksqu(k))
               gr(kgluon,ksqu(j),ksqu(k)) = gl(kgluon,ksqu(j),ksqu(k))
               gr(ksqu(j),kgluon,ksqu(k)) = gl(kgluon,ksqu(j),ksqu(k))
               gr(ksqu(j),ksqu(k),kgluon) = gl(kgluon,ksqu(j),ksqu(k))
            else
               gl(kgluon,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
               gl(ksqu(j),kgluon,ksqu(k)) = dcmplx(0.d0,0.d0)
               gl(ksqu(j),ksqu(k),kgluon) = dcmplx(0.d0,0.d0)
               gr(kgluon,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
               gr(ksqu(j),kgluon,ksqu(k)) = dcmplx(0.d0,0.d0)
               gr(ksqu(j),ksqu(k),kgluon) = dcmplx(0.d0,0.d0)
            endif
c gluon-d~-d~
            if (j.eq.k) then
               gl(kgluon,ksqd(j),ksqd(k)) = -g3stro
               gl(ksqd(j),kgluon,ksqd(k)) = gl(kgluon,ksqd(j),ksqd(k))
               gl(ksqd(j),ksqd(k),kgluon) = gl(kgluon,ksqd(j),ksqd(k))
               gr(kgluon,ksqd(j),ksqd(k)) = gl(kgluon,ksqd(j),ksqd(k))
               gr(ksqd(j),kgluon,ksqd(k)) = gl(kgluon,ksqd(j),ksqd(k))
               gr(ksqd(j),ksqd(k),kgluon) = gl(kgluon,ksqd(j),ksqd(k))
            else
               gl(kgluon,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
               gl(ksqd(j),kgluon,ksqd(k)) = dcmplx(0.d0,0.d0)
               gl(ksqd(j),ksqd(k),kgluon) = dcmplx(0.d0,0.d0)
               gr(kgluon,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
               gr(ksqd(j),kgluon,ksqd(k)) = dcmplx(0.d0,0.d0)
               gr(ksqd(j),ksqd(k),kgluon) = dcmplx(0.d0,0.d0)
            endif
c gluon-e~-e~
            gl(kgluon,ksl(j),ksl(k)) = dcmplx(0.d0,0.d0)
            gl(ksl(j),kgluon,ksl(k)) = dcmplx(0.d0,0.d0)
            gl(ksl(j),ksl(k),kgluon) = dcmplx(0.d0,0.d0)
            gr(kgluon,ksl(j),ksl(k)) = dcmplx(0.d0,0.d0)
            gr(ksl(j),kgluon,ksl(k)) = dcmplx(0.d0,0.d0)
            gr(ksl(j),ksl(k),kgluon) = dcmplx(0.d0,0.d0)
c W-u~-d~
            gl(kw,ksqu(j),ksqd(k)) = 
     &           -g2weak/sqrt(2.d0)*gmudll(j,k) 
            gl(ksqu(j),kw,ksqd(k)) = gl(kw,ksqu(j),ksqd(k))
            gl(ksqu(j),ksqd(k),kw) = gl(kw,ksqu(j),ksqd(k))
            gr(kw,ksqu(j),ksqd(k)) = gl(kw,ksqu(j),ksqd(k))
            gr(ksqu(j),kw,ksqd(k)) = gl(kw,ksqu(j),ksqd(k))
            gr(ksqu(j),ksqd(k),kw) = gl(kw,ksqu(j),ksqd(k))
c W-d~-u~
            gl(kw,ksqd(j),ksqu(k)) = 
     &           -g2weak/sqrt(2.d0)*conjg(gmudll(k,j))
            gl(ksqd(j),kw,ksqu(k)) = gl(kw,ksqd(j),ksqu(k))
            gl(ksqd(j),ksqu(k),kw) = gl(kw,ksqd(j),ksqu(k))
            gr(kw,ksqd(j),ksqu(k)) = gl(kw,ksqd(j),ksqu(k))
            gr(ksqd(j),kw,ksqu(k)) = gl(kw,ksqd(j),ksqu(k))
            gr(ksqd(j),ksqu(k),kw) = gl(kw,ksqd(j),ksqu(k))
         enddo
      enddo
      do j=1,3
         do k=1,3
c Z-nu~-nu~
            gl(kz,ksnu(j),ksnu(k)) = -g2o2cw*gmnll(j,k) 
            gl(ksnu(j),kz,ksnu(k)) = gl(kz,ksnu(j),ksnu(k))
            gl(ksnu(j),ksnu(k),kz) = gl(kz,ksnu(j),ksnu(k))
            gr(kz,ksnu(j),ksnu(k)) = gl(kz,ksnu(j),ksnu(k))
            gr(ksnu(j),kz,ksnu(k)) = gl(kz,ksnu(j),ksnu(k))
            gr(ksnu(j),ksnu(k),kz) = gl(kz,ksnu(j),ksnu(k))
c gamma-nu~-nu~
            gl(kgamma,ksnu(j),ksnu(k)) = dcmplx(0.d0,0.d0)
            gl(ksnu(j),kgamma,ksnu(k)) = dcmplx(0.d0,0.d0)
            gl(ksnu(j),ksnu(k),kgamma) = dcmplx(0.d0,0.d0)
            gr(kgamma,ksnu(j),ksnu(k)) = dcmplx(0.d0,0.d0)
            gr(ksnu(j),kgamma,ksnu(k)) = dcmplx(0.d0,0.d0)
            gr(ksnu(j),ksnu(k),kgamma) = dcmplx(0.d0,0.d0)
c gluon-nu~-nu~
            gl(kgluon,ksnu(j),ksnu(k)) = dcmplx(0.d0,0.d0)
            gl(ksnu(j),kgluon,ksnu(k)) = dcmplx(0.d0,0.d0)
            gl(ksnu(j),ksnu(k),kgluon) = dcmplx(0.d0,0.d0)
            gr(kgluon,ksnu(j),ksnu(k)) = dcmplx(0.d0,0.d0)
            gr(ksnu(j),kgluon,ksnu(k)) = dcmplx(0.d0,0.d0)
            gr(ksnu(j),ksnu(k),kgluon) = dcmplx(0.d0,0.d0)
         enddo
         do k=1,6
c W-nu~-e~
            gl(kw,ksnu(j),ksl(k)) = 
     &           -g2weak/sqrt(2.d0)*gmnell(j,k) 
            gl(ksnu(j),kw,ksl(k)) = gl(kw,ksnu(j),ksl(k))
            gl(ksnu(j),ksl(k),kw) = gl(kw,ksnu(j),ksl(k))
            gr(kw,ksnu(j),ksl(k)) = gl(kw,ksnu(j),ksl(k))
            gr(ksnu(j),kw,ksl(k)) = gl(kw,ksnu(j),ksl(k))
            gr(ksnu(j),ksl(k),kw) = gl(kw,ksnu(j),ksl(k))
c W-e~-nu~
            gl(kw,ksl(k),ksnu(j)) = 
     &           -g2weak/sqrt(2.d0)*conjg(gmnell(j,k))
            gl(ksl(k),kw,ksnu(j)) = gl(kw,ksl(k),ksnu(j))
            gl(ksl(k),ksnu(j),kw) = gl(kw,ksl(k),ksnu(j))
            gr(kw,ksl(k),ksnu(j)) = gl(kw,ksl(k),ksnu(j))
            gr(ksl(k),kw,ksnu(j)) = gl(kw,ksl(k),ksnu(j))
            gr(ksl(k),ksnu(j),kw) = gl(kw,ksl(k),ksnu(j))
         enddo
      enddo

c------------------------------------------------sfermion-sfermion-higgs
c     gl=gr ; L = gl(h,sfj,sfk) h sfj* sfk

      do j=1,6
         do k=1,6
c h1-u~-u~
c     in dsvertx1
c h1-d~-d~
c     in dsvertx1
c h1-e~-e~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh1,ksl(j),ksl(k)) = -g2o2cw*mass(kz)*cbpa*(
     &           (-1.d0-2.d0*echarg(ke)*sw2)*gmell(j,k)+
     &           (2.d0*echarg(ke)*sw2)*gmerr(j,k))
            do l=1,3
               gl(kh1,ksl(j),ksl(k)) = gl(kh1,ksl(j),ksl(k))
     &              -g2weak*ca/cb/mass(kw)*(
     &              sldrmx(j,l)*rmass(kl(l))**2*conjg(sldrmx(k,l))+
     &              sldlmx(j,l)*rmass(kl(l))**2*conjg(sldlmx(k,l)))
     &              -g2weak/cb/2.d0/mass(kw)*(
     &              sldrmx(j,l)*rmass(kl(l))*(-sa*mu+ca*asofte(l))*
     &              conjg(sldlmx(k,l))+
     &              conjg(sldrmx(k,l)*rmass(kl(l))*
     &                    (-sa*mu+ca*asofte(l))*
     &                    conjg(sldlmx(j,l))))
            enddo
            gl(ksl(j),kh1,ksl(k)) = gl(kh1,ksl(j),ksl(k))
            gl(ksl(j),ksl(k),kh1) = gl(kh1,ksl(j),ksl(k))
            gr(kh1,ksl(j),ksl(k)) = gl(kh1,ksl(j),ksl(k))
            gr(ksl(j),kh1,ksl(k)) = gl(kh1,ksl(j),ksl(k))
            gr(ksl(j),ksl(k),kh1) = gl(kh1,ksl(j),ksl(k))
c h2-u~-u~
c     in dsvertx1
c h2-d~-d~
c     in dsvertx1
c h2-e~-e~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh2,ksl(j),ksl(k)) = g2o2cw*mass(kz)*sbpa*(
     &           (-1.d0-2.d0*echarg(ke)*sw2)*gmell(j,k)+
     &           (2.d0*echarg(ke)*sw2)*gmerr(j,k))
            do l=1,3
               gl(kh2,ksl(j),ksl(k)) = gl(kh2,ksl(j),ksl(k))
     &              +g2weak*sa/cb/mass(kw)*(
     &              sldrmx(j,l)*rmass(kl(l))**2*conjg(sldrmx(k,l))+
     &              sldlmx(j,l)*rmass(kl(l))**2*conjg(sldlmx(k,l)))
     &              -g2weak/cb/2.d0/mass(kw)*(
     &              sldrmx(j,l)*rmass(kl(l))*(-ca*mu-sa*asofte(l))*
     &              conjg(sldlmx(k,l))+
     &              conjg(sldrmx(k,l)*rmass(kl(l))*
     &                    (-ca*mu-sa*asofte(l))*
     &                    conjg(sldlmx(j,l))))
            enddo
            gl(ksl(j),kh2,ksl(k)) = gl(kh2,ksl(j),ksl(k))
            gl(ksl(j),ksl(k),kh2) = gl(kh2,ksl(j),ksl(k))
            gr(kh2,ksl(j),ksl(k)) = gl(kh2,ksl(j),ksl(k))
            gr(ksl(j),kh2,ksl(k)) = gl(kh2,ksl(j),ksl(k))
            gr(ksl(j),ksl(k),kh2) = gl(kh2,ksl(j),ksl(k))
c h3-u~-u~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh3,ksqu(j),ksqu(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh3,ksqu(j),ksqu(k)) = gl(kh3,ksqu(j),ksqu(k)) +
     &              squrmx(j,l)*rmass(kqu(l))*(-mu-asoftu(l)*cb/sb)*
     &              conjg(squlmx(k,l)) - conjg(
     &              squrmx(k,l)*rmass(kqu(l))*(-mu-asoftu(l)*cb/sb)*
     &              conjg(squlmx(j,l)))
            enddo
            gl(kh3,ksqu(j),ksqu(k)) = gl(kh3,ksqu(j),ksqu(k))*
     &           dcmplx(0.d0,g2weak/2.d0/mass(kw))
            gl(ksqu(j),kh3,ksqu(k)) = gl(kh3,ksqu(j),ksqu(k))
            gl(ksqu(j),ksqu(k),kh3) = gl(kh3,ksqu(j),ksqu(k))
            gr(kh3,ksqu(j),ksqu(k)) = gl(kh3,ksqu(j),ksqu(k))
            gr(ksqu(j),kh3,ksqu(k)) = gl(kh3,ksqu(j),ksqu(k))
            gr(ksqu(j),ksqu(k),kh3) = gl(kh3,ksqu(j),ksqu(k))
c h3-d~-d~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh3,ksqd(j),ksqd(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh3,ksqd(j),ksqd(k)) = gl(kh3,ksqd(j),ksqd(k)) +
     &              sqdrmx(j,l)*rmass(kqd(l))*(-mu-asoftd(l)*sb/cb)*
     &              conjg(sqdlmx(k,l)) - conjg(
     &              sqdrmx(k,l)*rmass(kqd(l))*(-mu-asoftd(l)*sb/cb)*
     &              conjg(sqdlmx(j,l)))
            enddo
            gl(kh3,ksqd(j),ksqd(k)) = gl(kh3,ksqd(j),ksqd(k))*
     &           dcmplx(0.d0,g2weak/2.d0/mass(kw))
            gl(ksqd(j),kh3,ksqd(k)) = gl(kh3,ksqd(j),ksqd(k))
            gl(ksqd(j),ksqd(k),kh3) = gl(kh3,ksqd(j),ksqd(k))
            gr(kh3,ksqd(j),ksqd(k)) = gl(kh3,ksqd(j),ksqd(k))
            gr(ksqd(j),kh3,ksqd(k)) = gl(kh3,ksqd(j),ksqd(k))
            gr(ksqd(j),ksqd(k),kh3) = gl(kh3,ksqd(j),ksqd(k))
c h3-e~-e~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(kh3,ksl(j),ksl(k)) = dcmplx(0.d0,0.d0)
            do l=1,3
               gl(kh3,ksl(j),ksl(k)) = gl(kh3,ksl(j),ksl(k)) +
     &              sldrmx(j,l)*rmass(kl(l))*(-mu-asofte(l)*sb/cb)*
     &              conjg(sldlmx(k,l)) - conjg(
     &              sldrmx(k,l)*rmass(kl(l))*(-mu-asofte(l)*sb/cb)*
     &              conjg(sldlmx(j,l)))
            enddo
            gl(kh3,ksl(j),ksl(k)) = gl(kh3,ksl(j),ksl(k))*
     &           dcmplx(0.d0,g2weak/2.d0/mass(kw))
            gl(ksl(j),kh3,ksl(k)) = gl(kh3,ksl(j),ksl(k))
            gl(ksl(j),ksl(k),kh3) = gl(kh3,ksl(j),ksl(k))
            gr(kh3,ksl(j),ksl(k)) = gl(kh3,ksl(j),ksl(k))
            gr(ksl(j),kh3,ksl(k)) = gl(kh3,ksl(j),ksl(k))
            gr(ksl(j),ksl(k),kh3) = gl(kh3,ksl(j),ksl(k))
c hc-u~-d~
            gl(khc,ksqu(j),ksqd(k)) = dcmplx(0.d0,0.d0)
            do a=1,3
               do b=1,3
                  auxc1=dcmplx(0.d0,0.d0)
                  auxc2=dcmplx(0.d0,0.d0)
                  do c=1,3
                     auxc1=auxc1+ckm(a,c)*conjg(ad(b,c))*rmass(kqd(b))
                     auxc2=auxc2+rmass(kqu(a))*au(a,c)*ckm(c,b)
                  enddo
                  gl(khc,ksqu(j),ksqd(k)) = gl(khc,ksqu(j),ksqd(k)) +
     &                 squlmx(j,a)*
     &                 (rmass(kqu(a))**2*ckm(a,b)*cb/sb
     &                 +ckm(a,b)*rmass(kqd(b))**2*sb/cb
     &                 -ckm(a,b)*mass(kw)**2*2.d0*sb*cb
     &                 )*conjg(sqdlmx(k,b)) +
     &                 squlmx(j,a)*
     &                 (auxc1*sb/cb+mu*ckm(a,b)*rmass(kqd(b))
     &                 )*conjg(sqdrmx(k,b)) +
     &                 squrmx(j,a)*
     &                 (auxc2*cb/sb+mu*rmass(kqu(a))*ckm(a,b)
     &                 )*conjg(sqdlmx(k,b)) +
     &                 squrmx(j,a)*
     &                 (rmass(kqu(a))*ckm(a,b)*rmass(kqd(b))/(sb*cb)
     &                 )*conjg(sqdrmx(k,b))
               enddo
            enddo
            gl(khc,ksqu(j),ksqd(k)) = gl(khc,ksqu(j),ksqd(k))*
     &                 g2weak/sqrt(2.d0)/mass(kw)
            gl(ksqu(j),khc,ksqd(k)) = gl(khc,ksqu(j),ksqd(k))
            gl(ksqu(j),ksqd(k),khc) = gl(khc,ksqu(j),ksqd(k))
            gr(khc,ksqu(j),ksqd(k)) = gl(khc,ksqu(j),ksqd(k))
            gr(ksqu(j),khc,ksqd(k)) = gl(khc,ksqu(j),ksqd(k))
            gr(ksqu(j),ksqd(k),khc) = gl(khc,ksqu(j),ksqd(k))
         enddo
      enddo
      do j=1,6
         do k=1,6
c hc-d~-u~
            gl(khc,ksqd(j),ksqu(k)) = conjg(gl(khc,ksqu(k),ksqd(j)))
            gl(ksqd(j),khc,ksqu(k)) = gl(khc,ksqd(j),ksqu(k))
            gl(ksqd(j),ksqu(k),khc) = gl(khc,ksqd(j),ksqu(k))
            gr(khc,ksqd(j),ksqu(k)) = gl(khc,ksqd(j),ksqu(k))
            gr(ksqd(j),khc,ksqu(k)) = gl(khc,ksqd(j),ksqu(k))
            gr(ksqd(j),ksqu(k),khc) = gl(khc,ksqd(j),ksqu(k))
         enddo
      enddo
      do j=1,3
         do k=1,3
c h1-nu~-nu~
            gl(kh1,ksnu(j),ksnu(k)) = -g2o2cw*mass(kz)*cbpa*gmnll(j,k) 
            gl(ksnu(j),kh1,ksnu(k)) = gl(kh1,ksnu(j),ksnu(k))
            gl(ksnu(j),ksnu(k),kh1) = gl(kh1,ksnu(j),ksnu(k))
            gr(kh1,ksnu(j),ksnu(k)) = gl(kh1,ksnu(j),ksnu(k))
            gr(ksnu(j),kh1,ksnu(k)) = gl(kh1,ksnu(j),ksnu(k))
            gr(ksnu(j),ksnu(k),kh1) = gl(kh1,ksnu(j),ksnu(k))
c h2-nu~-nu~
            gl(kh2,ksnu(j),ksnu(k)) = g2o2cw*mass(kz)*sbpa*gmnll(j,k) 
            gl(ksnu(j),kh2,ksnu(k)) = gl(kh2,ksnu(j),ksnu(k))
            gl(ksnu(j),ksnu(k),kh2) = gl(kh2,ksnu(j),ksnu(k))
            gr(kh2,ksnu(j),ksnu(k)) = gl(kh2,ksnu(j),ksnu(k))
            gr(ksnu(j),kh2,ksnu(k)) = gl(kh2,ksnu(j),ksnu(k))
            gr(ksnu(j),ksnu(k),kh2) = gl(kh2,ksnu(j),ksnu(k))
c h3-nu~-nu~
            gl(kh3,ksnu(j),ksnu(k)) = dcmplx(0.d0,0.d0)
            gl(ksnu(j),kh3,ksnu(k)) = gl(kh3,ksnu(j),ksnu(k))
            gl(ksnu(j),ksnu(k),kh3) = gl(kh3,ksnu(j),ksnu(k))
            gr(kh3,ksnu(j),ksnu(k)) = gl(kh3,ksnu(j),ksnu(k))
            gr(ksnu(j),kh3,ksnu(k)) = gl(kh3,ksnu(j),ksnu(k))
            gr(ksnu(j),ksnu(k),kh3) = gl(kh3,ksnu(j),ksnu(k))
         enddo
      enddo
      do j=1,3
         do k=1,6
c hc-nu~-e~
c...Sign of A-term changed 020903 (Mia Schelke)
            gl(khc,ksnu(j),ksl(k)) = dcmplx(0.d0,0.d0)
            do a=1,3
               gl(khc,ksnu(j),ksl(k)) = gl(khc,ksnu(j),ksl(k)) +
     &              slulmx(j,a)*mu*rmass(kl(a))*conjg(sldrmx(k,a))+
     &              slulmx(j,a)*(
     &              -mass(kw)**2*2.d0*sb*cb+
     &              rmass(kl(a))**2*sb/cb
     &              )*conjg(sldlmx(k,a))
               do b=1,3
                  gl(khc,ksnu(j),ksl(k)) = gl(khc,ksnu(j),ksl(k)) +
     &                 slulmx(j,a)*(
     &                  rmass(kl(a))*conjg(ae(b,a))*sb/cb
     &                 )*conjg(sldrmx(k,b))
               enddo
            enddo
            gl(khc,ksnu(j),ksl(k)) = gl(khc,ksnu(j),ksl(k))*
     &           g2weak/sqrt(2.d0)/mass(kw)
            gl(ksnu(j),khc,ksl(k)) = gl(khc,ksnu(j),ksl(k))
            gl(ksnu(j),ksl(k),khc) = gl(khc,ksnu(j),ksl(k))
            gr(khc,ksnu(j),ksl(k)) = gl(khc,ksnu(j),ksl(k))
            gr(ksnu(j),khc,ksl(k)) = gl(khc,ksnu(j),ksl(k))
            gr(ksnu(j),ksl(k),khc) = gl(khc,ksnu(j),ksl(k))
         enddo
      enddo
      do j=1,6
         do k=1,3
c hc-e~-nu~
            gl(khc,ksl(j),ksnu(k)) = conjg(gl(khc,ksnu(k),ksl(j)))
            gl(ksl(j),khc,ksnu(k)) = gl(khc,ksl(j),ksnu(k))
            gl(ksl(j),ksnu(k),khc) = gl(khc,ksl(j),ksnu(k))
            gr(khc,ksl(j),ksnu(k)) = gl(khc,ksl(j),ksnu(k))
            gr(ksl(j),khc,ksnu(k)) = gl(khc,ksl(j),ksnu(k))
            gr(ksl(j),ksnu(k),khc) = gl(khc,ksl(j),ksnu(k))
         enddo
      enddo
            
c-----------------------------------------beginning of quartic couplings
c
c     for each set of four particles (p1,p2,p3,p4), three vertices must
c     be set: (p1,p2,p3,p4), (p2,p1,p3,p4), and (p1,p3,p2,p4), obtained
c     via (1<->2) and (2<->3), unless the permutation can be obtained
c     from a previous permutation via one of the following exchanges:
c     (1<->3), (2<->4), or (1<->2,3<->4). The vertices for the permuted
c     particles are obtained according to the following rules: let g be
c     the (p1,p2,p3,p4) vertex, then
c
c                               g(1<->2)       g(2<->3)
c     four neutral particles       g              g
c     two neutral: 1 and 2         g              0
c     two neutral: 1 and 3         0              g^*
c     two neutral: 2 and 3         0              g
c     two neutral: 3 and 4         g^*            0
c     one neutral: 1
c     one neutral: 2
c     one neutral: 3
c     one neutral: 4
c     none neutral              independent     independent
c
c     other combinations of neutral and charged are forbidden
c     

c--------------------------------------------------------W-W-higgs-higgs
c     Higgs Hunter's Guide p. 364
      aux=gsq/2.d0
      call dsg4set34(kw,kw,kh1,kh1,aux,0.d0)
      call dsg4set34(kw,kw,kh2,kh2,aux,0.d0)
      call dsg4set34(kw,kw,kh3,kh3,aux,0.d0)
      call dsg4set(kw,kw,khc,khc,aux,0.d0)
      call dsg4set(kw,khc,kw,khc,aux,0.d0)

c--------------------------------------------------------Z-Z-higgs-higgs
c     Higgs Hunter's Guide p. 364
      aux=gsq/2.d0/costhw**2
      call dsg4set1234(kz,kz,kh1,kh1,aux,0.d0)
      call dsg4set1234(kz,kz,kh2,kh2,aux,0.d0)
      call dsg4set1234(kz,kz,kh3,kh3,aux,0.d0)
      aux=gsq/2.d0/costhw**2*(costhw**2-sinthw**2)**2
      call dsg4set12(kz,kz,khc,khc,aux,0.d0)

c--------------------------------------------------------Z-W-higgs-higgs
c     Higgs Hunter's Guide p. 365
      aux=gsq*sinthw**2/2.d0/costhw*sbma
      call dsg4set12(kh1,kz,kw,khc,aux,0.d0)
      aux=-gsq*sinthw**2/2.d0/costhw*cbma
      call dsg4set12(kh2,kz,kw,khc,aux,0.d0)
      aux=gsq*sinthw**2/2.d0/costhw
      call dsg4set12(kh3,kz,kw,khc,0.d0,aux)

c------------------------------------------------higgs-higgs-higgs-higgs
c     Higgs Hunter's Guide p. 376
      call dsg4set(kgoldc,kgoldc,kgoldc,kgoldc,
     &     -gsq/2.d0/cw2*c2b**2,0.d0)
      call dsg4set(khc,khc,khc,khc,-gsq/2.d0/cw2*c2b**2,0.d0)
      call dsg4set(kgoldc,kgoldc,kgoldc,khc,gsq/2.d0/cw2*s2b*c2b,0.d0)
      call dsg4set(khc,khc,khc,kgoldc,-gsq/2.d0/cw2*s2b*c2b,0.d0)
      call dsg4set(khc,kgoldc,khc,kgoldc,-gsq/2.d0/cw2*s2b**2,0.d0)
      call dsg4set(khc,khc,kgoldc,kgoldc,
     &     gsq/4.d0/cw2*(c2b**2-s2b**2),0.d0)
      call dsg4set12(kh1,kh1,kgoldc,kgoldc,
     &     -gsq/4.d0*(1.d0-s2b*s2a+tw2*c2b*c2a),0.d0)
      call dsg4set12(kh1,kh1,khc,khc,
     &     -gsq/4.d0*(1.d0+s2b*s2a-tw2*c2b*c2a),0.d0)
      call dsg4set12(kh1,kh1,khc,kgoldc,
     &     gsq/4.d0*(c2b*s2a+tw2*s2b*c2a),0.d0)
      call dsg4set12(kh2,kh2,kgoldc,kgoldc,
     &     -gsq/4.d0*(1.d0+s2b*s2a-tw2*c2b*c2a),0.d0)
      call dsg4set12(kh2,kh2,khc,khc,
     &     -gsq/4.d0*(1.d0-s2b*s2a+tw2*c2b*c2a),0.d0)
      call dsg4set12(kh2,kh2,khc,kgoldc,
     &     -gsq/4.d0*(c2b*s2a+tw2*s2b*c2a),0.d0)
      call dsg4set12(kh1,kh2,kgoldc,kgoldc,
     &     gsq/4.d0*(s2b*c2a+tw2*c2b*s2a),0.d0)
      call dsg4set12(kh1,kh2,khc,khc,
     &     -gsq/4.d0*(s2b*c2a+tw2*c2b*s2a),0.d0)
      call dsg4set12(kh1,kh2,khc,kgoldc,
     &     gsq/4.d0*(c2b*c2a-tw2*s2b*s2a),0.d0)
      call dsg4set12(kh3,kh3,kgoldc,kgoldc,
     &     -gsq/4.d0*(1.d0+s2b**2-tw2*c2b**2),0.d0)
      call dsg4set12(kh3,kh3,khc,khc,
     &     -gsq/4.d0/cw2*c2b**2,0.d0)
      call dsg4set12(kh3,kh3,khc,kgoldc,
     &     -gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set12(kgold0,kgold0,kgoldc,kgoldc,
     &     -gsq/4.d0/cw2*c2b**2,0.d0)
      call dsg4set12(kgold0,kgold0,khc,khc,
     &     -gsq/4.d0*(1.d0+s2b**2-tw2*c2b**2),0.d0)
      call dsg4set12(kgold0,kgold0,khc,kgoldc,
     &     gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set12(kh3,kgold0,kgoldc,kgoldc,
     &     gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set12(kh3,kgold0,khc,khc,
     &     -gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set12(kh3,kgold0,khc,kgoldc,
     &     gsq/4.d0*(c2b**2-tw2*s2b**2),0.d0)
      call dsg4set12(kh1,kgold0,khc,kgoldc,
     &     0.d0,gsq/4.d0*sbma)
      call dsg4set12(kh2,kgold0,khc,kgoldc,
     &     0.d0,-gsq/4.d0*cbma)
      call dsg4set12(kh1,kh3,khc,kgoldc,
     &     0.d0,gsq/4.d0*cbma)
      call dsg4set12(kh2,kh3,khc,kgoldc,
     &     0.d0,gsq/4.d0*sbma)
      call dsg4set(kgold0,kgold0,kgold0,kgold0,
     &     -3.d0*gsq/4.d0/cw2*c2b**2,0.d0)
      call dsg4set(kh1,kh1,kh1,kh1,
     &     -3.d0*gsq/4.d0/cw2*c2a**2,0.d0)
      call dsg4set(kh2,kh2,kh2,kh2,
     &     -3.d0*gsq/4.d0/cw2*c2a**2,0.d0)
      call dsg4set(kh3,kh3,kh3,kh3,
     &     -3.d0*gsq/4.d0/cw2*c2b**2,0.d0)
      call dsg4set1234(kh2,kh1,kh1,kh1,
     &     3.d0*gsq/4.d0/cw2*s2a*c2a,0.d0)
      call dsg4set1234(kh1,kh2,kh2,kh2,
     &     -3.d0*gsq/4.d0/cw2*s2a*c2a,0.d0)
      call dsg4set1234(kgold0,kh3,kh3,kh3,
     &     -3.d0*gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set1234(kh3,kgold0,kgold0,kgold0,
     &     3.d0*gsq/4.d0/cw2*s2b*c2b,0.d0)
      call dsg4set1234(kh1,kh1,kh2,kh2,
     &     -gsq/4.d0/cw2*(3.d0*s2a**2-1.d0),0.d0)
      call dsg4set1234(kh3,kh3,kgold0,kgold0,
     &     -gsq/4.d0/cw2*(3.d0*s2b**2-1.d0),0.d0)
      call dsg4set1234(kh1,kh1,kgold0,kgold0,
     &     -gsq/4.d0/cw2*c2b*c2a,0.d0)
      call dsg4set1234(kh2,kh2,kgold0,kgold0,
     &     gsq/4.d0/cw2*c2b*c2a,0.d0)
      call dsg4set1234(kh1,kh2,kgold0,kgold0,
     &     gsq/4.d0/cw2*c2b*s2a,0.d0)
      call dsg4set1234(kh1,kh1,kh3,kgold0,
     &     gsq/4.d0/cw2*s2b*c2a,0.d0)
      call dsg4set1234(kh2,kh2,kh3,kgold0,
     &     -gsq/4.d0/cw2*s2b*c2a,0.d0)
      call dsg4set1234(kh1,kh2,kh3,kgold0,
     &     -gsq/4.d0/cw2*s2b*s2a,0.d0)
      call dsg4set1234(kh1,kh1,kh3,kh3,
     &     gsq/4.d0/cw2*c2b*c2a,0.d0)
      call dsg4set1234(kh2,kh2,kh3,kh3,
     &     -gsq/4.d0/cw2*c2b*c2a,0.d0)
      call dsg4set1234(kh1,kh2,kh3,kh3,
     &     -gsq/4.d0/cw2*c2b*s2a,0.d0)

c------------------------------------------sfermion-sfermion-gauge-gauge
c     Haber and Kane p. 227, Jungman et al p. 353
      do k=1,6
         do l=1,6
            ! gamma-gamma-sf-sf
            aux=2.d0*(e_charge*echarg(ku))**2
            if (k.eq.l) then
               call dsg4set(kgamma,kgamma,ksqu(k),ksqu(l),aux,0.d0)
            else
               call dsg4set(kgamma,kgamma,ksqu(k),ksqu(l),0.d0,0.d0)
            endif
            call dsg4set(kgamma,ksqu(k),kgamma,ksqu(l),0.d0,0.d0)
            aux=2.d0*(e_charge*echarg(kd))**2
            if (k.eq.l) then
               call dsg4set(kgamma,kgamma,ksqd(k),ksqd(l),aux,0.d0)
            else
               call dsg4set(kgamma,kgamma,ksqd(k),ksqd(l),0.d0,0.d0)
            endif
            call dsg4set(kgamma,ksqd(k),kgamma,ksqd(l),0.d0,0.d0)
            aux=2.d0*(e_charge*echarg(ke))**2
            if (k.eq.l) then
               call dsg4set(kgamma,kgamma,ksl(k),ksl(l),aux,0.d0)
            else
               call dsg4set(kgamma,kgamma,ksl(k),ksl(l),0.d0,0.d0)
            endif
            call dsg4set(kgamma,ksl(k),kgamma,ksl(l),0.d0,0.d0)
            ! gluon-gluon-sf-sf
            ! L = vrtx/2 (G_mu)^a (G^mu)^b q~^*_i C_{abij} q~_j
            ! with C_{abij} = delta_{ab} delta_{ij}/3 + d_{abc} T^c_{ij}
            aux=g3stro**2
            if (k.eq.l) then
               call dsg4set(kgluon,kgluon,ksqu(k),ksqu(l),aux,0.d0)
            else
               call dsg4set(kgluon,kgluon,ksqu(k),ksqu(l),0.d0,0.d0)
            endif
            call dsg4set(kgluon,ksqu(k),kgluon,ksqu(l),0.d0,0.d0)
            if (k.eq.l) then
               call dsg4set(kgluon,kgluon,ksqd(k),ksqd(l),aux,0.d0)
            else
               call dsg4set(kgluon,kgluon,ksqd(k),ksqd(l),0.d0,0.d0)
            endif
            call dsg4set(kgluon,ksqd(k),kgluon,ksqd(l),0.d0,0.d0)
            call dsg4set(kgluon,kgluon,ksl(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgluon,ksl(k),kgluon,ksl(l),0.d0,0.d0)
            ! gluon-gamma-sf-sf
            ! L = vrtx (G_mu)^a (A^mu) q~^*_i T^a_{ij} q~_j
            aux=2.d0*g3stro*e_charge*echarg(ku)
            if (k.eq.l) then
               call dsg4set(kgamma,kgluon,ksqu(k),ksqu(l),aux,0.d0)
               call dsg4set(kgluon,kgamma,ksqu(k),ksqu(l),aux,0.d0)
            else
               call dsg4set(kgamma,kgluon,ksqu(k),ksqu(l),0.d0,0.d0)
               call dsg4set(kgluon,kgamma,ksqu(k),ksqu(l),0.d0,0.d0)
            endif
            call dsg4set(kgamma,ksqu(k),kgluon,ksqu(l),0.d0,0.d0)
            aux=2.d0*g3stro*e_charge*echarg(kd)
            if (k.eq.l) then
               call dsg4set(kgamma,kgluon,ksqd(k),ksqd(l),aux,0.d0)
               call dsg4set(kgluon,kgamma,ksqd(k),ksqd(l),aux,0.d0)
            else
               call dsg4set(kgamma,kgluon,ksqd(k),ksqd(l),0.d0,0.d0)
               call dsg4set(kgluon,kgamma,ksqd(k),ksqd(l),0.d0,0.d0)
            endif
            call dsg4set(kgamma,ksqd(k),kgluon,ksqd(l),0.d0,0.d0)
            aux=0.d0
            call dsg4set(kgamma,kgluon,ksl(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgluon,kgamma,ksl(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgamma,ksl(k),kgluon,ksl(l),0.d0,0.d0)
            ! gluon-Z-sf-sf
            ! L = vrtx (G_mu)^a (Z^mu) q~^*_i T^a_{ij} q~_j
            auxl1=2.d0*g3stro*g2weak/costhw*(0.5d0-sw2*echarg(ku))
            auxr1=2.d0*g3stro*g2weak/costhw*(-sw2*echarg(ku))
            auxc=auxl1*gmull(k,l)+auxr1*gmurr(k,l)
            call dsg4setc(kgluon,kz,ksqu(k),ksqu(l),auxc)
            call dsg4setc(kz,kgluon,ksqu(k),ksqu(l),auxc)
            call dsg4set(kgluon,ksqu(k),kz,ksqu(l),0.d0,0.d0)
            auxl1=2.d0*g3stro*g2weak/costhw*(-0.5d0-sw2*echarg(kd))
            auxr1=2.d0*g3stro*g2weak/costhw*(-sw2*echarg(kd))
            auxc=auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)
            call dsg4setc(kgluon,kz,ksqd(k),ksqd(l),auxc)
            call dsg4setc(kz,kgluon,ksqd(k),ksqd(l),auxc)
            call dsg4set(kgluon,ksqd(k),kz,ksqd(l),0.d0,0.d0)
            aux=0.d0
            call dsg4set(kgluon,kz,ksl(k),ksl(l),0.d0,0.d0)
            call dsg4set(kz,kgluon,ksl(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgluon,ksl(k),kz,ksl(l),0.d0,0.d0)
            ! gamma-Z-sf-sf
            auxl1=2.d0*g2weak*e_charge/costhw*
     &           (0.5d0-sw2*echarg(ku))*echarg(ku)
            auxr1=2.d0*g2weak*e_charge/costhw*
     &           (-sw2*echarg(ku))*echarg(ku)
            auxc=auxl1*gmull(k,l)+auxr1*gmurr(k,l)
            call dsg4setc(kgamma,kz,ksqu(k),ksqu(l),auxc)
            call dsg4setc(kz,kgamma,ksqu(k),ksqu(l),auxc)
            call dsg4set(kgamma,ksqu(k),kz,ksqu(l),0.d0,0.d0)
            auxl1=2.d0*g2weak*e_charge/costhw*
     &           (-0.5d0-sw2*echarg(kd))*echarg(kd)
            auxr1=2.d0*g2weak*e_charge/costhw*
     &           (-sw2*echarg(kd))*echarg(kd)
            auxc=auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)
            call dsg4setc(kgamma,kz,ksqd(k),ksqd(l),auxc)
            call dsg4setc(kz,kgamma,ksqd(k),ksqd(l),auxc)
            call dsg4set(kgamma,ksqd(k),kz,ksqd(l),0.d0,0.d0)
            auxl1=2.d0*g2weak*e_charge/costhw*
     &           (-0.5d0-sw2*echarg(ke))*echarg(ke)
            auxr1=2.d0*g2weak*e_charge/costhw*
     &           (-sw2*echarg(ke))*echarg(ke)
            auxc=auxl1*gmell(k,l)+auxr1*gmerr(k,l)
            call dsg4setc(kgamma,kz,ksl(k),ksl(l),auxc)
            call dsg4setc(kz,kgamma,ksl(k),ksl(l),auxc)
            call dsg4set(kgamma,ksl(k),kz,ksl(l),0.d0,0.d0)
            ! Z-Z-sf-sf
            auxl1=2.d0*(g2weak*(0.5d0-sw2*echarg(ku))/costhw)**2
            auxr1=2.d0*(g2weak*(-sw2*echarg(ku))/costhw)**2
            auxc=auxl1*gmull(k,l)+auxr1*gmurr(k,l)
            call dsg4setc12(kz,kz,ksqu(k),ksqu(l),auxc)
            call dsg4set(kz,ksqu(k),kz,ksqu(l),0.d0,0.d0)
            auxl1=2.d0*(g2weak*(-0.5d0-sw2*echarg(kd))/costhw)**2
            auxr1=2.d0*(g2weak*(-sw2*echarg(kd))/costhw)**2
            auxc=auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)
            call dsg4setc12(kz,kz,ksqd(k),ksqd(l),auxc)
            call dsg4set(kz,ksqd(k),kz,ksqd(l),0.d0,0.d0)
            auxl1=2.d0*(g2weak*(-0.5d0-sw2*echarg(ke))/costhw)**2
            auxr1=2.d0*(g2weak*(-sw2*echarg(ke))/costhw)**2
            auxc=auxl1*gmell(k,l)+auxr1*gmerr(k,l)
            call dsg4setc12(kz,kz,ksl(k),ksl(l),auxc)
            call dsg4set(kz,ksl(k),kz,ksl(l),0.d0,0.d0)
            ! W-W-sf-sf
            aux=gsq/2.d0
            call dsg4setc(kw,kw,ksqu(k),ksqu(l),aux*gmull(k,l))
            call dsg4set(kw,ksqu(k),kw,ksqu(l),0.d0,0.d0)
            call dsg4setc(kw,kw,ksqd(k),ksqd(l),aux*gmdll(k,l))
            call dsg4set(kw,ksqd(k),kw,ksqd(l),0.d0,0.d0)
            call dsg4setc(kw,kw,ksl(k),ksl(l),aux*gmell(k,l))
            call dsg4set(kw,ksl(k),kw,ksl(l),0.d0,0.d0)
            ! gamma-W-su-sd
            auxc=dcmplx(0.d0,0.d0)
            do a=1,3
              do b=1,3
                auxc=auxc+squlmx(k,a)*ckm(a,b)*conjg(sqdlmx(l,b))
              enddo
            enddo
            aux=g2weak/sqrt(2.d0)*(-1.d0+2.d0*echarg(ku))*e_charge
            call dsg4setc(kgamma,kw,ksqu(k),ksqd(l),aux*auxc)
            call dsg4set(kw,kgamma,ksqu(k),ksqd(l),0.d0,0.d0)
            call dsg4set(kgamma,ksqu(k),kw,ksqd(l),0.d0,0.d0)
            call dsg4setc(kw,kgamma,ksqd(l),ksqu(k),aux*conjg(auxc))
            call dsg4set(kgamma,kw,ksqd(l),ksqu(k),0.d0,0.d0)
            call dsg4set(kw,ksqd(l),kgamma,ksqu(k),0.d0,0.d0)
            ! Z-W-su-sd
            aux=-gsq/sqrt(2.d0)*(-1.d0+2.d0*echarg(ku))*
     &           sw2/costhw
            call dsg4setc(kz,kw,ksqu(k),ksqd(l),aux*auxc)
            call dsg4set(kw,kz,ksqu(k),ksqd(l),0.d0,0.d0)
            call dsg4set(kz,ksqu(k),kw,ksqd(l),0.d0,0.d0)
            call dsg4setc(kw,kz,ksqd(l),ksqu(k),aux*conjg(auxc))
            call dsg4set(kz,kw,ksqd(l),ksqu(k),0.d0,0.d0)
            call dsg4set(kw,ksqd(l),kz,ksqu(k),0.d0,0.d0)
            ! gluon-W-su-sd
            ! L = vrtx (A^mu)_a (W^+)_mu u~^*_i T^a_{ij} d~_j
            aux=sqrt(2.d0)*g2weak*g3stro
            call dsg4setc(kgluon,kw,ksqu(k),ksqd(l),aux*auxc)
            call dsg4set(kw,kgluon,ksqu(k),ksqd(l),0.d0,0.d0)
            call dsg4set(kgluon,ksqu(k),kw,ksqd(l),0.d0,0.d0)
            call dsg4setc(kw,kgluon,ksqd(l),ksqu(k),aux*conjg(auxc))
            call dsg4set(kgluon,kw,ksqd(l),ksqu(k),0.d0,0.d0)
            call dsg4set(kw,ksqd(l),kgluon,ksqu(k),0.d0,0.d0)
         enddo
      enddo
      do k=1,3
         do l=1,6
            ! gamma-W-snu-se
            aux=g2weak/sqrt(2.d0)*(-1.d0)*e_charge
            call dsg4setc(kgamma,kw,ksnu(k),ksl(l),aux*gmnell(k,l))
            call dsg4set(kw,kgamma,ksnu(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgamma,ksnu(k),kw,ksl(l),0.d0,0.d0)
            call dsg4setc(kw,kgamma,ksl(l),ksnu(k),
     &           aux*conjg(gmnell(k,l)))
            call dsg4set(kgamma,kw,ksl(l),ksnu(k),0.d0,0.d0)
            call dsg4set(kw,ksl(l),kgamma,ksnu(k),0.d0,0.d0)
            ! gluon-W-snu-se
            call dsg4set(kgluon,kw,ksnu(k),ksl(l),0.d0,0.d0)
            call dsg4set(kw,kgluon,ksnu(k),ksl(l),0.d0,0.d0)
            call dsg4set(kgluon,ksnu(k),kw,ksl(l),0.d0,0.d0)
            call dsg4set(kw,kgluon,ksl(l),ksnu(k),0.d0,0.d0)
            call dsg4set(kgluon,kw,ksl(l),ksnu(k),0.d0,0.d0)
            call dsg4set(kw,ksl(l),kgluon,ksnu(k),0.d0,0.d0)
            ! Z-W-snu-se
            aux=-gsq/sqrt(2.d0)*(-1.d0)*sw2/costhw
            call dsg4setc(kz,kw,ksnu(k),ksl(l),aux*gmnell(k,l))
            call dsg4set(kw,kz,ksnu(k),ksl(l),0.d0,0.d0)
            call dsg4set(kz,ksnu(k),kw,ksl(l),0.d0,0.d0)
            call dsg4setc(kw,kz,ksl(l),ksnu(k),
     &           aux*conjg(gmnell(k,l)))
            call dsg4set(kz,kw,ksl(l),ksnu(k),0.d0,0.d0)
            call dsg4set(kw,ksl(l),kz,ksnu(k),0.d0,0.d0)
         enddo
      enddo
      do k=1,3
         do l=1,3
            ! W-W-snu-snu
            aux=gsq/2.d0
            call dsg4setc(kw,kw,ksnu(k),ksnu(l),aux*gmnll(k,l))
            call dsg4set(kw,ksnu(k),kw,ksnu(l),0.d0,0.d0)
            ! Z-Z-snu-snu, factor of two fix, Mia Schelke (020903)
            aux=2.d0*(g2weak*(0.5d0)/costhw)**2  ! Corr. 020903
            call dsg4setc(kz,kz,ksnu(k),ksnu(l),aux*gmnll(k,l))
            call dsg4set(kz,ksnu(k),kz,ksnu(l),0.d0,0.d0)
         enddo
      enddo

      
c------------------------------------------sfermion-sfermion-higgs-higgs
c     Higgs Hunter's Guide p. 402
      hj(1)=kh1
      hj(2)=kh2
      hj(3)=kh3
      hj(4)=kgold0
      cc(1)=-c2a
      cc(2)=c2a
      cc(3)=c2b
      cc(4)=-c2b
      dd(1,1)=sa2/sb2
      dd(2,1)=ca2/sb2
      dd(3,1)=cb2/sb2
      dd(4,1)=1.d0
      dd(1,2)=ca2/cb2
      dd(2,2)=sa2/cb2
      dd(3,2)=sb2/cb2
      dd(4,2)=1.d0
      rr(1)=dcmplx(sbpa,0.d0)
      rr(2)=dcmplx(cbpa,0.d0)
      rr(3)=dcmplx(0.d0,c2b)
      rr(4)=dcmplx(0.d0,s2b)
      ss(1)=dcmplx(sa*cb/sb2,0.d0)
      ss(2)=dcmplx(ca*cb/sb2,0.d0)
      ss(3)=dcmplx(0.d0,cb2/sb2)
      ss(4)=dcmplx(0.d0,cb/sb)
      tt(1)=dcmplx(ca*sb/cb2,0.d0)
      tt(2)=dcmplx(-sa*sb/cb2,0.d0)
      tt(3)=dcmplx(0.d0,-sb2/cb2)
      tt(4)=dcmplx(0.d0,sb/cb)
      vv(1)=dcmplx(cbma,0.d0)
      vv(2)=dcmplx(sbma,0.d0)
      vv(3)=dcmplx(0.d0,0.d0)
      vv(4)=dcmplx(0.d0,1.d0)
      rp(1)=dcmplx(-cbpa,0.d0)
      rp(2)=dcmplx(sbpa,0.d0)
      rp(3)=dcmplx(0.d0,s2b)
      rp(4)=dcmplx(0.d0,-c2b)
      sp(1)=dcmplx(sa/sb,0.d0)
      sp(2)=dcmplx(ca/sb,0.d0)
      sp(3)=dcmplx(0.d0,cb/sb)
      sp(4)=dcmplx(0.d0,1.d0)
      tp(1)=dcmplx(-ca/cb,0.d0)
      tp(2)=dcmplx(sa/cb,0.d0)
      tp(3)=dcmplx(0.d0,sb/cb)
      tp(4)=dcmplx(0.d0,-1.d0)
      vp(1)=dcmplx(sbma,0.d0)
      vp(2)=dcmplx(-cbma,0.d0)
      vp(3)=dcmplx(0.d0,-1.d0)
      vp(4)=dcmplx(0.d0,0.d0)
      do j=1,4 ! Hj=H1,H2,H3,G0
         ! Hj-Hj-u~-u~
         auxl1=gsq/2.d0*(0.5d0-echarg(ku)*sw2)/cw2*cc(j)
         auxl2=-gsq/2.d0*dd(j,1)/mass(kw)**2
         auxr1=gsq/2.d0*(echarg(ku)*sw2)/cw2*cc(j)
         auxr2=-gsq/2.d0*dd(j,1)/mass(kw)**2
         do k=1,6
            do l=1,6
               call dsg4setc12(hj(j),hj(j),ksqu(k),ksqu(l),
     &              auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &              auxl2*gm2ull(k,l)+auxr2*gm2urr(k,l))
            enddo
         enddo
         ! Hj-Hj-d~-d~
         auxl1=gsq/2.d0*(-0.5d0-echarg(kd)*sw2)/cw2*cc(j)
         auxl2=-gsq/2.d0*dd(j,2)/mass(kw)**2
         auxr1=gsq/2.d0*(echarg(kd)*sw2)/cw2*cc(j)
         auxr2=-gsq/2.d0*dd(j,2)/mass(kw)**2
         do k=1,6
            do l=1,6
               call dsg4setc12(hj(j),hj(j),ksqd(k),ksqd(l),
     &              auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &              auxl2*gm2dll(k,l)+auxr2*gm2drr(k,l))
            enddo
         enddo
         ! Hj-Hj-nu~-nu~
         auxl1=gsq/2.d0*(0.5d0)/cw2*cc(j)
         do k=1,3
            do l=1,3
               call dsg4setc12(hj(j),hj(j),ksnu(k),ksnu(l),
     &              auxl1*gmnll(k,l))
            enddo
         enddo
         ! Hj-Hj-e~-e~
         auxl1=gsq/2.d0*(-0.5d0-echarg(ke)*sw2)/cw2*cc(j)
         auxl2=-gsq/2.d0*dd(j,2)/mass(kw)**2
         auxr1=gsq/2.d0*(echarg(ke)*sw2)/cw2*cc(j)
         auxr2=-gsq/2.d0*dd(j,2)/mass(kw)**2
         do k=1,6
            do l=1,6
               call dsg4setc12(hj(j),hj(j),ksl(k),ksl(l),
     &              auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &              auxl2*gm2ell(k,l)+auxr2*gm2err(k,l))
            enddo
         enddo
      enddo
      ! H1-H2-u~-u~
      auxl1=gsq*s2a/2.d0*(
     &     (0.5d0-echarg(ku)*sw2)/cw2)
      auxl2=-gsq*s2a/4.d0/sb2/mass(kw)**2
      auxr1=gsq*s2a/2.d0*(echarg(ku)*sw2)/cw2
      auxr2=-gsq*s2a/4.d0/sb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kh1,kh2,ksqu(k),ksqu(l),
     &           auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &           auxl2*gm2ull(k,l)+auxr2*gm2urr(k,l))
         enddo
      enddo
      ! H1-H2-d~-d~
      auxl1=gsq*s2a/2.d0*(
     &     (-0.5d0-echarg(kd)*sw2)/cw2)
      auxl2=+gsq*s2a/4.d0/cb2/mass(kw)**2
      auxr1=gsq*s2a/2.d0*(echarg(kd)*sw2)/cw2
      auxr2=+gsq*s2a/4.d0/cb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kh1,kh2,ksqd(k),ksqd(l),
     &           auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &           auxl2*gm2dll(k,l)+auxr2*gm2drr(k,l))
         enddo
      enddo
      ! H1-H2-nu~-nu~
      auxl1=gsq*s2a/2.d0*((0.5d0)/cw2)
      do k=1,3
         do l=1,3
            call dsg4setc12(kh1,kh2,ksnu(k),ksnu(l),
     &           auxl1*gmnll(k,l))
         enddo
      enddo
      ! H1-H2-e~-e~
      auxl1=gsq*s2a/2.d0*(
     &     (-0.5d0-echarg(ke)*sw2)/cw2)
      auxl2=+gsq*s2a/4.d0/cb2/mass(kw)**2
      auxr1=gsq*s2a/2.d0*(echarg(ke)*sw2)/cw2
      auxr2=+gsq*s2a/4.d0/cb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kh1,kh2,ksl(k),ksl(l),
     &           auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &           auxl2*gm2ell(k,l)+auxr2*gm2err(k,l))
         enddo
      enddo
      ! G0-H3-u~-u~
      auxl1=gsq*s2b/2.d0*(0.5d0-echarg(ku)*sw2)/cw2
      auxl2=-gsq*s2b/4.d0/sb2/mass(kw)**2
      auxr1=gsq*s2b/2.d0*(echarg(ku)*sw2)/cw2
      auxr2=-gsq*s2b/4.d0/sb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kgold0,kh3,ksqu(k),ksqu(l),
     &           auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &           auxl2*gm2ull(k,l)+auxr2*gm2urr(k,l))
         enddo
      enddo
      ! G0-H3-d~-d~
      auxl1=gsq*s2b/2.d0*(-0.5d0-echarg(kd)*sw2)/cw2
      auxl2=+gsq*s2b/4.d0/cb2/mass(kw)**2
      auxr1=gsq*s2b/2.d0*(echarg(kd)*sw2)/cw2
      auxr2=+gsq*s2b/4.d0/cb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kgold0,kh3,ksqd(k),ksqd(l),
     &           auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &           auxl2*gm2dll(k,l)+auxr2*gm2drr(k,l))
         enddo
      enddo
      ! G0-H3-nu~-nu~
      auxl1=gsq*s2b/2.d0*(0.5d0)/cw2
      do k=1,3
         do l=1,3
            call dsg4setc12(kgold0,kh3,ksnu(k),ksnu(l),
     &           auxl1*gmnll(k,l))
         enddo
      enddo
      ! G0-H3-e~-e~
      auxl1=gsq*s2b/2.d0*(-0.5d0-echarg(ke)*sw2)/cw2
      auxl2=+gsq*s2b/4.d0/cb2/mass(kw)**2
      auxr1=gsq*s2b/2.d0*(echarg(ke)*sw2)/cw2
      auxr2=+gsq*s2b/4.d0/cb2/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc12(kgold0,kh3,ksl(k),ksl(l),
     &           auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &           auxl2*gm2ell(k,l)+auxr2*gm2err(k,l))
         enddo
      enddo
      ! Hc-Hc-u~-u~
      auxl1=gsq*c2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0-echarg(ku)*sw2)/cw2)
      auxl2=-gsq*(sb/cb)**2/2.d0/mass(kw)**2
      auxr1=gsq*c2b/2.d0*((+echarg(ku)*sw2)/cw2)
      auxr2=-gsq*(cb/sb)**2/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(khc,khc,ksqu(k),ksqu(l),
     &           auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &           auxl2*gmd2ull(k,l)+auxr2*gm2urr(k,l))
            call dsg4set(khc,ksqu(k),khc,ksqu(l),0.d0,0.d0)
         enddo
      enddo
      ! Hc-Hc-d~-d~
      auxl1=gsq*c2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(kd)*sw2)/cw2)
      auxl2=-gsq*(cb/sb)**2/2.d0/mass(kw)**2
      auxr1=gsq*c2b/2.d0*((+echarg(kd)*sw2)/cw2)
      auxr2=-gsq*(sb/cb)**2/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(khc,khc,ksqd(k),ksqd(l),
     &           auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &           auxl2*gmu2dll(k,l)+auxr2*gm2drr(k,l))
            call dsg4set(khc,ksqd(k),khc,ksqd(l),0.d0,0.d0)
         enddo
      enddo
      ! Hc-Hc-nu~-nu~
      auxl1=gsq*c2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0)/cw2)
      auxl2=-gsq*(sb/cb)**2/2.d0/mass(kw)**2
      do k=1,3
         do l=1,3
            call dsg4setc(khc,khc,ksnu(k),ksnu(l),
     &           auxl1*gmnll(k,l)+auxl2*gme2nll(k,l))
            call dsg4set(khc,ksnu(k),khc,ksnu(l),0.d0,0.d0)
         enddo
      enddo
      ! Hc-Hc-e~-e~
      auxl1=gsq*c2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(ke)*sw2)/cw2)
      auxr1=gsq*c2b/2.d0*((+echarg(ke)*sw2)/cw2)
      auxr2=-gsq*(sb/cb)**2/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(khc,khc,ksl(k),ksl(l),
     &           auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &           auxr2*gm2err(k,l))
            call dsg4set(khc,ksl(k),khc,ksl(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Gc-u~-u~
      auxl1=-gsq*c2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0-echarg(ku)*sw2)/cw2)
      auxl2=-gsq/2.d0/mass(kw)**2
      auxr1=-gsq*c2b/2.d0*((+echarg(ku)*sw2)/cw2)
      auxr2=-gsq/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(kgoldc,kgoldc,ksqu(k),ksqu(l),
     &           auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &           auxl2*gmd2ull(k,l)+auxr2*gm2urr(k,l))
            call dsg4set(kgoldc,ksqu(k),kgoldc,ksqu(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Gc-d~-d~
      auxl1=-gsq*c2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(kd)*sw2)/cw2)
      auxl2=-gsq/2.d0/mass(kw)**2
      auxr1=-gsq*c2b/2.d0*((+echarg(kd)*sw2)/cw2)
      auxr2=-gsq/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(kgoldc,kgoldc,ksqd(k),ksqd(l),
     &           auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &           auxl2*gmu2dll(k,l)+auxr2*gm2drr(k,l))
            call dsg4set(kgoldc,ksqd(k),kgoldc,ksqd(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Gc-nu~-nu~
      auxl1=-gsq*c2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0)/cw2)
      auxl2=-gsq/2.d0/mass(kw)**2
      do k=1,3
         do l=1,3
            call dsg4setc(kgoldc,kgoldc,ksnu(k),ksnu(l),
     &           auxl1*gmnll(k,l)+auxl2*gme2nll(k,l))
            call dsg4set(kgoldc,ksnu(k),kgoldc,ksnu(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Gc-e~-e~
      auxl1=-gsq*c2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(ke)*sw2)/cw2)
      auxr1=-gsq*c2b/2.d0*((+echarg(ke)*sw2)/cw2)
      auxr2=-gsq/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            call dsg4setc(kgoldc,kgoldc,ksl(k),ksl(l),
     &           auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &           auxr2*gm2err(k,l))
            call dsg4set(kgoldc,ksl(k),kgoldc,ksl(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Hc-u~-u~
      auxl1=gsq*s2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0-echarg(ku)*sw2)/cw2)
      auxl2=-gsq*(-sb/cb)/2.d0/mass(kw)**2
      auxr1=gsq*s2b/2.d0*((+echarg(ku)*sw2)/cw2)
      auxr2=-gsq*(cb/sb)/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            auxc=auxl1*gmull(k,l)+auxr1*gmurr(k,l)+
     &           auxl2*gmd2ull(k,l)+auxr2*gm2urr(k,l)
            call dsg4setc(kgoldc,khc,ksqu(k),ksqu(l),auxc)
            call dsg4setc(khc,kgoldc,ksqu(k),ksqu(l),auxc)
            call dsg4set(kgoldc,ksqu(k),khc,ksqu(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Hc-d~-d~
      auxl1=gsq*s2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(kd)*sw2)/cw2)
      auxl2=-gsq*(cb/sb)/2.d0/mass(kw)**2
      auxr1=gsq*s2b/2.d0*((+echarg(kd)*sw2)/cw2)
      auxr2=-gsq*(-sb/cb)/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            auxc=auxl1*gmdll(k,l)+auxr1*gmdrr(k,l)+
     &           auxl2*gmu2dll(k,l)+auxr2*gm2drr(k,l)
            call dsg4setc(kgoldc,khc,ksqd(k),ksqd(l),auxc)
            call dsg4setc(khc,kgoldc,ksqd(k),ksqd(l),auxc)
            call dsg4set(kgoldc,ksqd(k),khc,ksqd(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Hc-nu~-nu~
      auxl1=gsq*s2b/2.d0*(-2.d0*0.5d0+
     &     (0.5d0)/cw2)
      auxl2=-gsq*(-sb/cb)/2.d0/mass(kw)**2
      do k=1,3
         do l=1,3
            auxc=auxl1*gmnll(k,l)+auxl2*gme2nll(k,l)
            call dsg4setc(kgoldc,khc,ksnu(k),ksnu(l),auxc)
            call dsg4setc(khc,kgoldc,ksnu(k),ksnu(l),auxc)
            call dsg4set(kgoldc,ksnu(k),khc,ksnu(l),0.d0,0.d0)
         enddo
      enddo
      ! Gc-Hc-e~-e~
      auxl1=gsq*s2b/2.d0*(+2.d0*0.5d0+
     &     (-0.5d0-echarg(ke)*sw2)/cw2)
      auxr1=gsq*s2b/2.d0*((+echarg(ke)*sw2)/cw2)
      auxr2=-gsq*(-sb/cb)/2.d0/mass(kw)**2
      do k=1,6
         do l=1,6
            auxc=auxl1*gmell(k,l)+auxr1*gmerr(k,l)+
     &           auxr2*gm2err(k,l)
            call dsg4setc(kgoldc,khc,ksl(k),ksl(l),auxc)
            call dsg4setc(khc,kgoldc,ksl(k),ksl(l),auxc)
            call dsg4set(kgoldc,ksl(k),khc,ksl(l),0.d0,0.d0)
         enddo
      enddo
      do j=1,4 ! Hj=H1,H2,H3,G0
         ! Hj-Hc-u~-d~
         do k=1,6
            do l=1,6
               auxc=-conjg(rr(j))*gmudll(k,l)
               do m=1,3
                  do n=1,3
                     auxc=auxc
     &                    +conjg(ss(j))/mass(kw)**2*squlmx(k,m)*
     &                    rmass(kqu(m))**2*ckm(m,n)*conjg(sqdlmx(l,n))
     &                    +conjg(tt(j))/mass(kw)**2*squlmx(k,m)*
     &                    ckm(m,n)*rmass(kqd(n))**2*conjg(sqdlmx(l,n))
     &                    +2.d0*conjg(vv(j))/mass(kw)**2/s2b*
     &                    squrmx(k,m)*rmass(kqu(m))*ckm(m,n)*
     &                    rmass(kqd(n))*conjg(sqdrmx(l,n))
                  enddo
               enddo
               auxc=auxc*gsq/2.d0/sqrt(2.d0)
               call dsg4setc(hj(j),khc,ksqu(k),ksqd(l),auxc)
               call dsg4set(khc,hj(j),ksqu(k),ksqd(l),0.d0,0.d0)
               call dsg4set(hj(j),ksqu(k),khc,ksqd(l),0.d0,0.d0)
            enddo
         enddo
         ! Hj-Gc-u~-d~
         do k=1,6
            do l=1,6
               auxc=-conjg(rp(j))*gmudll(k,l)
               do m=1,3
                  do n=1,3
                     auxc=auxc
     &                    +conjg(sp(j))/mass(kw)**2*squlmx(k,m)*
     &                    rmass(kqu(m))**2*ckm(m,n)*conjg(sqdlmx(l,n))
     &                    +conjg(tp(j))/mass(kw)**2*squlmx(k,m)*
     &                    ckm(m,n)*rmass(kqd(n))**2*conjg(sqdlmx(l,n))
     &                    +2.d0*conjg(vp(j))/mass(kw)**2/s2b*
     &                    squrmx(k,m)*rmass(kqu(m))*ckm(m,n)*
     &                    rmass(kqd(n))*conjg(sqdrmx(l,n))
                  enddo
               enddo
               auxc=auxc*gsq/2.d0/sqrt(2.d0)
               call dsg4setc(hj(j),kgoldc,ksqu(k),ksqd(l),auxc)
               call dsg4set(kgoldc,hj(j),ksqu(k),ksqd(l),0.d0,0.d0)
               call dsg4set(hj(j),ksqu(k),kgoldc,ksqd(l),0.d0,0.d0)
            enddo
         enddo
      enddo
      do j=1,4 ! Hj=H1,H2,H3,G0
         ! Hj-Hc-nu~-e~
         do k=1,3
            do l=1,6
               auxc=-conjg(rr(j))*gmnell(k,l)
               do m=1,3
                  do n=1,3
                     auxc=auxc
     &                    +conjg(tt(j))/mass(kw)**2*slulmx(k,m)*
     &                    rmass(kl(n))**2*conjg(sldlmx(l,n))
                  enddo
               enddo
               auxc=auxc*gsq/2.d0/sqrt(2.d0)
               call dsg4setc(hj(j),khc,ksnu(k),ksl(l),auxc)
               call dsg4set(khc,hj(j),ksnu(k),ksl(l),0.d0,0.d0)
               call dsg4set(hj(j),ksnu(k),khc,ksl(l),0.d0,0.d0)
            enddo
         enddo
         ! Hj-Gc-nu~-e~
         do k=1,3
            do l=1,6
               auxc=-conjg(rp(j))*gmnell(k,l)
               do m=1,3
                  do n=1,3
                     auxc=auxc
     &                    +conjg(tp(j))/mass(kw)**2*slulmx(k,m)*
     &                    rmass(kl(n))**2*conjg(sldlmx(l,n))
                  enddo
               enddo
               auxc=auxc*gsq/2.d0/sqrt(2.d0)
               call dsg4setc(hj(j),kgoldc,ksnu(k),ksl(l),auxc)
               call dsg4set(kgoldc,hj(j),ksnu(k),ksl(l),0.d0,0.d0)
               call dsg4set(hj(j),ksnu(k),kgoldc,ksl(l),0.d0,0.d0)
            enddo
         enddo
      enddo

c--------------------------------------------------- the end of dsvertx3

      return
      end

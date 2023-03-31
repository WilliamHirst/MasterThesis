c====================================================================
c
c   this subroutine gives the imaginary part of the amplitude of the 
c   process of neutralino annihilation into two gluons in the limit 
c   of vanishing relative velocity of the neutralino pair
c
c   l. bergstrom & p. ullio, nucl. phys. b 504 (1997) 27
c
c   imres: imaginary part
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanglglim(imres)
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'
      integer kf,ksf,g,i
      complex*16 gaq,gan,gznc
      real*8 betaq,logg,mlsp,mf,mft,imres,sqre,dqre,sum,quarks,
     &  gaqre,gzn,tst,dsabsq,t3f,signm,rootab,aq,bq
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      mlsp=mass(kn(1))
      signm=mlsp/dabs(mlsp)
      sum=0.d0
c
cccccccccc  first contribution to imres  cccccccccccccccccccccccccccc
c
      do 20 g=1,3
c  up-type quarks
c  define particle code:
        if(g.eq.1) then
          kf=ku
        elseif(g.eq.2) then
          kf=kc
        else ! if(g.eq.3) then
          kf=kt
        endif
c  and define particle parameters
        t3f=.5d0
        mf=mass(kf)
        quarks=0.d0
        if (mf.ge.dabs(mlsp)) goto 20
        betaq=dsqrt(1.d0-(mf/mlsp)**2)
        logg=dlog((1.d0+betaq)/(1.d0-betaq))
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        gzn=dreal(gznc)*g2weak/costhw*t3f
        tst=dimag(gznc)
        if (abs(tst).ge.1.d-10) then
          write(*,*) 'complex couplings!!! glglim tst, model no. ',
     &      idtag
          stop
        endif
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! glglim gaq*gan, model no. '
     &      ,idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
c  box part:
        do 10 i=1,6
c  up-type squarks
          ksf=ksqu(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &            +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          rootab=signm*dsqrt(dabs(aq*bq))
          quarks=quarks+0.5*(sqre*bq+dqre*rootab)/(1.d0+aq-bq)
 10     continue
c  triangular part:
        quarks=quarks
     &    -mf*mlsp/4.d0*gaqre*(4.d0*mlsp**2-mass(kh3)**2)
     &     /((4.d0*mlsp**2-mass(kh3)**2)**2
     &     +width(kh3)**2*mass(kh3)**2)
     &    +mf**2/mass(kz)**2*gzn/4.d0
     &     *(4.0d0*mlsp**2-mass(kz)**2)**2
     &     /((4.0d0*mlsp**2-mass(kz)**2)**2
     &     +width(kz)**2*mass(kz)**2)
        sum=sum-pi*logg*quarks/mlsp**2
 20   continue
c
cccccccccc  second contribution to imres  ccccccccccccccccccccccccccc
c
      do 21 g=1,3
c  down-type quarks
c  define particle code:
        if(g.eq.1) then
          kf=kd
        elseif(g.eq.2) then
          kf=ks
        else ! if(g.eq.3) then
          kf=kb
        endif
c  and define particle parameters
        t3f=-.5d0
        mf=mass(kf)
        quarks=0.d0
        if (mf.ge.dabs(mlsp)) goto 21
        betaq=dsqrt(1.d0-(mf/mlsp)**2)
        logg=dlog((1.d0+betaq)/(1.d0-betaq))
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        tst=dimag(gznc)
        if (abs(tst).ge.1.d-10) then
          write(*,*) 'complex couplings!!! glglim tst, model no. ',
     &      idtag
          stop
        endif
        gzn=dreal(gznc)*g2weak/costhw*t3f
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! glglim gaq*gan, model no. ',
     &      idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
c  box part:
        do 11 i=1,6
c  down-type squarks
          ksf=ksqd(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &            +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          rootab=signm*dsqrt(dabs(aq*bq))
          quarks=quarks+0.5*(sqre*bq+dqre*rootab)/(1.d0+aq-bq)
 11     continue
c  triangular part:
        quarks=quarks
     &    -mf*mlsp/4.d0*gaqre*(4.d0*mlsp**2-mass(kh3)**2)
     &     /((4.d0*mlsp**2-mass(kh3)**2)**2
     &     +width(kh3)**2*mass(kh3)**2)
     &    +mf**2/mass(kz)**2*gzn/4.d0
     &     *(4.0d0*mlsp**2-mass(kz)**2)**2
     &     /((4.0d0*mlsp**2-mass(kz)**2)**2
     &     +width(kz)**2*mass(kz)**2)
        sum=sum-pi*logg*quarks/mlsp**2
 21   continue
      imres=sum
      end













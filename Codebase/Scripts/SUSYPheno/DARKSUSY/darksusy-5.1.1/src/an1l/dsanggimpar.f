c====================================================================
c
c   this subroutine gives the imaginary part of the amplitude of the 
c   process of neutralino annihilation into two photons in the limit 
c   of vanishing relative velocity of the neutralino pair
c
c   l. bergstrom & p. ullio, nucl. phys. b 504 (1997) 27
c
c   see header of dsanggim.f for details
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanggimpar(imres,imfbx,imftx,imgbx)
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'
      real*8 imres,imfbx,imftx,imgbx
      integer g,i
      integer kf,ksf,kkc
      complex*16 gaq,gan,gznc
      real*8 mlsp,mf,mft,betaq,logg,sqre,dqre,gaqre,gzn
     &  ,dsabsq,t3f,signm,rootab,aq,bq,eqs,col,mw,mc,par1,par2
     &  ,e7temp,e7temp2,h3wid,zwid
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      mlsp=mass(kn(1))
      signm=mlsp/dabs(mlsp)
      h3wid=width(kh3)
      if(h3wid.lt.1.d0) h3wid=1.d0
      zwid=width(kz)
      if(zwid.lt.1.d0) zwid=1.d0
      imfbx=0.d0
      imftx=0.d0
      imgbx=0.d0
c
cccccccccc  first contribution to imfbx and imftx  cccccccccccccccccc
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
        eqs=4.d0/9.d0
        col=3.d0
        mf=mass(kf)
c
        par1=0.d0
        par2=0.d0
        if (mf.ge.dabs(mlsp)) goto 20
        betaq=dsqrt(1.d0-(mf/mlsp)**2)
        logg=dlog((1.d0+betaq)/(1.d0-betaq))
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        gzn=dreal(gznc)*g2weak/costhw*t3f
        if (abs(dimag(gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gznc, model no. ',idtag
          stop
        endif
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gaq*gan, model no. ',idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
c  box part:
        do 10 i=1,6
c  up-type squarks
          ksf=ksqu(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &          +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          rootab=signm*dsqrt(dabs(aq*bq))
          par1=par1+.5d0*(sqre*bq+dqre*rootab)/(1.d0+aq-bq)
 10       continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=-mf*mlsp/4.d0*gaqre*e7temp
     &    +mf**2/mass(kz)**2*gzn/4.d0*e7temp2
        imfbx=imfbx-pi*logg*eqs*col/mlsp**2*par1
        imftx=imftx-pi*logg*eqs*col/mlsp**2*par2
 20   continue
c
cccccccccc  second contribution to imfbx and imftx  ccccccccccccccccc
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
        eqs=1.d0/9.d0
        col=3.d0
        mf=mass(kf)
c
        par1=0.d0
        par2=0.d0
        if (mf.ge.dabs(mlsp)) goto 21
        betaq=dsqrt(1.d0-(mf/mlsp)**2)
        logg=dlog((1.d0+betaq)/(1.d0-betaq))
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        gzn=dreal(gznc)*g2weak/costhw*t3f
        if (abs(dimag(gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gznc, model no. ',idtag
          stop
        endif
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gaq*gan, model no. ',idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
c  box part:
        do 11 i=1,6
c  down-type squarks
          ksf=ksqd(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &          +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          rootab=signm*dsqrt(dabs(aq*bq))
          par1=par1+.5d0*(sqre*bq+dqre*rootab)/(1.d0+aq-bq)
 11     continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=-mf*mlsp/4.d0*gaqre*e7temp+
     &    mf**2/mass(kz)**2*gzn/4.d0*e7temp2
        imfbx=imfbx-pi*logg*eqs*col/mlsp**2*par1
        imftx=imftx-pi*logg*eqs*col/mlsp**2*par2
 21   continue
c
cccccccccc  third contribution to imfbx and imftx  cccccccccccccccccc
c
      do 22 g=1,3
c  charged leptons
c  define particle code:
        if(g.eq.1) then
          kf=ke
        elseif(g.eq.2) then
          kf=kmu
        else ! if(g.eq.3) then
          kf=ktau
        endif
c  and define particle parameters
        t3f=-.5d0
        eqs=1.d0
        col=1.d0
        mf=mass(kf)
c
        par1=0.d0
        par2=0.d0
        if (mf.ge.dabs(mlsp)) goto 22
        betaq=dsqrt(1.d0-(mf/mlsp)**2)
        logg=dlog((1.d0+betaq)/(1.d0-betaq))
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        gzn=dreal(gznc)*g2weak/costhw*t3f
        if (abs(dimag(gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gznc, model no. ',idtag
          stop
        endif
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gaq*gan, model no. ',idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
c  box part:
        do 12 i=1,6
c  charged sleptons
          ksf=ksl(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &          +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          rootab=signm*dsqrt(dabs(aq*bq))
          par1=par1+.5d0*(sqre*bq+dqre*rootab)/(1.d0+aq-bq)
 12     continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=-mf*mlsp/4.d0*gaqre*e7temp+
     &    mf**2/mass(kz)**2*gzn/4.d0*e7temp2
        imfbx=imfbx-pi*logg*eqs*col/mlsp**2*par1
        imftx=imftx-pi*logg*eqs*col/mlsp**2*par2
 22   continue
c
cccccccccc  contribution to imgbx  cccccccccccccccccccccccccccccccccc
c
      mw=mass(kw)
      par1=0.0d0
      if (mw.ge.dabs(mlsp)) goto 23
      betaq=dsqrt(1.d0-(mw/mlsp)**2)
      logg=dlog((1.d0+betaq)/(1.d0-betaq))
      do 13 i=1,2
        kkc=kcha(i)
        mc=mass(kkc)
        sqre=(dsabsq(gl(kw,kn(1),kkc))
     &        +dsabsq(gr(kw,kn(1),kkc)))/2.d0
        dqre=dreal(gl(kw,kn(1),kkc)*conjg(gr(kw,kn(1),kkc)))
        aq=(mlsp/mc)**2
        bq=(mw/mc)**2
        par1=par1+2.d0*(aq-bq)*sqre/(1.d0+aq-bq)
 13   continue
      imgbx=imgbx-pi*logg/mlsp**2*par1
 23   continue
c
cccccccccc  total and partial results  cccccccccccccccccccccccccccccc
c
      imres=imfbx+imftx+imgbx
      imfbx=imfbx/imres
      imftx=imftx/imres
      imgbx=imgbx/imres
      return
      end

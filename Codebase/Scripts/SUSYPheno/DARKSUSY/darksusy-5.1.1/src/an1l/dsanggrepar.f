c====================================================================
c
c   this subroutine gives the imaginary part of the amplitude of the 
c   process of neutralino annihilation into two photons in the limit 
c   of vanishing relative velocity of the neutralino pair
c
c   l. bergstrom & p. ullio, nucl. phys. b 504 (1997) 27
c
c   see header of dsanggre.f for details
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanggrepar(reres,refbx,reftx,rehbx,rehtx,regbx)
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'
      real*8 reres,refbx,reftx,rehbx,rehtx,regbx
      integer g,i
      integer kf,ksf,kkc
      complex*16 gaq,gan,gznc,gzq
      real*8 mlsp,mf,mft,mhc,mc,mw,signm,aq,bq,sqre,dqre,
     &  gzn,gaqre,t3f,col,eqs,gzqre,dsrepfbox,dsrepw,dsrepgh,
     &  dspi1,dsabsq,e7temp,e7temp2,h3wid,zwid,par1,par2
c
      mlsp=mass(kn(1))
      signm=mlsp/dabs(mlsp)
      h3wid=width(kh3)
      if(h3wid.lt.1.d0) h3wid=1.d0
      zwid=width(kz)
      if(zwid.lt.1.d0) zwid=1.d0
      reres=0.d0
      refbx=0.d0
      reftx=0.d0
      rehbx=0.d0
      rehtx=0.d0
      regbx=0.d0
c
cccccccccc  first contribution to refbx and reftx  cccccccccccccccccc
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
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        if (abs(dimag(gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gznc, model no. ',idtag
          stop
        endif
        gzn=dreal(gznc)*g2weak/costhw*t3f
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
          par1=par1+dsrepfbox(aq,bq,sqre,dqre,signm)
 10       continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=(-mf*mlsp/4.d0*gaqre*e7temp
     &    +mf**2/mass(kz)**2*gzn/4.d0*e7temp2)*dspi1(aq,bq)
        refbx=refbx+eqs*col/mlsp**2*par1
        reftx=reftx+eqs*col/mlsp**2*par2
 20   continue
c
cccccccccc  second contribution to refbx and reftx  ccccccccccccccccc
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
        gaq=gl(kh3,kf,kf)-gr(kh3,kf,kf)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        gznc=gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1))
        if (abs(dimag(gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gznc, model no. ',idtag
          stop
        endif
        gzn=dreal(gznc)*g2weak/costhw*t3f
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
          par1=par1+dsrepfbox(aq,bq,sqre,dqre,signm)
 11       continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=(-mf*mlsp/4.d0*gaqre*e7temp
     &    +mf**2/mass(kz)**2*gzn/4.d0*e7temp2)*dspi1(aq,bq)
        refbx=refbx+eqs*col/mlsp**2*par1
        reftx=reftx+eqs*col/mlsp**2*par2
 21   continue
c
cccccccccc  third contribution to refbx and reftx  cccccccccccccccccc
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
c  charged slepton
          ksf=ksl(i)
          mft=mass(ksf)
          sqre=(dsabsq(gl(ksf,kf,kn(1)))
     &          +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          dqre=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          par1=par1+dsrepfbox(aq,bq,sqre,dqre,signm)
 12     continue
c  triangular part:
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=(-mf*mlsp/4.d0*gaqre*e7temp
     &    +mf**2/mass(kz)**2*gzn/4.d0*e7temp2)*dspi1(aq,bq)
        refbx=refbx+eqs*col/mlsp**2*par1
        reftx=reftx+eqs*col/mlsp**2*par2
 22   continue
c
cccccccccc  contribution to rehbx and rehtx  cccccccccccccccccccccccc
c
      mhc=mass(khc)
      par1=0.d0
      par2=0.d0
      do 23 i=1,2
        kkc=kcha(i)
        mc=mass(kkc)
        signm=mlsp/dabs(mlsp)*mc/dabs(mc)
        sqre=(dsabsq(gl(khc,kn(1),kkc))
     &        +dsabsq(gr(khc,kn(1),kkc)))/2.d0
        dqre=dreal(gl(khc,kn(1),kkc)*conjg(gr(khc,kn(1),kkc)))
        gaq=gl(kh3,kkc,kkc)-gr(kh3,kkc,kkc)
        gan=gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1))
        if (abs(dimag(gaq*gan)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gaq*gan, model no. ',idtag
          stop
        endif
        gaqre=dreal(gaq*gan)
        gzq=gl(kz,kkc,kkc)-gr(kz,kkc,kkc)
        gznc=gr(kz,kn(1),kn(1))-gl(kz,kn(1),kn(1))
        if (abs(dimag(gzq*gznc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gaq*gan, model no. ',idtag
          stop
        endif
        gzqre=dreal(gzq*gznc)
        aq=(mlsp/mhc)**2
        bq=(mc/mhc)**2
c  box part:
        par1=par1+dsrepfbox(aq,bq,sqre,dqre,signm)
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mass(kz)**2)**2/
     &    ((4.0d0*mlsp**2-mass(kz)**2)**2+zwid**2*mass(kz)**2)
        par2=par2+(-mc*mlsp/4.d0*gaqre*e7temp
     &    +mc**2/mass(kz)**2*gzqre/4.d0*e7temp2)*dspi1(aq,bq)
 23   continue
      rehbx=rehbx+par1/mlsp**2
      rehtx=rehtx+par2/mlsp**2
c
cccccccccc  first contribution to regbx  cccccccccccccccccccccccccccc
c
      mw=mass(kw)
      par1=0.d0
      do 24 i=1,2
        kkc=kcha(i)
        mc=mass(kkc)
        signm=mlsp/dabs(mlsp)*mc/dabs(mc)
        sqre=(dsabsq(gl(kw,kn(1),kkc))
     &        +dsabsq(gr(kw,kn(1),kkc)))/2.d0
        dqre=dreal(gl(kw,kn(1),kkc)*conjg(gr(kw,kn(1),kkc)))
        aq=(mlsp/mc)**2
        bq=(mw/mc)**2
        par1=par1+dsrepw(aq,bq,sqre,dqre,signm)
 24   continue
      regbx=regbx+par1/mlsp**2
c
cccccccccc  second contribution to regbx  ccccccccccccccccccccccccccc
c
      par1=0.d0
      do 25 i=1,2
        kkc=kcha(i)
        mc=mass(kkc)
        signm=mlsp/dabs(mlsp)*mc/dabs(mc)
        sqre=(dsabsq(gl(kgoldc,kn(1),kkc))
     &        +dsabsq(gr(kgoldc,kn(1),kkc)))/2.d0
        dqre=dreal(gl(kgoldc,kn(1),kkc)*conjg(gr(kgoldc,kn(1),kkc)))
        aq=(mlsp/mc)**2
        bq=(mw/mc)**2
        par1=par1+dsrepgh(aq,bq,sqre,dqre,signm)
 25   continue
      regbx=regbx+par1/mlsp**2
c
cccccccccc  total and partial results  cccccccccccccccccccccccccccccc
c
      reres=refbx+reftx+rehbx+rehtx+regbx
      refbx=refbx/reres
      reftx=reftx/reres
      rehbx=rehbx/reres
      rehtx=rehtx/reres
      regbx=regbx/reres
      return
      end



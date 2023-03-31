c====================================================================
c
c   this subroutine gives the real and imaginary parts of the 
c   amplitude of the process of neutralino annihilation into 
c   one photon and one z boson in the limit of vanishing relative 
c   velocity of the neutralino pair
c
c   p. ullio & l. bergstrom, phys. rev. d 57 (1998) 1962
c
c   see header of dsanzg.f for details
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanzgpar(imres,imfbx,imftx,imgbx,reres,refbx,reftx,
     &  rehbx,rehtx,regbx)
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'
      real*8 imres,imfbx,imftx,imgbx,reres,refbx,reftx,
     &  rehbx,rehtx,regbx,t3f
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c  functions
      real*8 dsi_14,dsi_13,dsi_24,dsi_23,dsi_34,dsi_33,
     &  dsti_214,dsti_224,
     &  dsti_23,dsti_33,dsi_41,dsi_42,dsj_1,dsj_2,dsj_3
c  couplings
      real*8 gzf,ghf,df,saf,sbf,scf,sah,sbh,dbh,sch,dch,
     &  sz,dz,sh,dh,saw,daw,sbw,dbw,sag,sdg,ddg,sbg,dbg,scg,dcg,dsabsq
      real*8 zcoup(6)
      complex*16 gzfc,ghfc,dfc,safc,sbfc,scfc,sahc,sbhc,dbhc,
     &  schc,dchc,szc,dzc,shc,dhc,sawc,dawc,sbwc,dbwc,sagc,sdgc,
     &  ddgc,sbgc,dbgc,scgc,dcgc
c  coefficients
      real*8 e1,e2,e3,e4,e5,e6,e7,f1,f2,f3,f4,f5,f6,f7,f8,
     &  f9,f10,f11,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,h1,h2,h3,h4,h5,
     &  h6,h7,h8,h9,h10,e7temp,e7temp2,f10temp,f10temp2
c  masses
      real*8 mlsp,mf,mft,mz,mhc,mci,mcj,mw
c  others
      real*8 eq,col,aq,bq,cq,dq,signm,rootab,signmi,
     &  signmj,roota,rootd,rootad,rootc,h3wid,zwid
      integer kf,ksf,kkci,kkcj
      integer g,i,j
      real*8 par1r,par1i,par2r,par2i
c
c
      imfbx=0.d0
      imftx=0.d0
      imgbx=0.d0
      refbx=0.d0
      reftx=0.d0
      rehbx=0.d0
      rehtx=0.d0
      regbx=0.d0
      mlsp=mass(kn(1))
      signm=mlsp/(dabs(mlsp))
      h3wid=width(kh3)
      if(h3wid.lt.1.d0) h3wid=1.d0
      mz=mass(kz)
      zwid=width(kz)
      if(zwid.lt.1.d0) zwid=1.d0
c
cccc  first contribution to imfbx, imftx, refbx and reftx  cccccccccc
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
        eq=2.d0/3.d0
        col=3.d0
        mf=mass(kf)
        par1r=0.d0
        par1i=0.d0
        par2r=0.d0
        par2i=0.d0
c  box part:
        do 10 i=1,6
c  up-type squarks
          ksf=ksqu(i)
          mft=mass(ksf)
c  up-type quarks
          if(i.le.3) then
	    zcoup(i)=dreal(gl(kz,kqu(i),kqu(i)))
          else
	    zcoup(i)=dreal(gr(kz,kqu(i-3),kqu(i-3)))
          endif
          dfc=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
     &      *(gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
          df=dreal(dfc)
          if (abs(dimag(dfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! dfc, model no. ',idtag
            stop
          endif
          safc=zcoup(i)*(dsabsq(gl(ksf,kf,kn(1)))
     &      +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          saf=dreal(safc)
          if (abs(dimag(safc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! saf, model no. ',idtag
            stop
          endif
          sbfc=(dsabsq(gl(ksf,kf,kn(1)))*gr(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gl(kz,kf,kf))/2.d0
          sbf=dreal(sbfc)
          if (abs(dimag(sbfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! sbf, model no. ',idtag
            stop
          endif
          scfc=(dsabsq(gl(ksf,kf,kn(1)))*gl(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gr(kz,kf,kf))/2.d0
          scf=dreal(scfc)
          if (abs(dimag(scfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! scf, model no. ',idtag
            stop
          endif
c  definition of a,b,c
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          cq=(mz/mft)**2
c  definition of rootab
          rootab=signm*dsqrt(dabs(aq*bq))
c  definitions of the coefficients for box diagrams
          e1=-((bq*(sbf+scf)+2.d0*rootab*df)/(4.d0*(1.d0+aq-bq)))
          e2=-(((1.d0-cq/4.d0)*saf)/(1.d0-bq-cq/4.d0)+
     &      sbf/(1.d0-bq+cq/4.d0))/4.d0
          e3=((bq*saf)/(1.d0-bq-cq/4.d0))/4.d0
          e4=(((bq-cq/4.d0)*sbf)/(1.d0-bq+cq/4.d0))/4.d0
          e5=-(((1.d0-bq-cq/4.d0)*saf)/(aq-cq/4.d0)**2)*cq/16.d0
          e6=-(((1.d0-bq+cq/4.d0)*sbf)/(aq-cq/4.d0)**2)*cq/16.d0
c  real part:
          par1r=par1r+(e1*dsi_13(aq,bq,cq/4.d0)+
     &      (e2+e5)*dsti_23(aq,bq,cq/4.d0)+(e1+e3)*
     &      dsi_33(aq,bq,cq/4.d0)+
     &      (e1+e4+e6)*dsti_33(aq,bq,cq/4.d0)+
     &      saf*dsi_42(aq,bq,1.d0,cq)/4.d0-sbf*dsi_41(aq,bq,cq)/4.d0)
c  imaginary part:
          if (mf.ge.dabs(mlsp)) goto 30
          par1i=par1i-e1*dsj_1(aq,bq)
 30       if (2.d0*mf.ge.mz) goto 10
          par1i=par1i+(e1*dsj_2(bq,cq)+(e1+e4+e6)*dsj_3(aq,bq,cq)+
     &      cq/8.d0*(sqrt(1.d0-4.d0*bq/cq)*sbf)/(aq-cq/4.d0))
 10     continue
c  triangular part:
        gzfc=(gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1)))*
     &    ((gl(kz,kf,kf))**2-(gr(kz,kf,kf))**2)/2.d0
        gzf=dreal(gzfc)
        if (abs(dimag(gzfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gzf, model no. ',idtag
          stop
        endif
        ghfc=(gl(kh3,kn(1),kn(1))-gr(kh3,kn(1),kn(1)))*
     &    (gl(kh3,kf,kf)-gr(kh3,kf,kf))*
     &    (gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
        ghf=dreal(ghfc)
        if (abs(dimag(ghfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! ghf, model no. ',idtag
          stop
        endif
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mz**2)**2/
     &    ((4.0d0*mlsp**2-mz**2)**2+zwid**2*mz**2)
        e7=(mf**2*gzf)/(4.d0*mz**2)*e7temp2
     &     -(mf*mlsp*ghf)*e7temp/4.d0
c  real part:
        par2r=par2r+e7*dsi_13(aq,bq,cq/4.d0)
c  imaginary part:
        if (mf.ge.dabs(mlsp)) goto 40
        par2i=par2i-e7*dsj_1(aq,bq)
 40     if ((2.d0*mf).ge.mz) goto 50
        par2i=par2i+e7*dsj_2(bq,cq)
 50     continue
        refbx=refbx+eq*col*par1r
        imfbx=imfbx+eq*col*par1i
        reftx=reftx+eq*col*par2r
        imftx=imftx+eq*col*par2i
 20   continue
c
cccc  second contribution to imfbx, imftx, refbx and reftx  ccccccccc
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
        eq=-1.d0/3.d0
        col=3.d0
        mf=mass(kf)
        par1r=0.d0
        par1i=0.d0
        par2r=0.d0
        par2i=0.d0
c  box part:
        do 11 i=1,6
c  down-type squarks
          ksf=ksqd(i)
          mft=mass(ksf)
c  down-type quarks
          if(i.le.3) then
	    zcoup(i)=dreal(gl(kz,kqd(i),kqd(i)))
          else
	    zcoup(i)=dreal(gr(kz,kqd(i-3),kqd(i-3)))
          endif
          dfc=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
     &      *(gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
          df=dreal(dfc)
          if (abs(dimag(dfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! dfc, model no. ',idtag
            stop
          endif
          safc=zcoup(i)*(dsabsq(gl(ksf,kf,kn(1)))
     &      +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          saf=dreal(safc)
          if (abs(dimag(safc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! saf, model no. ',idtag
            stop
          endif
          sbfc=(dsabsq(gl(ksf,kf,kn(1)))*gr(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gl(kz,kf,kf))/2.d0
          sbf=dreal(sbfc)
          if (abs(dimag(sbfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! sbf, model no. ',idtag
            stop
          endif
          scfc=(dsabsq(gl(ksf,kf,kn(1)))*gl(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gr(kz,kf,kf))/2.d0
          scf=dreal(scfc)
          if (abs(dimag(scfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! scf, model no. ',idtag
            stop
          endif
c  definition of a,b,c
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          cq=(mz/mft)**2
c  definition of rootab
          rootab=signm*dsqrt(dabs(aq*bq))
c  definitions of the coefficients for box diagrams
          e1=-((bq*(sbf+scf)+2.d0*rootab*df)/(4.d0*(1.d0+aq-bq)))
          e2=-(((1.d0-cq/4.d0)*saf)/(1.d0-bq-cq/4.d0)+
     &      sbf/(1.d0-bq+cq/4.d0))/4.d0
          e3=((bq*saf)/(1.d0-bq-cq/4.d0))/4.d0
          e4=(((bq-cq/4.d0)*sbf)/(1.d0-bq+cq/4.d0))/4.d0
          e5=-(((1.d0-bq-cq/4.d0)*saf)/(aq-cq/4.d0)**2)*cq/16.d0
          e6=-(((1.d0-bq+cq/4.d0)*sbf)/(aq-cq/4.d0)**2)*cq/16.d0
c  real part:
          par1r=par1r+(e1*dsi_13(aq,bq,cq/4.d0)+
     &      (e2+e5)*dsti_23(aq,bq,cq/4.d0)+(e1+e3)*
     &      dsi_33(aq,bq,cq/4.d0)+
     &      (e1+e4+e6)*dsti_33(aq,bq,cq/4.d0)+
     &      saf*dsi_42(aq,bq,1.d0,cq)/4.d0-sbf*dsi_41(aq,bq,cq)/4.d0)
c  imaginary part:
          if (mf.ge.dabs(mlsp)) goto 31
          par1i=par1i-e1*dsj_1(aq,bq)
 31       if (2.d0*mf.ge.mz) goto 11
          par1i=par1i+(e1*dsj_2(bq,cq)+(e1+e4+e6)*dsj_3(aq,bq,cq)+
     &      cq/8.d0*(sqrt(1.d0-4.d0*bq/cq)*sbf)/(aq-cq/4.d0))
 11     continue
c  triangular part:
        gzfc=(gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1)))*
     &    ((gl(kz,kf,kf))**2-(gr(kz,kf,kf))**2)/2.d0
        gzf=dreal(gzfc)
        if (abs(dimag(gzfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gzf, model no. ',idtag
          stop
        endif
        ghfc=(gl(kh3,kn(1),kn(1))-gr(kh3,kn(1),kn(1)))*
     &    (gl(kh3,kf,kf)-gr(kh3,kf,kf))*
     &    (gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
        ghf=dreal(ghfc)
        if (abs(dimag(ghfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! ghf, model no. ',idtag
          stop
        endif
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mz**2)**2/
     &    ((4.0d0*mlsp**2-mz**2)**2+zwid**2*mz**2)
        e7=(mf**2*gzf)/(4.d0*mz**2)*e7temp2
     &     -(mf*mlsp*ghf)*e7temp/4.d0
c  real part:
        par2r=par2r+e7*dsi_13(aq,bq,cq/4.d0)
c  imaginary part:
        if (mf.ge.dabs(mlsp)) goto 41
        par2i=par2i-e7*dsj_1(aq,bq)
 41     if ((2.d0*mf).ge.mz) goto 51
        par2i=par2i+e7*dsj_2(bq,cq)
 51     continue
        refbx=refbx+eq*col*par1r
        imfbx=imfbx+eq*col*par1i
        reftx=reftx+eq*col*par2r
        imftx=imftx+eq*col*par2i
 21   continue
c
cccc  third contribution to imfbx, imftx, refbx and reftx  cccccccccc
c
      do 22 g=1,3
c  charged lepton
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
c bug fixed: eq=-1 not eq=1        !Piero Ullio 00-12-04
c        eq=1.d0
        eq=-1.d0
        col=1.d0
        mf=mass(kf)
        par1r=0.d0
        par1i=0.d0
        par2r=0.d0
        par2i=0.d0
c  box part:
        do 12 i=1,6
c  charged slepton
          ksf=ksl(i)
          mft=mass(ksf)
c  charged leptons
          if(i.le.3) then
	    zcoup(i)=dreal(gl(kz,kl(i),kl(i)))
          else
	    zcoup(i)=dreal(gr(kz,kl(i-3),kl(i-3)))
          endif
          dfc=dreal(gl(ksf,kf,kn(1))*conjg(gr(ksf,kf,kn(1))))
     &      *(gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
          df=dreal(dfc)
          if (abs(dimag(dfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! dfc, model no. ',idtag
            stop
          endif
          safc=zcoup(i)*(dsabsq(gl(ksf,kf,kn(1)))
     &      +dsabsq(gr(ksf,kf,kn(1))))/2.d0
          saf=dreal(safc)
          if (abs(dimag(safc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! saf, model no. ',idtag
            stop
          endif
          sbfc=(dsabsq(gl(ksf,kf,kn(1)))*gr(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gl(kz,kf,kf))/2.d0
          sbf=dreal(sbfc)
          if (abs(dimag(sbfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! sbf, model no. ',idtag
            stop
          endif
          scfc=(dsabsq(gl(ksf,kf,kn(1)))*gl(kz,kf,kf)+
     &      dsabsq(gr(ksf,kf,kn(1)))*gr(kz,kf,kf))/2.d0
          scf=dreal(scfc)
          if (abs(dimag(scfc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!! scf, model no. ',idtag
            stop
          endif
c  definition of a,b,c
          aq=(mlsp/mft)**2
          bq=(mf/mft)**2
          cq=(mz/mft)**2
c  definition of rootab
          rootab=signm*dsqrt(dabs(aq*bq))
c  definitions of the coefficients for box diagrams
          e1=-((bq*(sbf+scf)+2.d0*rootab*df)/(4.d0*(1.d0+aq-bq)))
          e2=-(((1.d0-cq/4.d0)*saf)/(1.d0-bq-cq/4.d0)+
     &      sbf/(1.d0-bq+cq/4.d0))/4.d0
          e3=((bq*saf)/(1.d0-bq-cq/4.d0))/4.d0
          e4=(((bq-cq/4.d0)*sbf)/(1.d0-bq+cq/4.d0))/4.d0
          e5=-(((1.d0-bq-cq/4.d0)*saf)/(aq-cq/4.d0)**2)*cq/16.d0
          e6=-(((1.d0-bq+cq/4.d0)*sbf)/(aq-cq/4.d0)**2)*cq/16.d0
c  real part:
          par1r=par1r+(e1*dsi_13(aq,bq,cq/4.d0)+
     &      (e2+e5)*dsti_23(aq,bq,cq/4.d0)+
     &      (e1+e3)*dsi_33(aq,bq,cq/4.d0)+
     &      (e1+e4+e6)*dsti_33(aq,bq,cq/4.d0)+
     &      saf*dsi_42(aq,bq,1.d0,cq)/4.d0-sbf*dsi_41(aq,bq,cq)/4.d0)
c  imaginary part:
          if (mf.ge.dabs(mlsp)) goto 32
          par1i=par1i-e1*dsj_1(aq,bq)
 32       if (2.d0*mf.ge.mz) goto 12
          par1i=par1i+(e1*dsj_2(bq,cq)+(e1+e4+e6)*dsj_3(aq,bq,cq)+
     &      cq/8.d0*(sqrt(1.d0-4.d0*bq/cq)*sbf)/(aq-cq/4.d0))
 12       continue
c  triangular part:
        gzfc=(gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1)))*
     &    ((gl(kz,kf,kf))**2-(gr(kz,kf,kf))**2)/2.d0
        gzf=dreal(gzfc)
        if (abs(dimag(gzfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! gzf, model no. ',idtag
          stop
        endif
        ghfc=(gl(kh3,kn(1),kn(1))-gr(kh3,kn(1),kn(1)))*
     &    (gl(kh3,kf,kf)-gr(kh3,kf,kf))*
     &    (gl(kz,kf,kf)+gr(kz,kf,kf))/2.d0
        ghf=dreal(ghfc)
        if (abs(dimag(ghfc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!! ghf, model no. ',idtag
          stop
        endif
        e7temp=(4.d0*mlsp**2-mass(kh3)**2)/
     &    ((4.d0*mlsp**2-mass(kh3)**2)**2+mass(kh3)**2*h3wid**2)
        e7temp2=(4.0d0*mlsp**2-mz**2)**2/
     &    ((4.0d0*mlsp**2-mz**2)**2+zwid**2*mz**2)
        e7=(mf**2*gzf)/(4.d0*mz**2)*e7temp2
     &     -(mf*mlsp*ghf)*e7temp/4.d0
c  real part:
        par2r=par2r+e7*dsi_13(aq,bq,cq/4.d0)
c  imaginary part:
        if (mf.ge.dabs(mlsp)) goto 42
        par2i=par2i-e7*dsj_1(aq,bq)
 42     if ((2.d0*mf).ge.mz) goto 52
        par2i=par2i+e7*dsj_2(bq,cq)
 52     continue
        refbx=refbx+eq*col*par1r
        imfbx=imfbx+eq*col*par1i
        reftx=reftx+eq*col*par2r
        imftx=imftx+eq*col*par2i
 22   continue
c
cccccccccc  contribution to rehbx and rehtx  cccccccccccccccccccccccc
c
      mhc=mass(khc)
      do 23 i=1,2
        par1r=0.d0
        par2r=0.d0
        kkci=kcha(i)
        mci=mass(kkci)
c  diagonal coupling
        sahc=gl(kz,khc,khc)*(dsabsq(gl(khc,kn(1),kkci))+
     &     dsabsq(gr(khc,kn(1),kkci)))/2.d0
        sah=dreal(sahc)
        if (abs(dimag(sahc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
c  definition of a,b,c and signmi
        aq=(mlsp/mci)**2
        bq=(mhc/mci)**2
        cq=(mz/mci)**2
        signmi=mci/(dabs(mci))
c  diagonal coefficients
        f1=-(sah/(1.d0-bq+cq/4.d0))/4.d0
        f2=((sah*(bq-cq/4.d0))/(1.d0-bq+cq/4.d0))/4.d0
        f3=((sah*(1.d0-bq+cq/4.d0))/(aq-cq/4.d0)**2)*cq/16.d0
c  diagonal contributions
c  box part:
        par1r=par1r+(f1*dsi_23(aq,bq,cq/4.d0)+
     &    (f2+f3)*dsti_33(aq,bq,cq/4.d0)+sah/4.d0*dsi_41(aq,bq,cq))
        do 13 j=1,2
          kkcj=kcha(j)
          mcj=mass(kkcj)
c  non-diagonal couplings
          sbhc=(gl(khc,kn(1),kkcj)*conjg(gl(khc,kn(1),kkci))*
     &      gl(kz,kkcj,kkci)+gr(khc,kn(1),kkcj)*
     &      conjg(gr(khc,kn(1),kkci))*gr(kz,kkcj,kkci))/2.d0
          sbh=dreal(sbhc)
          if (abs(dimag(sbhc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dbhc=(gl(khc,kn(1),kkcj)*conjg(gr(khc,kn(1),kkci))*
     &      gl(kz,kkcj,kkci)+gr(khc,kn(1),kkcj)*
     &      conjg(gl(khc,kn(1),kkci))*gr(kz,kkcj,kkci))/2.d0
          dbh=dreal(dbhc)
          if (abs(dimag(dbhc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          schc=(gl(khc,kn(1),kkcj)*conjg(gl(khc,kn(1),kkci))*
     &      gr(kz,kkcj,kkci)+gr(khc,kn(1),kkcj)*
     &      conjg(gr(khc,kn(1),kkci))*gl(kz,kkcj,kkci))/2.d0
          sch=dreal(schc)
          if (abs(dimag(schc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dchc=(gl(khc,kn(1),kkcj)*conjg(gr(khc,kn(1),kkci))*
     &      gr(kz,kkcj,kkci)+gr(khc,kn(1),kkcj)*
     &      conjg(gl(khc,kn(1),kkci))*gl(kz,kkcj,kkci))/2.d0
          dch=dreal(dchc)
          if (abs(dimag(dchc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          szc=(gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1)))*
     &      (gl(kz,kkci,kkcj)*gl(kz,kkcj,kkci)-
     &      gr(kz,kkci,kkcj)*gr(kz,kkcj,kkci))
          sz=dreal(szc)
          if (abs(dimag(szc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dzc=(gl(kz,kn(1),kn(1))-gr(kz,kn(1),kn(1)))*
     &      (gl(kz,kkci,kkcj)*gr(kz,kkcj,kkci)-
     &      gr(kz,kkci,kkcj)*gl(kz,kkcj,kkci))
          dz=dreal(dzc)
          if (abs(dimag(dzc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          shc=(gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1)))*
     &      (gl(kh3,kkci,kkcj)*gl(kz,kkcj,kkci)-
     &      gr(kh3,kkci,kkcj)*gr(kz,kkcj,kkci))
          sh=dreal(shc)
          if (abs(dimag(shc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dhc=(gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1)))*
     &      (gl(kh3,kkci,kkcj)*gr(kz,kkcj,kkci)-
     &      gr(kh3,kkci,kkcj)*gl(kz,kkcj,kkci))
          dh=dreal(dhc)
          if (abs(dimag(dhc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
c  definition of d, signmj, roota, rootad, rootd
          dq=(mcj/mci)**2
          signmj=mcj/(dabs(mcj))
          roota=signm*signmi*dsqrt(dabs(aq))
          rootad=signm*signmj*dsqrt(dabs(aq*dq))
          rootd=signmi*signmj*dsqrt(dabs(dq))
c  non diagonal coefficients
          f4=((2.d0*roota*dbh+rootd*sch+sbh)
     &      /(0.5d0+dq/2.d0-aq-bq))/8.d0
          f5=((2.d0*rootad*dch+rootd*sch+dq*sbh)
     &      /(0.5d0+dq/2.d0-aq-bq))/8.d0
          f6=-((roota*dbh-rootad*dch)/(aq-cq/4.d0))/8.d0
     &      -(((1.d0-dq)*sbh)/(aq-cq/4.d0))/16.d0
          f7=(((0.5d0+dq/2.d0-bq-cq/4.d0)*sbh)/(aq-cq/4.d0)**2)
     &      *cq/32.d0
          f8=-(((0.5d0+dq/2.d0-cq/4.d0)*sbh)
     &      /(0.5d0+dq/2.d0-bq-cq/4.d0))/8.d0
          f9=((bq*sbh)/(0.5d0+dq/2.d0-bq-cq/4.d0))/8.d0
          f10temp=(4.d0*aq-mass(kh3)**2/mci**2)/
     &      ((4.d0*aq-mass(kh3)**2/mci**2)**2+
     &      mass(kh3)**2*h3wid**2/mci**4)
          f10temp2=(4.0d0*mlsp**2-mz**2)**2/
     &    ((4.0d0*mlsp**2-mz**2)**2+zwid**2*mz**2)
          f10=((rootd*dz+sz)/cq*f10temp2)/16.d0+
     &      (roota*sh*f10temp)/8.d0
          f11=((rootd*dz+dq*sz)/cq*f10temp2)/16.d0+
     &      ((rootad*dh)*f10temp)/8.d0
c  non diagonal contribution
c  box part:
          par1r=par1r+f4*dsi_14(aq,dq,1.d0,cq/4.d0)
     &      +f5*dsi_14(aq,1.d0,dq,cq/4.d0)
     &      +f4*dsi_23(aq,bq,cq/4.d0)
     &      +f5*dsi_24(aq,bq,dq,cq/4.d0)
     &      +(f5+f6+f7+f8)*dsti_214(aq,bq,dq,cq/4.d0)
     &      +(f4-f6+f7+f8)*dsti_224(aq,bq,dq,cq/4.d0)
     &      +f9*dsi_33(aq,bq,cq/4.d0)
     &      +f9*dsi_34(aq,bq,dq,cq/4.d0)
     &      -sbh/4.d0*dsi_42(aq,bq,dq,cq)
c  triangular part:
          par2r=par2r+f10*dsi_14(aq,dq,1.d0,cq/4.d0)
     &      +f11*dsi_14(aq,1.d0,dq,cq/4.d0)
 13     continue
        rehbx=rehbx+par1r
        rehtx=rehtx+par2r
 23   continue
c
cccccccccc  first contribution to regbx and imgbx  cccccccccccccccccc
c
      mw=mass(kw)
      do 24 i=1,2
        par1r=0.d0
        par1i=0.d0
        kkci=kcha(i)
        mci=mass(kkci)
c diagonal coupling
c bug fixed: khc for kw        !Piero Ullio 00-11-29
c        sawc=(dsabsq(gl(kw,kn(1),kkci))+dsabsq(gr(khc,kn(1),kkci)))
c     &    *g2weak*costhw/2.d0
        sawc=(dsabsq(gl(kw,kn(1),kkci))+dsabsq(gr(kw,kn(1),kkci)))
     &    *g2weak*costhw/2.d0
        saw=dreal(sawc)
        if (abs(dimag(sawc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
        dawc=dreal(gl(kw,kn(1),kkci)*conjg(gr(kw,kn(1),kkci)))
     &    *g2weak*costhw
        daw=dreal(dawc)
        if (abs(dimag(dawc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
c definition of a,b,c,signmi,roota
        aq=(mlsp/mci)**2
        bq=(mw/mci)**2
        cq=(mz/mci)**2
        signmi=mci/(dabs(mci))
        roota=signm*signmi*dsqrt(dabs(aq))
c diagonal coefficients
        g1=2.d0*(((aq-bq)*saw)/(1.d0+aq-bq))
        g2=-((3*saw-4*roota*daw)/(1.d0-bq+cq/4.d0))/2.d0
        g3=-(((1.d0-aq+bq)*saw)/(1.d0+aq-bq))
        g4=(((2.d0+bq-cq/4.d0)*saw-4.d0*roota*daw)
     &    /(1.d0-bq+cq/4.d0))/2.d0
        g9=-(((1.d0-bq+cq/4.d0)*saw)/(aq-cq/4.d0)**2)*cq/8.d0
c diagonal contributions
        par1r=par1r+(g1*dsi_13(aq,bq,cq/4.d0)+
     &    g2*dsi_23(aq,bq,cq/4.d0)+g3*dsi_33(aq,bq,cq/4.d0)+
     &    (g4+g3+g9)*dsti_33(aq,bq,cq/4.d0)
     &    -saw/2.d0*dsi_41(aq,bq,cq))
        if(mw.ge.dabs(mlsp)) goto 54
        par1i=par1i-g1*dsj_1(aq,bq)
 54     continue
        do 14 j=1,2
          kkcj=kcha(j)
          mcj=mass(kkcj)
c non-diagonal couplings
          sbwc=(gl(kw,kn(1),kkcj)*conjg(gl(kw,kn(1),kkci))+
     &      gr(kw,kn(1),kkcj)*conjg(gr(kw,kn(1),kkci)))*
     &      (gl(kz,kkcj,kkci)+gr(kz,kkcj,kkci))/4.d0
          sbw=dreal(sbwc)
          if (abs(dimag(sbwc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dbwc=(gl(kw,kn(1),kkcj)*conjg(gr(kw,kn(1),kkci))+
     &      gr(kw,kn(1),kkcj)*conjg(gl(kw,kn(1),kkci)))*
     &      (gl(kz,kkcj,kkci)+gr(kz,kkcj,kkci))/4.d0
          dbw=dreal(dbwc)
          if (abs(dimag(dbwc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
c definition of d and signmj
          dq=(mcj/mci)**2
          signmj=mcj/(dabs(mcj))
          rootad=signm*signmj*dsqrt(dabs(aq*dq))
          rootd=signmi*signmj*dsqrt(dabs(dq))
c non diagonal coefficients
          g5=((2.d0*roota*dbw-(rootd/2.d0+0.5d0)*sbw)
     &      /(0.5d0+dq/2.d0-aq-bq))/2.d0
          g6=((2.d0*rootad*dbw-(rootd/2.d0+dq/2.d0)*sbw)
     &      /(0.5d0+dq/2.d0-aq-bq))/2.d0
          g7=-(((rootad+roota)*dbw
     &      -(rootd+dq/4.d0+0.25d0+cq/8.d0)*sbw)
     &      /(0.5d0+dq/2.d0-bq-cq/4.d0))/2.d0
          g8=(((rootad+roota)*dbw-(rootd+bq/2.d0+cq/4.d0)*sbw)
     &      /(0.5d0+dq/2.d0-bq-cq/4.d0))/2.d0
          g10=(((0.5d0+dq/2.d0-bq-cq/4.d0)*sbw)
     &      /(aq-cq/4.d0)**2)*cq/16.d0
          g11=-(((1.d0-dq)*sbw)/(aq-cq/4.d0))/8.d0
c contribution to the real part
          par1r=par1r+g5*dsi_14(aq,dq,1.d0,cq/4.d0)+
     &      g6*dsi_14(aq,1.d0,dq,cq/4.d0)+
     &      g5*dsi_23(aq,bq,cq/4.d0)+g6*dsi_24(aq,bq,dq,cq/4.d0)+
     &      (g6+g7+g10+g11)*dsti_214(aq,bq,dq,cq/4.d0)+
     &      (g5+g7+g10-g11)*dsti_224(aq,bq,dq,cq/4.d0)+
     &      g8*dsi_33(aq,bq,cq/4.d0)+
     &      g8*dsi_34(aq,bq,dq,cq/4.d0)
     &      -sbw/2.d0*dsi_42(aq,bq,dq,cq)
 14     continue
        regbx=regbx+par1r
        imgbx=imgbx+par1i
 24   continue
c
cccccccccc  second contribution to regbx and imgbx  ccccccccccccccccc
c
      do 25 i=1,2
        par1r=0.d0
        kkci=kcha(i)
        mci=mass(kkci)
c diagonal coupling
        sagc=gl(kz,khc,khc)*(dsabsq(gl(kgoldc,kn(1),kkci))+
     &     dsabsq(gr(kgoldc,kn(1),kkci)))/2.d0
        sag=dreal(sagc)
        if (abs(dimag(sagc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
        sdgc=g2weak*(gl(kw,kn(1),kkci)*conjg(gl(kgoldc,kn(1),kkci))+
     &    gr(kw,kn(1),kkci)*conjg(gr(kgoldc,kn(1),kkci)))/2.d0
        sdg=dreal(sdgc)
        if (abs(dimag(sdgc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
        ddgc=g2weak*(conjg(gr(kw,kn(1),kkci))*gl(kgoldc,kn(1),kkci)+
     &    conjg(gl(kw,kn(1),kkci))*gr(kgoldc,kn(1),kkci))/2.d0
        ddg=dreal(ddgc)
        if (abs(dimag(ddgc)).ge.1.d-10) then
          write(*,*) 'complex couplings!!!, model no. ',idtag
          stop
        endif
c definition of a,b,c
        aq=(mlsp/mci)**2
        bq=(mw/mci)**2
        cq=(mz/mci)**2
        signmi=mci/(dabs(mci))
        roota=signm*signmi*dsqrt(dabs(aq))
        rootc=signmi*dsqrt(dabs(cq))
c diagonal coefficients
        h1=-(sag/(1.d0-bq+cq/4.d0))/4.d0
        h2=((sag*(bq-cq/4.d0))/(1.d0-bq+cq/4.d0))/4.d0
        h3=((sag*(1.d0-bq+cq/4.d0))/(aq-cq/4.d0)**2)*cq/16.d0
        h10=-((sdg+ddg)/(1.d0-bq+cq/4.d0))*rootc/4.d0
c diagonal contributions
        par1r=par1r+((h1+h10)*dsi_23(aq,bq,cq/4.d0)+
     &    (h2+h3-h10)*dsti_33(aq,bq,cq/4.d0)+sag/4.d0*dsi_41(aq,bq,cq))
        do 15 j=1,2
          kkcj=kcha(j)
          mcj=mass(kkcj)
c non-diagonal couplings
          sbgc=(gl(kgoldc,kn(1),kkcj)*conjg(gl(kgoldc,kn(1),kkci))*
     &      gl(kz,kkcj,kkci)+gr(kgoldc,kn(1),kkcj)*
     &      conjg(gr(kgoldc,kn(1),kkci))*gr(kz,kkcj,kkci))/2.d0
          sbg=dreal(sbgc)
          if (abs(dimag(sbgc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dbgc=(gl(kgoldc,kn(1),kkcj)*conjg(gr(kgoldc,kn(1),kkci))*
     &      gl(kz,kkcj,kkci)+gr(kgoldc,kn(1),kkcj)*
     &      conjg(gl(kgoldc,kn(1),kkci))*gr(kz,kkcj,kkci))/2.d0
          dbg=dreal(dbgc)
          if (abs(dimag(dbgc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          scgc=(gl(kgoldc,kn(1),kkcj)*conjg(gl(kgoldc,kn(1),kkci))*
     &      gr(kz,kkcj,kkci)+gr(kgoldc,kn(1),kkcj)*
     &      conjg(gr(kgoldc,kn(1),kkci))*gl(kz,kkcj,kkci))/2.d0
          scg=dreal(scgc)
          if (abs(dimag(scgc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
          dcgc=(gl(kgoldc,kn(1),kkcj)*conjg(gr(kgoldc,kn(1),kkci))*
     &      gr(kz,kkcj,kkci)+gr(kgoldc,kn(1),kkcj)*
     &      conjg(gl(kgoldc,kn(1),kkci))*gl(kz,kkcj,kkci))/2.d0
          dcg=dreal(dcgc)
          if (abs(dimag(dcgc)).ge.1.d-10) then
            write(*,*) 'complex couplings!!!, model no. ',idtag
            stop
          endif
c definition of d and signmj
          dq=(mcj/mci)**2
          signmj=mcj/(dabs(mcj))
          rootad=signm*signmj*dsqrt(dabs(aq*dq))
          rootd=signmi*signmj*dsqrt(dabs(dq))
c non diagonal coefficients
          h4=((2.d0*roota*dbg+rootd*scg+sbg)
     &      /(0.5d0+dq/2.d0-aq-bq))/8.d0
          h5=((2.d0*rootad*dcg+rootd*scg+dq*sbg)
     &      /(0.5d0+dq/2.d0-aq-bq))/8.d0
          h6=-((roota*dbg-rootad*dcg)/(aq-cq/4.d0))/8.d0
     &      -(((1.d0-dq)*sbg)/(aq-cq/4.d0))/16.d0
          h7=(((0.5d0+dq/2.d0-bq-cq/4.d0)*sbg)/(aq-cq/4.d0)**2)
     &      *cq/32.d0
          h8=-(((0.5d0+dq/2.d0-cq/4.d0)*sbg)
     &      /(0.5d0+dq/2.d0-bq-cq/4.d0))/8.d0
          h9=((bq*sbg)/(0.5d0+dq/2.d0-bq-cq/4.d0))/8.d0
c contribution to the real part
          par1r=par1r+((h4)*dsi_14(aq,dq,1.d0,cq/4.d0)+
     &      h5*dsi_14(aq,1.d0,dq,cq/4.d0)+h4*dsi_23(aq,bq,cq/4.d0)+
     &      h5*dsi_24(aq,bq,dq,cq/4.d0)+
     &      (h5+h6+h7+h8)*dsti_214(aq,bq,dq,cq/4.d0)+
     &      (h4-h6+h7+h8)*dsti_224(aq,bq,dq,cq/4.d0)+
     &      h9*dsi_33(aq,bq,cq/4.d0)+
     &      h9*dsi_34(aq,bq,dq,cq/4.d0)-sbg/4.d0*dsi_42(aq,bq,dq,cq))
 15     continue
        regbx=regbx+par1r
 25   continue
c
cccccccccc  total and partial results  cccccccccccccccccccccccccccccc
c
      reres=refbx+reftx+rehbx+rehtx+regbx
      refbx=refbx/reres
      reftx=reftx/reres
      rehbx=rehbx/reres
      rehtx=rehtx/reres
      regbx=regbx/reres
      imres=imfbx+imftx+imgbx
      imfbx=imfbx/imres
      imftx=imftx/imres
      imgbx=imgbx/imres
      imres=imres*pi
      return
      end




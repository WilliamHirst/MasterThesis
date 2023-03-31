************************************************************************
*** SUBROUTINE dsasdwdcossfsf                                        ***
*** computes dW_{ij}/dcostheta                                       ***
*** for sfermion(1) + antisfermion(2) plus sfermion(1) + sfermion(2) ***
***                                                                  *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                              ***
*** Date: 01-08-10                                                   ***
*** modified: Joakim Edsjo, Mia Schelke to include squarks with      ***
***   gauge and Higgs boson final states and gluons                  ***
***   02-05-22                                                       ***
***   bug with switching of initial states fixed 020613 (edsjo)      ***
*** modified: Piero Ullio                                            ***
***   02-03-22                                                       ***
*** modified: Piero Ullio                                            ***
***   02-07-01                                                       ***
*** modified: Mia Schelke to include squark mixing:                  ***
***   diff family i.s. sq's with bosons in final state               ***
***   Jan 2007                                                       ***
*** modified: Piero Ullio to include mixing in fermion final states  ***
***   and a new labelling of states                                  ***
***   08-05-30                                                       ***
************************************************************************

      real*8 function dsasdwdcossfsf(p,costhe,kp1,kp2)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      real*8 p,costhe
      integer kp1,kp2,kp3,kp4,itmp,icase
      integer kp1i,kp2i
      integer f,f1,f2,i,kfer,kfer1,kfer2,fsave,f2save
      real*8 w,result,par
      logical wok

************************************************************************
      if (p.lt.0.0d0) then
         dsasdwdcossfsf=0.0d0
      return
      endif

***** tollerance on the position of a threshold and on a negative
***** fermion channel
      thstep=1.d-15
      fertoll=1.d-10

c...Since we want to be able to switch kp1 and kp2 without having
c...them changed on return, let's use some temporary variables
      kp1i=kp1
      kp2i=kp2

************************************************************************
*****
***** set the askin variables
      p12=p
      costheta=costhe

 100  continue
***** family indices:
      iifam(1)=ivfam(kp1i)
      iifam(2)=ivfam(kp2i)

***** decide on type for kp1i and kp2i, i.e.:
*****   itype(iii)=ivfam(ku) for up-type (s)quark
*****   itype(iii)=ivfam(kd) for down-type (s)quark
*****   itype(iii)=ivfam(iii) for (s)leptons
      itype(1)=ivtype(kp1i)
      itype(2)=ivtype(kp2i)

***** NOTE: boson final states are computed just for kp1i up-type and 
***** kp2i down-type, not viceversa, so in the latter case switch them
***** reordering this squark first if it appears in case 3
      if(itype(1).eq.(itype(2)+1)) then
        itmp=kp1i
        kp1i=kp2i
        kp2i=itmp
        goto 100
      endif 
***** reordering this squark first if it appears in case 3
      if(itype(2).ge.ivtype(ku).and.itype(1).lt.ivtype(ku)) then
        itmp=kp1i
        kp1i=kp2i
        kp2i=itmp
        goto 100
      endif 

***** masses in initial state:
      mass1=mass(kp1i)
      mass2=mass(kp2i)

***** kinematic variables:
      call dsaskinset1
***** 
*****    3 cases:
***** 1) sfermions with same type index, i.e.:
*****    slepton + slepton* (same family, both charged or both neutral) 
*****    up-type-squark1 + up-type-squark2  or 
*****    down-type-squark1 + down-type-squark2 
      if(itype(1).eq.itype(2)) then
        do i=1,54
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then
        kp1s=kp1i
        kp2s=kp2i
ccc
ccc do the resetting in the various cases:
ccc
        do i=1,30
          chon(i)=.false.
        enddo
        gluonin=.false.
        gammain=.false.
        neutcurr=.false.
        nosneutrinov=.true.
ccc
        if(itype(1).ge.ivtype(knue).and.itype(1).le.ivtype(ktau)) then
ccc input:   \tilde{l}_1,2(i) \tilde{l}^*_1,2(i)
ccc
ccc fermion channels:
ccc
        i=0
***** f fbar
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
        enddo
***** no up-type quark(i), up-type antiquark(j) i.ne.j 
***** no down-type quark(i), down-type antiquark(j) i.ne.j 
***** lepton1,  lepton2
        i=25
        chon(i)=.true.
        do f=knue,ktau
          if(itype(1).eq.ivtype(f)) then
            kp3in(i)=f
            kp4in(i)=f
            fsave=f
          endif
        enddo      
ccc
ccc gauge boson channels:
ccc
        if(dabs(echarg(fsave)).gt.1.d-15) gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        if(dabs(echarg(fsave)).gt.1.d-15) then ! for charged sleptons:
          ksfertc(1)=ksnu(fsave/2)
          nsfertc=1
          nsferuc=0
          ksfertn(1)=ksl(fsave/2)
          ksfertn(2)=ksl(fsave/2+3)
          nsfertn=2
          ksferun(1)=ksl(fsave/2)
          ksferun(2)=ksl(fsave/2+3)
          nsferun=2
        else ! for sneutrinos:
          nsfertc=0
          ksferuc(1)=ksl((fsave+1)/2)
          ksferuc(2)=ksl((fsave+1)/2+3)
          nsferuc=2
          ksfertn(1)=ksnu((fsave+1)/2)
          nsfertn=1
          ksferun(1)=ksnu((fsave+1)/2)
          nsferun=1
          nosneutrinov=.false.
        endif
ccc
        gg1=1.d0 
        gg2=1.d0
        chcol=1
        cfactini=1.d0
        cfactfin=3.d0
        endif
ccc
        if(iifam(1).eq.iifam(2).and.itype(1).eq.ivtype(ku)) then
ccc input:   \tilde{u}_1,2(i) \tilde{u}^*_1,2(i)
ccc
ccc fermion channels:
ccc
        i=0
***** f fbar
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
        enddo
***** no up-type quark(i), up-type antiquark(j) i.ne.j 
        i=i+6
***** down-type quark(i), down-type antiquark(j) i.ne.j 
        do f1=1,3
        do f2=1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo
***** up-type quark1, up-type quark2 
        do f1=1,3
        do f2=f1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          i=i+1  
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or.
     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
            fsave=f1
          endif
        enddo
        enddo
ccc
ccc boson channels:
ccc
        gluonin=.true.
        gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        nsfertc=0 
        do i=1,6
          ksferuc(i)=ksqd(i)
        enddo
        nsferuc=6
        ksfertn(1)=ksqu(fsave)
        ksfertn(2)=ksqu(fsave+3)
        nsfertn=2
        ksferun(1)=ksqu(fsave)
        ksferun(2)=ksqu(fsave+3)
        nsferun=2
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
        if(iifam(1).eq.iifam(2).and.itype(1).eq.ivtype(kd)) then
ccc input:   \tilde{d}_1,2(i) \tilde{d}^*_1,2(i)
        i=0
***** f fbar
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
        enddo
***** up-type quark(i), up-type antiquark(j) i.ne.j 
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo
***** no down-type quark(i), down-type antiquark(j) i.ne.j 
        i=i+6
***** down-type quark1, down-type quark2 
        do f1=1,3
        do f2=f1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          i=i+1  
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or.
     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
            fsave=f1
          endif
        enddo
        enddo  
ccc
ccc boson channels:
ccc
        gluonin=.true.
        gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        do i=1,6
          ksfertc(i)=ksqu(i)
        enddo
        nsfertc=6
        nsferuc=0
        ksfertn(1)=ksqd(fsave)
        ksfertn(2)=ksqd(fsave+3)
        nsfertn=2
        ksferun(1)=ksqd(fsave)
        ksferun(2)=ksqd(fsave+3)
        nsferun=2
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
        if(iifam(1).ne.iifam(2).and.itype(1).eq.ivtype(ku)) then
ccc input:   \tilde{u}_1,2(i) \tilde{u}^*_1,2(j)  i.ne.j
ccc
ccc fermion channels:
ccc
        i=0
***** down-type quark(i), down-type antiquark(i) 
        do f=1,3
          i=kd+2*(f-1)
          kp3in(i)=i
          kp4in(i)=i
          chon(i)=.true.
        enddo
***** up-type quark(i), up-type antiquark(j) i,j matching initial state
        i=12
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4)) then
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif
          endif         
        enddo
        enddo
***** down-type quark(i), down-type antiquark(j) i.ne.j 
        do f1=1,3
        do f2=1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo
***** up-type quark1, up-type quark2 
        do f1=1,3
        do f2=f1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          i=i+1  
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or.
     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
          endif
        enddo
        enddo
ccc
ccc only w gauge boson channels:
ccc
ccc
ccc t- and u-channels:
ccc
        nsfertc=0 
        do i=1,6
          ksferuc(i)=ksqd(i)
        enddo
        nsferuc=6
        nsfertn=0
        nsferun=0
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
        if(iifam(1).ne.iifam(2).and.itype(1).eq.ivtype(kd)) then
ccc input:   \tilde{d}_1,2(i) \tilde{d}^*_1,2(j)  i.ne.j
ccc
ccc fermion channels:
ccc
        i=0
***** up-type quark(i), up-type antiquark(i) 
        do f=1,3
          i=ku+2*(f-1)
          kp3in(i)=i
          kp4in(i)=i
          chon(i)=.true.
        enddo
***** up-type quark(i), up-type antiquark(j) i.ne.j
        i=12
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo
***** down-type quark(i), down-type antiquark(j) i,j matching initial state 
        do f1=1,3
        do f2=1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4)) then
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif
          endif         
        enddo
        enddo
***** down-type quark1, down-type quark2 
        do f1=1,3
        do f2=f1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          i=i+1  
          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or.
     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
          endif
        enddo
        enddo
ccc
ccc only w gauge boson channels:
ccc
ccc
ccc t- and u-channels:
ccc
        do i=1,6
          ksfertc(i)=ksqu(i)
        enddo
        nsfertc=6
        nsferuc=0
        nsfertn=0
        nsferun=0
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
        endif
*****
******************************************************* sf sf* -> f fbar
        do i=1,12
          if(chon(i)) then
          kp3=kp3in(i)
          kp4=kp4in(i)
          if(chcol.eq.1) then
            icase=1
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            if(ncolor(i).gt.2.d0) then
              result=result*cfactfin 
            endif
          elseif(chcol.eq.2) then
            icase=1  
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
          else
            write(*,*) 'chcol not set correctly, = ',chcol
            stop
          endif  
          prtial(i)=result
          endif
        enddo  
***************************************************** sf sf* -> q1 q2bar
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state, not in the final state:
*****
        if(chcol.eq.2) then  
        do i=13,24
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=1
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result
          endif
        enddo
        endif
************************************************* sl sl -> lepton lepton
        if(chcol.eq.1) then  
          i=25
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=1
            call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result
          endif
        endif 
********************************************************* sf sf -> q1 q2
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state, not in the final state:
*****
        if(chcol.eq.2) then  
        do i=25,30
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=1
            call dsasferecol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result
          endif
        enddo
        endif
******************************************************* sf sf* ->  w+ w-
        kp3=kw
        kp4=kw
        if(nsferuc.ne.0) then
          icase=2
        else
          icase=3
        endif
        call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(31)=par
******************************************************** sf sf* -> h- w+
        kp3=kw
        kp4=khc
        if(nsferuc.ne.0) then
          icase=5
        else
          icase=6
        endif
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(47)=2.d0*par    !factor of 2 to include sf sf* -> h+ w-
******************************************************** sf sf* -> h+ h-
        kp3=khc
        kp4=khc
        if(nsferuc.ne.0) then
          icase=6
        else
          icase=7
        endif
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(48)=par
****************************************************** sf sf* -> z gamma
        if(gammain) then 
          kp3=kz
          kp4=kgamma
          icase=2
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(33)=par
************************************************** sf sf* -> gamma gamma
          kp3=kgamma
          kp4=kgamma
          icase=1
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(34)=par
****************************************************** sf sf* -> gamma h
          kp3=kgamma
          kp4=kh2
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(38)=par
****************************************************** sf sf* -> gamma H
          kp3=kgamma
          kp4=kh1
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(39)=par
****************************************************** sf sf* -> gamma a
          kp3=kgamma
          kp4=kh3
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(40)=par
        endif  
********************************************************** sf sf* -> z z
        if(neutcurr) then
          kp3=kz
          kp4=kz
          icase=4
          call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(32)=par
********************************************************** sf sf* -> z h
          kp3=kz
          kp4=kh2
          icase=8
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(35)=par
********************************************************** sf sf* -> z H
          kp3=kz
          kp4=kh1
          icase=8
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(36)=par
********************************************************** sf sf* -> z a
          kp3=kz
          kp4=kh3
          icase=9
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(37)=par
********************************************************** sf sf* -> h h
          kp3=kh2
          kp4=kh2
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(41)=par
********************************************************** sf sf* -> h a
          kp3=kh2
          kp4=kh3
          icase=4
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(42)=par
********************************************************** sf sf* -> H a
          kp3=kh1
          kp4=kh3
          icase=4
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(43)=par
********************************************************** sf sf* -> a a
          kp3=kh3
          kp4=kh3
          icase=5
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(44)=par
********************************************************** sf sf* -> H h
          kp3=kh1
          kp4=kh2
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(45)=par
********************************************************** sf sf* -> H H
          kp3=kh1
          kp4=kh1
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(46)=par
        endif
************************************************** sf sf* -> gluon gluon
        if(gluonin) then
          kp3=kgluon
          kp4=kgluon
          icase=2
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(49)=par  
************************************************** sf sf* -> Z gluon
          kp3=kz
          kp4=kgluon
          icase=3
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(50)=par
************************************************** sf sf* -> gamma gluon
          kp3=kgamma
          kp4=kgluon
          icase=3
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(51)=par
************************************************** sf sf* -> gluon H1
          kp3=kgluon
          kp4=kh1
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(52)=par
************************************************** sf sf* -> gluon H2
          kp3=kgluon
          kp4=kh2
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(53)=par
************************************************** sf sf* -> gluon H3
          kp3=kgluon
          kp4=kh3
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(54)=par
        endif  
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sf sf* -> nue nuebar
        w=w+prtial(2)            ! sf sf* -> e+ e-
        w=w+prtial(3)            ! sf sf* -> numu numubar
        w=w+prtial(4)            ! sf sf* -> mu+ mu-
        w=w+prtial(5)            ! sf sf* -> nutau nutaubar
        w=w+prtial(6)            ! sf sf* -> tau+ tau-
        w=w+prtial(7)            ! sf sf* -> u ubar
        w=w+prtial(8)            ! sf sf* -> d dbar
        w=w+prtial(9)            ! sf sf* -> c cbar
        w=w+prtial(10)           ! sf sf* -> s sbar
        w=w+prtial(11)           ! sf sf* -> t tbar
        w=w+prtial(12)           ! sf sf* -> b bbar
        w=w+prtial(13)           ! sf sf* -> u cbar
        w=w+prtial(14)           ! sf sf* -> u tbar
        w=w+prtial(15)           ! sf sf* -> c ubar
        w=w+prtial(16)           ! sf sf* -> c tbar
        w=w+prtial(17)           ! sf sf* -> t ubar
        w=w+prtial(18)           ! sf sf* -> t cbar
        w=w+prtial(19)           ! sf sf* -> d sbar
        w=w+prtial(20)           ! sf sf* -> d bbar
        w=w+prtial(21)           ! sf sf* -> s dbar
        w=w+prtial(22)           ! sf sf* -> s bbar
        w=w+prtial(23)           ! sf sf* -> b dbar
        w=w+prtial(24)           ! sf sf* -> b sbar
        w=w+prtial(25)           ! sf sf -> lepton lepton or u u or d d
        w=w+prtial(26)           ! sf sf -> u c  or  d s 
        w=w+prtial(27)           ! sf sf -> u t  or  d b
        w=w+prtial(28)           ! sf sf -> c c  or  s s
        w=w+prtial(29)           ! sf sf -> c t  or  s b
        w=w+prtial(30)           ! sf sf -> t t  or  b b
        w=w+prtial(31)           ! sf sf* ->  w+ w-
        w=w+prtial(32)           ! sf sf* -> z z
        w=w+prtial(33)           ! sf sf* -> z gamma
        w=w+prtial(34)           ! sf sf* -> gamma gamma
        w=w+prtial(35)           ! sf sf* -> z h
        w=w+prtial(36)           ! sf sf* -> z H
        w=w+prtial(37)           ! sf sf* -> z a
        w=w+prtial(38)           ! sf sf* -> gamma h
        w=w+prtial(39)           ! sf sf* -> gamma H
        w=w+prtial(40)           ! sf sf* -> gamma a
        w=w+prtial(41)           ! sf sf* -> h h
        w=w+prtial(42)           ! sf sf* -> h a
        w=w+prtial(43)           ! sf sf* -> H a
        w=w+prtial(44)           ! sf sf* -> a a
        w=w+prtial(45)           ! sf sf* -> H h
        w=w+prtial(46)           ! sf sf* -> H H
        w=w+prtial(47)           ! sf sf* -> h- w+ and h+ w-
        w=w+prtial(48)           ! sf sf* -> h+ h-
        w=w+prtial(49)           ! sf sf* -> gluon gluon
        w=w+prtial(50)           ! sf sf* -> Z gluon
        w=w+prtial(51)           ! sf sf* -> gamma gluon
        w=w+prtial(52)           ! sf sf* -> gluon H
        w=w+prtial(53)           ! sf sf* -> gluon h
        w=w+prtial(54)           ! sf sf* -> gluon A

c
c check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,30
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo
c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,54
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
************************************************************************
*****
***** 2) sneutrino + slepton in same family or
*****    up-type squark + down-type squark
      elseif(abs(itype(1)-itype(2)).eq.1) then
        do i=1,33
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then
        kp1s=kp1i
        kp2s=kp2i
*****
***** check whether to include gluon final states
***** 
        gluonin=.false.
ccc
        if(itype(1).ge.ivtype(knue).and.itype(1).le.ivtype(ktau)) then
ccc input:   \tilde{nu}(i) \tilde{l}^*_1,2(i)
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- same family j \in (1,2,3)
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    neutrino(i) lepton^-(i) -- same family as in the initial state
        i=13
        do f=1,3
          kp3=knue+2*(f-1)
          kp4=ke+2*(f-1)
          if(ivtype(kp3).eq.itype(1)) then
            kp3in(i)=kp3
            kp4in(i)=kp4
            fsave=f
          endif
        enddo
ccc
ccc gauge boson channels:
ccc
ccc
ccc t- and u-channels:
ccc
        ksfertc(1)=ksl(fsave)
        ksfertc(2)=ksl(fsave+3)
        nsfertc=2
        ksferuc(1)=ksnu(fsave)
        nsferuc=1
        ksfertn(1)=ksnu(fsave)
        nsfertn=1
        ksferun(1)=ksl(fsave)
        ksferun(2)=ksl(fsave+3)
        nsferun=2
ccc
        gg1=1.d0 
        gg2=1.d0
        chcol=1
        cfactini=1.d0
        cfactfin=3.d0
        endif
ccc
        if(itype(1).ge.ivtype(ku).and.(iifam(2)-iifam(1)).eq.1) then
ccc input:   \tilde{u}_1,2(i) \tilde{d}^*_1,2(i)
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- same family j \in (1,2,3)
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    up-type quark(i) down-type antiquark(i) -- same family 
ccc       as in the initial state
        do f=1,3
          kp3=ku+2*(f-1) ! up-type quark
          if(ivfam(kp3).eq.iifam(1)) then
            fsave=f
            f2save=f
          endif
        enddo
ccc
ccc gauge boson channels:
ccc
        gluonin=.true.
ccc
ccc t- and u-channels:
ccc
        ksfertc(1)=ksqd(fsave)
        ksfertc(2)=ksqd(fsave+3)
        nsfertc=2
        ksferuc(1)=ksqu(fsave)
        ksferuc(2)=ksqu(fsave+3)
        nsferuc=2
        ksfertn(1)=ksqu(fsave)
        ksfertn(2)=ksqu(fsave+3)
        nsfertn=2
        ksferun(1)=ksqd(fsave)
        ksferun(2)=ksqd(fsave+3)
        nsferun=2
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
        endif
ccc
        if(itype(1).ge.ivtype(ku).and.(iifam(2)-iifam(1)).ne.1) then
ccc input:   \tilde{u}(i) \tilde{d}^*_1,2(j)  i.ne.j
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- same family j \in (1,2,3)
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    up-type quark(i) down-type antiquark(i) -- same family 
ccc       as in the initial state
        do f1=1,3
          kp3=ku+2*(f1-1) ! up-type quark
          if(ivfam(kp3).eq.iifam(1)) then
            fsave=f1
          endif
        enddo
        do f2=1,3
          kp4=kd+2*(f2-1) ! down-type quark
          if(ivfam(kp4).eq.iifam(2)) then
            f2save=f2
          endif
        enddo
ccc
ccc gauge boson channels:
ccc
        gluonin=.true.
ccc
ccc t- and u-channels:
ccc
        ksfertc(1)=ksqd(f2save)
        ksfertc(2)=ksqd(f2save+3)
        nsfertc=2
        ksferuc(1)=ksqu(fsave)
        ksferuc(2)=ksqu(fsave+3)
        nsferuc=2
        ksfertn(1)=ksqu(fsave)
        ksfertn(2)=ksqu(fsave+3)
        nsfertn=2
        ksferun(1)=ksqd(f2save)
        ksferun(2)=ksqd(f2save+3)
        nsferun=2
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
        endif
ccc
        endif
ccc

*************************************************** sfu sfd* -> fu fdbar
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state and in the final state; for sleptons flavour mixing for 
***** quarks in the final state:
*****
        if(chcol.eq.1) then
***** first the lepton final states:
          do f=1,3
            kp3=knue+2*(f-1)
            kp4=ke+2*(f-1)
            icase=2
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(f)=result
          enddo
***** then the quark final states:
          i=0
          do f1=1,3
          do f2=1,3
            kp3=ku+2*(f1-1) ! up-type quark
            kp4=kd+2*(f2-1) ! down-type quark
            icase=2
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            result=result*cfactfin
            i=i+1  
            prtial(3+i)=result
          enddo
          enddo
        else
***** first the lepton final states:
          do f=1,3
            kp3=knue+2*(f-1)
            kp4=ke+2*(f-1)
            icase=2
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(f)=result
          enddo
***** then the quark final states:
          i=0
          do f1=1,3
          do f2=1,3
            kp3=ku+2*(f1-1) ! up-type quark
            kp4=kd+2*(f2-1) ! down-type quark
            icase=2
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            i=i+1  
            prtial(3+i)=result
          enddo
          enddo
        endif
******************************************************* sfu sfd -> fu fd
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state and in the final state
*****
***** first the lepton final states:
        if(chcol.eq.1) then
          i=13
          kp3=kp3in(i)
          kp4=kp4in(i)
          icase=2
          call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
          prtial(13)=result
*****
***** then the quark final states:
        else
        i=0
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1) ! up-type quark
          kp4=kd+2*(f2-1) ! down-type quark
          icase=2
          call dsasferecol(kp1i,kp2i,kp3,kp4,icase,result)
          i=i+1  
          prtial(12+i)=result
        enddo
        enddo
        endif
******************************************************* sfu sfd* -> w+ z
        kp3=kw
        kp4=kz
        icase=1
        call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(22)=par
*************************************************** sfu sfd* -> w+ gamma
        kp3=kw
        kp4=kgamma
        icase=1
        call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(23)=par
******************************************************* sfu sfd* -> w+ h
        kp3=kw
        kp4=kh2
        icase=3
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(24)=par
******************************************************* sfu sfd* -> w+ H
        kp3=kw
        kp4=kh1
        icase=3
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(25)=par
******************************************************* sfu sfd* -> w+ a
        kp3=kw
        kp4=kh3
        icase=4
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(26)=par
*************************************************** sfu sfd* -> h+ gamma
        kp3=kgamma
        kp4=khc
        icase=2
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(27)=par
******************************************************* sfu sfd* -> h+ z
        kp3=kz
        kp4=khc
        icase=1
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(28)=par
******************************************************* sfu sfd* -> h+ h
        kp3=khc
        kp4=kh2
        icase=1
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(29)=par
******************************************************* sfu sfd* -> h+ H
        kp3=khc
        kp4=kh1
        icase=1
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(30)=par
******************************************************* sfu sfd* -> h+ a
        kp3=khc
        kp4=kh3
        icase=2
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(31)=par
        if(gluonin) then
******************************************************* sfu sfd* -> W+ g
          kp3=kw
          kp4=kgluon
          icase=4
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(32)=par
******************************************************* sfu sfd* -> g H+
          kp3=kgluon
          kp4=khc
          icase=11
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(33)=par
        endif
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sfu sfd* -> nue e+
        w=w+prtial(2)            ! sfu sfd* -> numu mu+
        w=w+prtial(3)            ! sfu sfd* -> nutau tau+
        w=w+prtial(4)            ! sfu sfd* -> u dbar
        w=w+prtial(5)            ! sfu sfd* -> u sbar
        w=w+prtial(6)            ! sfu sfd* -> u bbar
        w=w+prtial(7)            ! sfu sfd* -> c ubar
        w=w+prtial(8)            ! sfu sfd* -> c sbar
        w=w+prtial(9)            ! sfu sfd* -> c bbar
        w=w+prtial(10)           ! sfu sfd* -> t dbar
        w=w+prtial(11)           ! sfu sfd* -> t sbar
        w=w+prtial(12)           ! sfu sfd* -> t bbar
        w=w+prtial(13)           ! sfu sfd -> u d
        w=w+prtial(14)           ! sfu sfd -> u s
        w=w+prtial(15)           ! sfu sfd -> u b
        w=w+prtial(16)           ! sfu sfd -> c d
        w=w+prtial(17)           ! sfu sfd -> c s
        w=w+prtial(18)           ! sfu sfd -> c b
        w=w+prtial(19)           ! sfu sfd -> t d
        w=w+prtial(20)           ! sfu sfd -> t s
        w=w+prtial(21)           ! sfu sfd -> t b
        w=w+prtial(22)           ! sfu sfd* -> w+ z
        w=w+prtial(23)           ! sfu sfd* -> w+ gamma
        w=w+prtial(24)           ! sfu sfd* -> w+ h
        w=w+prtial(25)           ! sfu sfd* -> w+ H
        w=w+prtial(26)           ! sfu sfd* -> w+ a
        w=w+prtial(27)           ! sfu sfd* -> h+ gamma
        w=w+prtial(28)           ! sfu sfd* -> h+ z
        w=w+prtial(29)           ! sfu sfd* -> h+ h
        w=w+prtial(30)           ! sfu sfd* -> h+ H
        w=w+prtial(31)           ! sfu sfd* -> h+ a
        w=w+prtial(32)           ! sfu sfd* -> W+ g
        w=w+prtial(33)           ! sfu sfd* -> g H+
c
c check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,21
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo
c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,33
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
************************************************************************
*****
***** 3) 2 sleptons in different families or
*****    one squark + one slepton
      else   
c        write(*,*) 'Now in case 3...'
        do i=1,12
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then
        kp1s=kp1i
        kp2s=kp2i
        do i=1,12
          chon(i)=.false.
        enddo
***** identify particles in initial and final state
***** the first only can be a up-type squark:
        if(itype(1).eq.ivtype(ku)) then
          i=0
          do f=1,3
            kp3=ku+2*(f-1)
            if(iifam(1).eq.ivfam(kp3)) then
              i=i+1
              kp3in(i)=kp3
              chon(i)=.true.
              kp3in(i+3)=kp3
              chon(i+3)=.true.
            endif
          enddo
          i=6
          do f=1,3
            kp3=kd+2*(f-1)
            i=i+1
            kp3in(i)=kp3
            chon(i)=.true.
          enddo
          ick1=1
***** 1 squark and 1 slepton
          gg1=3.d0
          gg2=1.d0
          chcol=1
          cfactini=3.d0
          cfactfin=1.d0
***** ... or a down-type squark:
        elseif(itype(1).eq.ivtype(kd)) then
          i=0
          do f=1,3
            kp3=kd+2*(f-1)
            if(iifam(1).eq.ivfam(kp3)) then
              i=i+1
              kp3in(i)=kp3
              chon(i)=.true.
              kp3in(i+3)=kp3
              chon(i+3)=.true.
            endif
          enddo
          i=6
          do f=1,3
            kp3=ku+2*(f-1)
            i=i+1
            kp3in(i)=kp3
            chon(i)=.true.
          enddo
          ick1=2
***** 1 squark and 1 slepton
          gg1=3.d0
          gg2=1.d0
          chcol=1
          cfactini=3.d0
          cfactfin=1.d0
***** ... or it is slepton:
        else
          i=0
          do f=1,3
            if(itype(1).eq.ivtype(knue+2*(f-1))) then
              kp3=knue+2*(f-1)
              i=i+1
              kp3in(i)=kp3
              chon(i)=.true.
              kp3in(i+3)=kp3
              chon(i+3)=.true.
              kp3=ke+2*(f-1)
              kp3in(i+6)=kp3
              chon(i+6)=.true.
              ick1=1
            elseif(itype(1).eq.ivtype(ke+2*(f-1))) then
              kp3=ke+2*(f-1)
              i=i+1
              kp3in(i)=kp3
              chon(i)=.true.
              kp3in(i+3)=kp3
              chon(i+3)=.true.
              kp3=knue+2*(f-1)
              kp3in(i+6)=kp3
              chon(i+6)=.true.
              ick1=2
            endif
          enddo 
***** 2 sleptons
          gg1=1.d0 
          gg2=1.d0
          chcol=1
          cfactini=1.d0
          cfactfin=3.d0
        endif
***** identify the second slepton:
        do f=1,3
          if(itype(2).eq.ivtype(knue+2*(f-1))) then
            kf2=knue+2*(f-1)
            kf2o=ke+2*(f-1)
            ick2=1
          elseif(itype(2).eq.ivtype(ke+2*(f-1))) then
            kf2=ke+2*(f-1)
            kf2o=knue+2*(f-1)
            ick2=2
          endif
        enddo
ccc
        endif
*****
*************************************************** sf1 sf2* -> f1 f2bar
        do i=1,3
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kf2
            icase=3
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result*cfactini
          endif
        enddo  
******************************************************* sf1 sf2 -> f1 f2
        do i=4,6
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kf2
            icase=3
            call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result*cfactini
          endif
        enddo  
*************************************************** sf1 sf2* -> f3 f4bar
**************************************************** or sf1 sf2 -> f3 f4
*****  f1 and f3 in the same doublet,  f2 and f4 in the same doublet
        do i=6,9
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kf2o
            if(ick1.eq.ick2) then
              icase=3
              call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            else
              icase=3
              call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
            endif
            prtial(i)=result*cfactini
          endif  
        enddo
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sq1(sl1) sl2* -> q1_1(l1) l2bar
        w=w+prtial(2)            ! sq1 sl2* -> q1_2 l2bar
        w=w+prtial(3)            ! sq1 sl2* -> q1_3 l2bar
        w=w+prtial(4)            ! sq1(sl1) sl2 -> q1_1(l1) l2
        w=w+prtial(5)            ! sq1 sl2 -> q1_2 l2
        w=w+prtial(6)            ! sq1 sl2 -> q1_3 l2
        w=w+prtial(7)            ! sf1(sl1) sf2*(sf2) -> f3_1(l3) f2bar(f2) 
        w=w+prtial(8)            ! sf1 sf2*(sf2) -> f3_2 f2bar(f2) 
        w=w+prtial(9)            ! sf1 sf2*(sf2) -> f3_3 f2bar(f2) 
c
c check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,9
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo
c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,12
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
************************************************************************
***** if statement corresponding to the 3 cases closed
      endif

      if (w.lt.0.0d0) then
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: dsasdwdcossfsf called with kp1i=',kp1i,
     &    ' kp2i=',kp2i
        write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
        write(*,*) 'DS: negative w, program stopped!'
        stop
      endif  

      return
      end









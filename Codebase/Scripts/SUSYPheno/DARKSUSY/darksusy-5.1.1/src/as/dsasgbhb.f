**************************************************************
*** SUBROUTINE dsasgbhb                                    ***
*** computes dW_{ij}/dcostheta                             ***
*** up-sfermion(i) + down-antisfermion(j) ->               ***
***    gauge boson + Higgs boson                           ***
*** ampl2 obtained summing over physical polarizations     ***
***                                                        *** 
*** input askin variables: p12,costheta                    ***
***    itype(1),itype(2),mass1,mass2                       ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 02-06-13                                         ***
*** This routine has been compared with the routine of     ***
*** Mia Schelke and the agreement is perfect, except in    ***
*** the low-p limit (due to widths in propagators)         ***
*** as expected.                                           ***
*** added flavour changing charged exchange for W^-H^+:    ***
*** added by Mia Schelke 2006-06-07  + Jan 2007            ***
*** added flavour changing charged exchange for:           ***
*** W^+ H^0_1,2,3 :  H^+ gamma,Z,gluon                     ***
*** added by Mia Schelke Jan 2007                          ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      SUBROUTINE dsasgbhb(kp1,kp2,kp3,kp4,icase,result)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,icase
      integer kgb,khb
      integer khs,kgs
      integer i
      real*8 ampl2,result,dsabsq,colfact
      double complex dsasdepro,A1,A2,Aaux
      real*8 pi
      parameter (pi=3.141592653589793238d0)


c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        result=0.d0
        return
      endif      

***** check initial state is ok and set type  
c      if(abs(itype(1)-itype(2)).gt.1) then
c        write(*,*) 'DS: dsasgbhb called with wrong particles'   
c        write(*,*) 'DS: in the initial state : ',pname(kp1),pname(kp2)  
c        stop  
c      endif    

***** set and check particles in the final state:
      kgb=kp3
      khb=kp4
c      if(kgb.ne.kz.and.kgb.ne.kw.and.
c     &   kgb.ne.kgamma.and.kgb.ne.kgluon) then
c        write(*,*) 'DS: dsasgbhb called with kp3 = ',pname(kgb)   
c        write(*,*) 'DS: rather than a gauge boson'
c        stop
c      endif 
c      if(khb.ne.kh1.and.khb.ne.kh2.and.khb.ne.kh3.and.
c     &   khb.ne.khc) then
c        write(*,*) 'DS: dsasgbhb called with kp4 = ',pname(khb)   
c        write(*,*) 'DS: rather than a Higgs boson'
c        stop
c      endif 

      colfact=1.d0
      if(itype(1).eq.ivtype(ku).or.itype(1).eq.ivtype(kd))
     &  colfact=3.d0

***** non equal particle final state: 
      s34=1.d0
***** masses in final state:
      mass3=mass(kp3)
      mass4=mass(kp4)
***** define the kinematic variables  
      call dsaskinset2
*****
      A1=(0.d0,0.d0)
      A2=(0.d0,0.d0)

*****
*****
******************* the first case **************************
***** sneutrino_i + slepton_i^+ -> Z-boson + H^+
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> Z-boson + H^+  with i,j = 1, 2, ... 6
      if(icase.eq.1) then 
***** t- & u-channel, up to 6 squarks or 1,2 sleptons **********
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
***** h^+ in s-channel ****************************************
        khs=khc
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        goto 600
      endif
*****
*****
******************* the second case *************************
***** sneutrino_i + slepton_i^+ -> photon + H^+
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> photon + H^+  with i,j = 1, 2, ... 6
      if(icase.eq.2) then
***** t- & u-channel, 6 squarks or 1,2 sleptons ***************
        if(itype(1).ge.ivtype(ku)) then
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        endif
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
***** h^+ in s-channel ****************************************
        khs=khc
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        goto 600
      endif
*****
*****
******************* the third case **************************
***** sneutrino_i + slepton_i^+ -> W^+ + H^0_1/H^0_2
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> W^+ + H^0_1/H^0_2  with i,j = 1, 2, ... 6
      if(icase.eq.3) then
***** t- & u-channel, up to 6 squarks or 1,2 sleptons ********
        do i=1,nsfertc
          A1=A1-2.d0*dsasdepro(Tvar,ksfertc(i))*gl(khb,kp2,ksfertc(i))
     &          *gl(kgb,ksfertc(i),kp1)
        enddo
        do i=1,nsferuc
          A2=A2+2.d0*dsasdepro(Uvar,ksferuc(i))*gl(khb,ksferuc(i),kp1)
     &          *gl(kgb,kp2,ksferuc(i))
        enddo
***** h^+ in s-channel ****************************************
        if(kp1-kp2.ge.-3) then
        khs=khc
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        endif
***** W^+ in s-channel ****************************************
        kgs=kw
        Aaux=gl(kgs,kp2,kp1)*mass(kw)*gl(khb,kgb,kgs)
     &        *dsasdepro(Svar,kgs)
        A1=A1+(1.d0-(mass1**2-mass2**2)/mass(kgs)**2)*Aaux
        A2=A2+(-1.d0-(mass1**2-mass2**2)/mass(kgs)**2)*Aaux
        goto 600
      endif
*****
*****
******************* the forth case *************************
***** sneutrino_i + slepton_i^+ -> W^+ + H^0_3
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> W^+ + H^0_3  with i,j = 1, 2, ... 6
      if(icase.eq.4) then
***** t- & u-channel, up to 6 squarks or 1,2 sleptons **********
        do i=1,nsfertc
          A1=A1-2.d0*dsasdepro(Tvar,ksfertc(i))*gl(khb,kp2,ksfertc(i))
     &          *gl(kgb,ksfertc(i),kp1)
        enddo
        if(itype(1).ge.ivtype(ku)) then
        do i=1,nsferuc
          A2=A2+2.d0*dsasdepro(Uvar,ksferuc(i))*gl(khb,ksferuc(i),kp1)
     &          *gl(kgb,kp2,ksferuc(i))
        enddo
        endif
***** h^+ in s-channel ****************************************
        if(kp1-kp2.ge.-3) then
        khs=khc
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        endif
        goto 600
      endif
*****
*****
******************* the fifth case **************************
***** up-type-squark(i) + anti-up-type-squark(j) -> W^- + H^+ 
*****   with i,j = 1, 2, ... 6
***** sneutrino_i + sneutrino_i^* -> W^- + H^+ 
      if(icase.eq.5) then
***** u-channel, 6 squarks or 2 sleptons ********************** 
        do i=1,nsferuc
          A2=A2+2.d0*dsasdepro(Uvar,ksferuc(i))*gl(khb,kp1,ksferuc(i))
     &          *gl(kgb,ksferuc(i),kp2)       
        enddo   
        if(abs(kp1-kp2).le.1) then
***** h1 in s-channel  ****************************************
        khs=kh1
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
***** h2 in s-channel  ****************************************
        khs=kh2
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        if(itype(1).eq.ivtype(ku)) then
***** h3 in s-channel  ****************************************
          khs=kh3
          A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
          A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        endif
        endif
        goto 600
      endif
***** 
*****
******************* the sixth case **************************
***** down-type-squark_i + anti-down-type-squark_j -> W^- + H^+ 
*****   with i,j = 1, 2, ... 6
***** slepton_i^- + slepton_i^+ -> W^- + H^+ 
      if(icase.eq.6) then
***** t-channel, 6 squarks or 1 sleptons **********************
        do i=1,nsfertc
          A1=A1-2.d0*dsasdepro(Tvar,ksfertc(i))*gl(khb,kp2,ksfertc(i))
     &          *gl(kgb,ksfertc(i),kp1)
        enddo
        if(abs(kp1-kp2).le.1) then
***** h1 in s-channel  ****************************************
        khs=kh1
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
***** h2 in s-channel  ****************************************
        khs=kh2
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
***** h3 in s-channel  ****************************************
        khs=kh3
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &        *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &        *gl(kgb,khb,khs)
        endif
        goto 600
      endif
*****
*****
******************* the seventh case ************************
***** slepton_i^- + slepton_i^+ -> photon + H^0_1/H^0_2/H^0_3
***** up-type-squark_i + anti-up-type-squark_j  
*****     -> photon + H^0_1/H^0_2/H^0_3  with i,j = 1, 2, ... 6
***** down-type-squark_i + anti-down-type-squark_j 
*****     -> photon + H^0_1/H^0_2/H^0_3  with i,j = 1, 2, ... 6
      if(icase.eq.7) then
***** t- & u-channel, 6 squarks or 2 sleptons *****************
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
        goto 600
      endif
***** 
*****
******************* the eighth case *************************
***** slepton_i + slepton_i^* -> Z-boson + H^0_1/H^0_2
***** up-type-squark_i + anti-up-type-squark_j  
*****     -> Z-boson + H^0_1/H^0_2  with i,j = 1, 2, ... 6
***** down-type-squark_i + anti-down-type-squark_j 
*****     -> Z-boson + H^0_1/H^0_2  with i,j = 1, 2, ... 6
      if(icase.eq.8) then
***** t- & u-channel, 6 squarks or 2 sleptons *****************
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
        if(abs(kp1-kp2).le.1) then
        if(nosneutrinov) then
***** (H^0_3 can not couple to a pair of sneutrinos as far as there 
*****  are no right-handed sneutrinos)
***** h3 in s-channel  ****************************************
          khs=kh3
          A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &        *gl(kgb,khb,khs)
          A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &        *gl(kgb,khb,khs)
        endif
***** Z^0 in s-channel ****************************************
        kgs=kz
        Aaux=gl(kgs,kp2,kp1)*mass(kw)*gl(khb,kgb,kgs)
     &        *dsasdepro(Svar,kgs)
        A1=A1+(1.d0-(mass1**2-mass2**2)/mass(kgs)**2)*Aaux
        A2=A2+(-1.d0-(mass1**2-mass2**2)/mass(kgs)**2)*Aaux
        endif
        goto 600
      endif
*****
*****
******************* the ninth case **************************
***** slepton_i + slepton_i^* -> Z-boson + H^0_3
***** up-type-squark_i + anti-up-type-squark_j  
*****     -> Z-boson + H^0_3  with i,j = 1, 2, ... 6
***** down-type-squark_i + anti-down-type-squark_j 
*****     -> Z-boson + H^0_3  with i,j = 1, 2, ... 6
      if(icase.eq.9) then
***** t- & u-channel, 6 squarks or 2 sleptons *****************
        if(nosneutrinov) then
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
        endif
        if(abs(kp1-kp2).le.1) then
***** h1 in s-channel  ****************************************
        khs=kh1
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
***** h2 in s-channel  ****************************************
        khs=kh2
        A1=A1-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        A2=A2-2.d0*dsasdepro(Svar,khs)*gl(khs,kp2,kp1)
     &          *gl(kgb,khb,khs)
        endif
        goto 600
      endif
*****
*****
******************* the tenth case **************************
***** up-type-squark_i + anti-up-type-squark_j 
*****     -> gluon + H^0_{1,2,3}  with i,j = 1, 2, ... 6
***** down-type-squark_i + anti-down-type-squark_j 
*****     -> gluon + H^0_{1,2,3}  with i,j = 1, 2, ... 6
      if(icase.eq.10) then
***** t- & u-channel, 6 squarks *******************************
        colfact=4.d0
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
        goto 600
      endif
*****
*****
******************* the eleventh case ***********************
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> gluon + H^+  with i,j = 1, 2, ... 6
      if(icase.eq.11) then
***** t- & u-channel, up to 6 squarks **************************
        colfact=4.d0
        do i=1,nsfertn
          A1=A1-2.d0*dsasdepro(Tvar,ksfertn(i))*gl(khb,kp2,ksfertn(i))
     &          *gl(kgb,ksfertn(i),kp1)
        enddo
        do i=1,nsferun
          A2=A2+2.d0*dsasdepro(Uvar,ksferun(i))*gl(khb,ksferun(i),kp1)
     &          *gl(kgb,kp2,ksferun(i))
        enddo
        goto 600
      endif
*****
      write(*,*) 'DS: dsasgbhb called with wrong icase :',icase   
      write(*,*) 'DS: initial and final states : ',  
     &             pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop
 600  continue
***** amplitude squared  **************************************
***** ampl2 obtained summing over physical polarizations    ***
      ampl2=p12**2*(1.d0-costheta**2)*dsabsq(A1-A2)
      if(dabs(mass3).gt.1.d-15) then
        Aaux=k34/mass3*(ep1*A1+ep2*A2)
     &       -costheta*ek3/mass3*p12*(A1-A2)
        ampl2=ampl2+dsabsq(Aaux)
      endif  
*****
      result=ampl2*k34/(8.d0*pi*gg1*gg2*s34*dsqrt(Svar))*colfact
      return
      end

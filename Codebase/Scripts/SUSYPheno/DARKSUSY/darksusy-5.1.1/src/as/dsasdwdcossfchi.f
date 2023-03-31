************************************************************************
*** SUBROUTINE dsasdwdcossfchi                                       ***
*** computes dW_{ij}/dcostheta                                       ***
*** for sfermion(1) + neutralino(2) (or chargino(2))                 ***
*** plus sfermion(1) + neutralino(2) (or chargino(2))                ***
***                                                                  *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                              ***
*** Date: 01-11-04                                                   ***
*** Modified: Joakim Edsjo, Mia Schelke                              ***
***   to include gluon final states, 2002-03-21                      ***
*** Modified: Piero Ullio                                            ***
***   to switch to ampl2 with physical polarizations, 02-07-01       ***
*** Modified: Mia Schelke                                            ***
***   to allow for initial state squarks of different family than    ***
***   final state quarks in selected cases, June 2006, Jan 2007      ***
*** modified: Piero Ullio to include a new labelling of states,      ***
***   08-05-30                                                       ***
************************************************************************

      real*8 function dsasdwdcossfchi(p,costhe,kp1,kp2)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      real*8 p,costhe
      integer kp1,kp2,kp3,kp4,kp2int,i,f,fsave,icase
      real*8 w,par
      logical wok
************************************************************************
      if (p.lt.0.0d0) then
         dsasdwdcossfchi=0.0d0
      return
      endif

***** tollerance on the position of a threshold
      thstep=1.d-15

      if(kp2.ge.kn1.and.kp2.le.kn4) kp2int=kn1
      if(kp2.ge.kcha1.and.kp2.le.kcha2) kp2int=kcha1
************************************************************************
*****
***** set the askin variables
      p12=p
      costheta=costhe
***** masses in initial state:
      mass1=mass(kp1)
      mass2=mass(kp2)
***** kinematic variables:
      call dsaskinset1
***** family indices:
      iifam(1)=ivfam(kp1)
      iifam(2)=0
***** internal degrees of freedom of the initial particles:
***** note: here we do not include particle anti-particle degrees of
*****       freedom. Only spin and colour are included here.
*****       The particle anti-particle degrees of freedom are
*****       taken into account when w is summed below.
      gg1c=1.d0
      if(iifam(1).ge.ivfam(ksu1).and.iifam(1).le.ivfam(ksb2)) then
        gg1c=3.d0
      endif
      gg2c=2.d0  ! both for neutralinos and charginos since particle
                 ! antiparticle states are not included here.
***** 
*****    2 cases:
***** 1) sfermion + neutralino
      if(kp2int.eq.kn1) then
        do i=1,54
          prtial(i)=0.0d0
        enddo
        if(kp1.ne.kp1c.or.kp2int.ne.kp2c) then
        kp1c=kp1
        kp2c=kp2int
        cgammain=.false.
        cgluonin=.false.
        ciaux=iifam(1)/10
        ciaux=iifam(1)-ciaux*10
        do f=knue,kb
          if(iifam(1).eq.ivfam(f)) then 
            fsave=f
          endif
        enddo
        kcfers=fsave
        if(kcfers.le.ktau.and.kcfers.ge.knue) then
          if(ciaux.eq.1) then
            ncsfert=1
            kcsfertn(1)=ksnu((fsave+1)/2)
            ncferd=1
            kcferd(1)=fsave+1
          endif
          if(ciaux.eq.2) then
            ncsfert=2
            kcsfertn(1)=ksl(fsave/2)
            kcsfertn(2)=ksl(fsave/2+3)
            ncferd=1
            kcferd(1)=fsave-1
            cgammain=.true.
          endif
        endif
        if(kcfers.le.kb.and.kcfers.ge.ku) then
          if(ciaux.eq.1) then
            ncsfert=2
            kcsfertn(1)=ksqu((fsave-6+1)/2)
            kcsfertn(2)=ksqu((fsave-6+1)/2+3)
            ncferd=3
            kcferd(1)=kd
            kcferd(2)=ks
            kcferd(3)=kb
            cgammain=.true.
            cgluonin=.true.
          endif
          if(ciaux.eq.2) then
            ncsfert=2
            kcsfertn(1)=ksqd((fsave-6)/2)
            kcsfertn(2)=ksqd((fsave-6)/2+3)
            ncferd=3
            kcferd(1)=ku
            kcferd(2)=kc
            kcferd(3)=kt
            cgammain=.true.
            cgluonin=.true.
          endif
        endif
********
        endif
******************************************************** sf chi^0 -> Z f
        kp3=kz
        kp4=kcfers
        icase=1
        call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
        prtial(1)=par
**************************************************** sf chi^0 -> gamma f
        if(cgammain) then
          kp3=kgamma
          kp4=kcfers
          icase=7
          call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
          prtial(2)=par
        endif
******************************************************** sf chi^0 -> h f
        kp3=kh2
        kp4=kcfers
        icase=1
        call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
        prtial(3)=par
******************************************************** sf chi^0 -> H f
        kp3=kh1
        kp4=kcfers
        icase=1
        call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
        prtial(4)=par
******************************************************** sf chi^0 -> a f
        kp3=kh3
        kp4=kcfers
        icase=2
        call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
        prtial(5)=par
**************************************************** sf chi^0 -> gluon f
        if(cgluonin) then
          kp3=kgluon
          kp4=kcfers
          icase=9
          call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
          prtial(8)=par
        endif
********
        if(ciaux.eq.1) then
          prtial(6)=0.d0
          prtial(7)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            if(kp4.le.ktau.and.kp4.ge.knue) then
              ncsfertc=2
              kcsfertc(1)=ksl(kp4/2)
              kcsfertc(2)=ksl(kp4/2+3)
            else
              ncsfertc=2
              kcsfertc(1)=ksqd((kp4-6)/2)
              kcsfertc(2)=ksqd((kp4-6)/2+3)
            endif
*********************************************** sfu^* chi^0 -> W^- fdbar
            kp3=kw
            icase=2
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par   
*********************************************** sfu^* chi^0 -> h^- fdbar
            kp3=khc
            icase=3
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
        elseif(ciaux.eq.2) then  
********
          prtial(6)=0.d0
          prtial(7)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            if(kp4.le.ktau.and.kp4.ge.knue) then
              ncsfertc=1
              kcsfertc(1)=ksnu((kp4+1)/2)
            else
              ncsfertc=2
              kcsfertc(1)=ksqu((kp4-6+1)/2)
              kcsfertc(2)=ksqu((kp4-6+1)/2+3)
            endif
**************************************************** sfd chi^0 -> W^- fu
            kp3=kw
            icase=2
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par
**************************************************** sfd chi^0 -> h^- fu
            kp3=khc
            icase=4
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
        else
          write(*,*) 'DS: wrong ciaux value =',ciaux
          write(*,*) 'DS: in dsasdwdcossfchi. program stopped.'
          stop
        endif
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sf chi^0 -> Z f
        w=w+prtial(2)            ! sf chi^0 -> gamma f
        w=w+prtial(3)            ! sf chi^0 -> h f
        w=w+prtial(4)            ! sf chi^0 -> H f
        w=w+prtial(5)            ! sf chi^0 -> a f
        w=w+prtial(6)            ! sfu^* chi^0 -> W^- fdbar
                                 ! or  sfd chi^0 -> W^- + fu
        w=w+prtial(7)            ! sfu^* chi^0 -> h^- fdbar
                                 ! or  sfd chi^0 -> h^- fu
        w=w+prtial(8)            ! sf chi^0 -> gluon f
        dsasdwdcossfchi = w      ! The weighting factor is 1

************************************************************************
*****
***** 2) sfermion + chargino
      elseif(kp2int.eq.kcha1) then
        do i=1,54
          prtial(i)=0.0d0
        enddo
        if(kp1.ne.kp1c.or.kp2int.ne.kp2c) then
        kp1c=kp1
        kp2c=kp2int
        cgluonin=.false.
        ciaux=iifam(1)/10
        ciaux=iifam(1)-ciaux*10
        do f=knue,kb
          if(iifam(1).eq.ivfam(f)) then 
            fsave=f
          endif
        enddo
        kcfers=fsave
        if(kcfers.le.ktau.and.kcfers.ge.knue) then
          if(ciaux.eq.1) then
            ncsfert=1
            kcsfertn(1)=ksnu((fsave+1)/2)
            ncfers=1
            kcfersv(1)=kcfers
            ncferd=1
            kcferd(1)=fsave+1
            ncsfertc=2
            kcsfertc(1)=ksl((fsave+1)/2)
            kcsfertc(2)=ksl((fsave+1)/2+3)
          endif
          if(ciaux.eq.2) then
            ncsfert=2
            kcsfertn(1)=ksl(fsave/2)
            kcsfertn(2)=ksl(fsave/2+3)
            ncfers=1
            kcfersv(1)=kcfers
            ncferd=1
            kcferd(1)=fsave-1
            ncsfertc=1
            kcsfertc(1)=ksnu(fsave/2)
          endif
        endif
        if(kcfers.le.kb.and.kcfers.ge.ku) then
          if(ciaux.eq.1) then
            ncsfert=2
            kcsfertn(1)=ksqu((fsave-6+1)/2)
            kcsfertn(2)=ksqu((fsave-6+1)/2+3)
            ncfers=3
            kcfersv(1)=ku
            kcfersv(2)=kc
            kcfersv(3)=kt
            ncferd=3
            kcferd(1)=kd
            kcferd(2)=ks
            kcferd(3)=kb
            ncsfertc=6
            do i=1,6
              kcsfertc(i)=ksqd(i)
            enddo
            cgluonin=.true.
          endif
          if(ciaux.eq.2) then
            ncsfert=2
            kcsfertn(1)=ksqd((fsave-6)/2)
            kcsfertn(2)=ksqd((fsave-6)/2+3)
            ncfers=3
            kcfersv(1)=kd
            kcfersv(2)=ks
            kcfersv(3)=kb
            ncferd=3
            kcferd(1)=ku
            kcferd(2)=kc
            kcferd(3)=kt
            ncsfertc=6
            do i=1,6
              kcsfertc(i)=ksqu(i)
            enddo
            cgluonin=.true.
          endif
        endif
********
        endif
********
        if(ciaux.eq.1) then
************************************************* sfu^* chi^+ -> Z fdbar
          kp3=kz
          prtial(1)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(1)=prtial(1)+par
          enddo
********************************************* sfu^* chi^+ -> gamma fdbar
          kp3=kgamma
          prtial(2)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=5
            call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
            prtial(2)=prtial(2)+par
          enddo
************************************************* sfu^* chi^+ -> h fdbar
          kp3=kh2
          prtial(3)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(3)=prtial(3)+par
          enddo
************************************************* sfu^* chi^+ -> H fdbar
          kp3=kh1
          prtial(4)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(4)=prtial(4)+par
          enddo
************************************************* sfu^* chi^+ -> a fdbar
          kp3=kh3
          prtial(5)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=2
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(5)=prtial(5)+par
          enddo
**************************************************** sfu chi^+ -> W^+ fu
          kp3=kw
          prtial(6)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par
          enddo
*********************************************** sfu^* chi^+ -> W^+ fubar
          kp3=kw
          prtial(7)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=3
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
**************************************************** sfu chi^+ -> H^+ fu
          kp3=khc
          prtial(8)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=5
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(8)=prtial(8)+par
          enddo
*********************************************** sfu^* chi^+ -> H^+ fubar
          kp3=khc
          prtial(9)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(9)=prtial(9)+par
          enddo
*********************************************** sfu^* chi^+ -> g fdbar
          if(cgluonin) then
            kp3=kgluon
            prtial(10)=0.d0
            do i=1,ncferd
              kp4=kcferd(i)
              icase=6
              call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
              prtial(10)=prtial(10)+par
            enddo
          else
            prtial(10)=0.0d0
          endif
**************************************************** sum partial results
          w=0.d0
          w=w+prtial(1)            ! sfu^* chi^+ -> Z fdbar
          w=w+prtial(2)            ! sfu^* chi^+ -> gamma fdbar
          w=w+prtial(3)            ! sfu^* chi^+ -> h fdbar
          w=w+prtial(4)            ! sfu^* chi^+ -> H fdbar
          w=w+prtial(5)            ! sfu^* chi^+ -> a fdbar
          w=w+prtial(6)            ! sfu chi^+ -> W^+ fu
          w=w+prtial(7)            ! sfu^* chi^+ -> W^+ fubar
          w=w+prtial(8)            ! sfu chi^+ -> H^+ fu
          w=w+prtial(9)            ! sfu^* chi^+ -> H^+ fubar
          w=w+prtial(10)           ! sfu^* chi^+ -> g fdbar
          dsasdwdcossfchi = w*0.5d0 ! 0.5d0 <- we should return weff for sf chi+ ann
                                    ! with sf and chi+ combined part. and anti-particle states
********
        elseif(ciaux.eq.2) then
****************************************************** sfd chi^+ -> Z fu
          kp3=kz
          prtial(1)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(1)=prtial(1)+par
          enddo
************************************************** sfd chi^+ -> gamma fu
          kp3=kgamma
          prtial(2)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=8
            call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
            prtial(2)=prtial(2)+par
          enddo
****************************************************** sfd chi^+ -> h fu
          kp3=kh2
          prtial(3)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(3)=prtial(3)+par
          enddo
****************************************************** sfd chi^+ -> H fu
          kp3=kh1
          prtial(4)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(4)=prtial(4)+par
          enddo  
****************************************************** sfd chi^+ -> a fu
          kp3=kh3
          prtial(5)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(5)=prtial(5)+par
          enddo
**************************************************** sfd chi^+ -> W^+ fd
          kp3=kw
          prtial(6)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
	    icase=5  
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par
          enddo
**************************************************** sfd chi^+ -> h^+ fd
          kp3=khc
          prtial(7)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=6
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
*********************************************** sfd^* chi^+ -> W^+ fdbar
          kp3=kw
          prtial(8)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(8)=prtial(8)+par
          enddo
*********************************************** sfd^* chi^+ -> h^+ fdbar
          kp3=khc
          prtial(9)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=5
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(9)=prtial(9)+par
          enddo
*********************************************** sfd chi^+ -> g fu
          if(cgluonin) then
            kp3=kgluon
            prtial(10)=0.d0
            do i=1,ncferd
              kp4=kcferd(i)
              icase=10
              call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
              prtial(10)=prtial(10)+par
            enddo 
          else
            prtial(10)=0.0d0
          endif
**************************************************** sum partial results
          w=0.d0
          w=w+prtial(1)            ! sfd chi^+ -> Z fu
          w=w+prtial(2)            ! sfd chi^+ -> gamma fu
          w=w+prtial(3)            ! sfd chi^+ -> h fu
          w=w+prtial(4)            ! sfd chi^+ -> H fu
          w=w+prtial(5)            ! sfd chi^+ -> a fu
          w=w+prtial(6)            ! sfd chi^+ -> W^+ fd
          w=w+prtial(7)            ! sfd chi^+ -> h^+ fd
          w=w+prtial(8)            ! sfd^* chi^+ -> W^+ fdbar
          w=w+prtial(9)            ! sfd^* chi^+ -> h^+ fdbar
          w=w+prtial(10)           ! sfd chi^+ -> g fu
          dsasdwdcossfchi = w*0.5d0 
               ! 0.5d0 <- we should return weff for sf chi+ ann
               ! with sf and chi+ combined part. and anti-particle states
********
        else
          write(*,*) 'DS: wrong ciaux value =',ciaux
          write(*,*) 'DS: in dsasdwdcossfchi. program stopped.'
          stop
        endif
      else 
        write(*,*) 'DS: in dsasdwdcossfchi wrong part. in init. state:'
        write(*,*) pname(kp2),' instead of neutralino or chargino.'
        write(*,*) 'DS: program stopped.'
        stop
      endif 


      if(aszeroprint) then
      wok=.true.
      do i=1,10
        if (prtial(i).lt.0.0d0) wok=.false.
      enddo
      if (.not.wok) then
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: dsasdwdcossfchi called with kp1=',kp1,'kp2=',kp2
        write(*,*) 'DS:  p=',p,' costh=',costhe,' w=',w
        do i=1,10
          write(*,*) 'DS: i= ',i,' prtial= ',prtial(i)
        enddo
      endif
      endif
  
      if (w.lt.0.0d0) w=abs(w)
      return
      end









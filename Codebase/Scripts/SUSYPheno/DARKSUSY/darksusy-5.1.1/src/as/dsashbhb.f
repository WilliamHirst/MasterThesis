c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.35, May 23, 2002, edsjo@physto.se)
c....Template file for dsashbhb begins here

**************************************************************
*** SUBROUTINE dsashbhb                                    ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + anti-sfermion(j)                         ***
*** -> higgs-boson + higgs-boson                           ***
***                                                        ***
*** The first mentioned particle (kp1) will be taken as    ***
*** a sfermion and the second particle (kp2) as an         ***
*** anti-sfermion -- not the opposite.                     ***
***                                                        ***
*** When kp1 and kp2 have different                        ***
*** weak isospin (T^3=+,-1/2), then kp1 must be an         ***
*** up-type-sfermion and kp2 a down-type-anti-sfermion.    ***
***                                                        ***
*** For the cases with one charged and one neutral higgs   ***
*** in the final state, the charged higgs must be          ***
*** mentioned first (i.e. kp3) and next the neutral        ***
*** higgs-boson (kp4) -- not the opposite.                 *** 
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-19                                         ***
*** Rewritten: 02-07-03 (because FORM input file now       ***
***  only gives the amplitude)                             ***
*** added flavour changing charged exchange for H^+H^-:    ***
*** added by Mia Schelke 2006-06-08  + Jan 2007            ***
*** added flavour changing charged exchange for:H^+H_1,2,3 ***
*** added by Mia Schelke Jan 2007                          ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************
***                                                        ***
*** NOTE: The FORM input file only gives the amplitude     ***
***       not the amplitude squared                        ***
*** THE FOLLOWING THEREFORE HAS TO BE CHANGED BY HAND      ***
*** after running of the PERL script                       ***
***                                                        ***
*** the form result should be denoted amplitude instead    ***
*** of dsashbhbas                                          ***
*** the amplitude squared calculation should be added      ***
*** -- the lines are written in the end of                 *** 
*** the template file                                      ***
**************************************************************

      subroutine dsashbhb(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      
      complex*16  amplitude
      complex*16  g4p
      complex*16  tsf,usf,sgba
      complex*16  sgbb,shb,point
      real*8      dsashbhbas
      real*8      dsabsq
      real*8      msferi2,msferj2,mhb12,mhb22
      real*8      morga(2)
      real*8      col
      integer     k,i
      integer     khb1,khb2
      integer     ksfert(6),ksferu(6),khbs(2),kgbs(2)
      integer     nsfert,nsferu,nhbs,ngbs,npoint
      logical nosneutrino
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      

***** check initial state is ok and set type  
c      if(abs(itype(1)-itype(2)).gt.1) then
c        write(*,*) 'DS: dsashbhb called with wrong particles'   
c        write(*,*) 'DS: in the initial state : ',pname(kp1),pname(kp2)  
c        stop  
c      endif    

***** set and check particles in the final state:
      khb1=kp3
      khb2=kp4
c      if(khb1.ne.kh1.and.khb1.ne.kh2.and.khb1.ne.kh3.and.
c     &   khb1.ne.khc) then
c        write(*,*) 'DS: dsashbhb called with kp3 = ',pname(khb1)  
c        write(*,*) 'DS: rather than a Higgs boson'
c        stop
c      endif 
c      if(khb2.ne.kh1.and.khb2.ne.kh2.and.khb2.ne.kh3.and.
c     &   khb2.ne.khc) then
c        write(*,*) 'DS: dsashbhb called with kp4 = ',pname(khb2)  
c        write(*,*) 'DS: rather than a Higgs boson'
c        stop
c      endif 

***** masses in final state:  
      mass3=mass(khb1)  
      mass4=mass(khb2)  
***** define the kinematic variables  
      call dsaskinset2

c....mass symbols used in form-code
      msferi2=mass(kp1)**2
      msferj2=mass(kp2)**2
      mhb12=mass(khb1)**2
      mhb22=mass(khb2)**2
      s=Svar
      t=Tvar
      u=Uvar

***** initially setting the arrays of exchange sfermions 
***** to be empty
      nsfert=0
      nsferu=0

*****
*****
******************* the first case **************************
***** sneutrino_i + slepton_i^+ -> H^+ + H^0_1/H^0_2
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> H^+ + H^0_1/H^0_2  with i,j = 1, 2, ... 6
      if(icase.eq.1) then
***** symmetry factor, for non-identical final state particles
        s34=1.d0
***** t- & u-channel, 6 squarks or 1,2 sleptons ***************
        nsfert=nsfertc
        do i=1,nsfert
          ksfert(i)=ksfertc(i)
        enddo
c        nsferu=nsfertc
        nsferu=nsferuc ! JE Correction 2009-04-03
        do i=1,nsferu
          ksferu(i)=ksferuc(i)
        enddo
***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kw
        kgbs(2)=0 
        ngbs=1
        morga(1)=1.d0/mass(kw)**2 
        morga(2)=0.d0
***** in the propagator, morga(k) is a factor that we have
***** defined for the term that is inversely proportional 
***** to the mass squared
***** morga(k) is mass(k)**-2 for a massive gauge boson in the 
***** exchange and 0 for the photon (or no exchange)
*****
***** in the s-channel with higgs exchange:
        khbs(1)=khc
        khbs(2)=0
        nhbs=1  
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,khb1,khb2)=g4p(f~_d,f~_u,H^+,H^0_1,2)
***** g4p interpretes g4p(1,2,3,4) as g4p(out,in,out,in)
***** (and it takes H^+ (not H^-) to give the direction of khc)
***** The vertex factor we want is constructed in g4p in  
***** the following way:
***** g4p(H^0_1,2,H^+,f~_u,f~_d) is defined in dsvertx3.f, 
***** the hermitian conjugated of this is 
***** g4p(H^+,H^0_1,2,f~_d,f~_u), i.e. (1<->2,3<->4)
***** which is identical to what we want (1<->3)+(2<->4)
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles 
*****
        goto 777
      endif
*****
*****
******************* the second case *************************
***** sneutrino_i + slepton_i^+ -> H^+ + H^0_3
***** up-type-squark_i + anti-down-type-squark_j 
*****     -> H^+ + H^0_3  with i,j = 1, 2, ... 6
      if(icase.eq.2) then
***** symmetry factor, for non-identical final state particles
        s34=1.d0  
***** t- & u-channel, 6 squarks or 1,2 sleptons ***************
        nsfert=nsfertc
        do i=1,nsfert
          ksfert(i)=ksfertc(i)
        enddo
        if(itype(1).eq.ivtype(ku)) then
c          nsferu=nsfertc
          nsferu=nsferuc ! JE Correction 2009-04-3
          do i=1,nsferu
            ksferu(i)=ksferuc(i)
          enddo
        else
***** (H^0_3 cannot couple to a pair of sneutrinos) 
           nsferu=0
        endif

***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kw
        kgbs(2)=0 
        ngbs=1
        morga(1)=1.d0/mass(kw)**2 
        morga(2)=0.d0
***** no s-channel with higgs exchange:
        khbs(1)=0
        khbs(2)=0
        nhbs=0
***** for the point interaction
        npoint=1
***** for g4p, a comment similar to the one for 
***** the previous case applies
*****
        goto 777
      endif  
*****
*****
******************* the third case **************************
***** slepton_i + slepton_i^* 
*****   -> H^0_1+H^0_1,H^0_2+H^0_2,H^0_1+H^0_2
***** up-type-squark_i + anti-up-type-squark_j 
*****   -> H^0_1+H^0_1,H^0_2+H^0_2,H^0_1+H^0_2
***** down-type-squark_i + anti-down-type-squark_j 
*****   -> H^0_1+H^0_1,H^0_2+H^0_2,H^0_1+H^0_2
*****     with i,j = 1, 2, ... 6
      if(icase.eq.3) then
***** symmetry factor, for identical or non-identical
***** final state particles
        if(khb1.eq.khb2) then
          s34=2.d0
        else
          s34=1.d0
        endif 
        nsfert=nsfertn
        do i=1,nsfert
          ksfert(i)=ksfertn(i)
        enddo
        nsferu=nsferun
        do i=1,nsferu
          ksferu(i)=ksferun(i)
        enddo
***** no s-channel with gauge-boson exchange:
        kgbs(1)=0
        kgbs(2)=0 
        ngbs=0
        morga(1)=0.d0 
        morga(2)=0.d0
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,khb1,khb2)=g4p(f~,f~,H^0_1,2,H^0_1,2)
***** These vertex factors are identical to the one set in dsvertx3.f:
***** for two identical higgs's we have g4p(H^0_i,H^0_i,f~,f~), which
***** is connected to what we have by (1<->3)+(2<->4)
***** for the case of two different higgs's, 
***** g4p(H^0_1,H^0_2,f~,f~) is set, and it must have identical 
***** vertex factor to g4p(H^0_2,H^0_1,f~,f~), since the two
***** higgs's are there own anti-particles
*****
        goto 777
      endif  
*****
*****
******************* the fourth case *************************
***** slepton_i + slepton_i^* 
*****   -> H^0_1+H^0_3, H^0_2+H^0_3
***** up-type-squark_i + anti-up-type-squark_j 
*****   -> H^0_1+H^0_3, H^0_2+H^0_3
***** down-type-squark_i + anti-down-type-squark_j 
*****   -> H^0_1+H^0_3, H^0_2+H^0_3
*****     with i,j = 1, 2, ... 6
      if(icase.eq.4) then
***** symmetry factor, for non-identical final state particles
        s34=1.d0  
        nosneutrino=.true.
        nsfert=nsfertn
        do i=1,nsfert
          ksfert(i)=ksfertn(i)
        enddo
        nsferu=nsferun
        do i=1,nsferu
          ksferu(i)=ksferun(i)
        enddo
        do i=1,3
          if(itype(1).eq.ivtype(2*i-1)) then
***** (H^0_3 cannot couple to a pair of sneutrinos) 
            nsfert=0
            nsferu=0
            nosneutrino=.false.
	  endif
        enddo
***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kz
        kgbs(2)=0 
        ngbs=1
        morga(1)=1.d0/mass(kz)**2
        morga(2)=0.d0
***** in the s-channel with higgs exchange:
        if(nosneutrino) then
          khbs(1)=kh3
          nhbs=1
        else
          khbs(1)=0
          nhbs=0
        endif
        khbs(2)=0
***** no point interaction
        npoint=0
*****
        goto 777
      endif 
*****
*****
******************* the fifth case **************************
***** slepton_i + slepton_i^* -> H^0_3+H^0_3
***** up-type-squark_i + anti-up-type-squark_j -> H^0_3+H^0_3
***** down-type-squark_i + anti-down-type-squark_j 
*****   -> H^0_3+H^0_3
*****     with i,j = 1, 2, ... 6
      if(icase.eq.5) then
***** symmetry factor, for identical final state particles
        s34=2.d0
        nsfert=nsfertn
        do i=1,nsfert
          ksfert(i)=ksfertn(i)
        enddo
        nsferu=nsferun
        do i=1,nsferu
          ksferu(i)=ksferun(i)
        enddo
        do i=1,3
          if(itype(1).eq.ivtype(2*i-1)) then
***** (H^0_3 cannot couple to a pair of sneutrinos) 
            nsfert=0
            nsferu=0
	  endif
        enddo
***** no s-channel with gauge-boson exchange:
        kgbs(1)=0
        kgbs(2)=0 
        ngbs=0
        morga(1)=0.d0 
        morga(2)=0.d0
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,khb1,khb2)=g4p(f~,f~,H^0_3,H^0_3)
***** The vertex factor we want is constructed in g4p in  
***** the following way:
***** g4p(H^0_3,H^0_3,f~,f~) is defined in dsvertx3.f,
***** and it is identical to the one we want, which is  
***** obtained by the particle interchange (1<->3)+(2<->4)
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles
***** 
        goto 777
      endif 
*****
*****
******************* the sixth case **************************
***** sneutrino_i + sneutrino_i^* -> H^+ + H^-
***** up-type-squark_i + anti-up-type-squark_j -> H^+ + H^-
*****     with i,j = 1, 2, ... 6
***** we take the first higgs boson 
***** (kp3, with four-momentum k1 in the form code) 
***** to be H^+ and the other one is then H^-
      if(icase.eq.6) then
***** symmetry factor, for non-identical final state particles
        s34=1.d0  
        nosneutrino=.true.
        nsfert=nsferuc
        do i=1,nsfert
          ksfert(i)=ksferuc(i)
        enddo
        do i=1,3
          if(itype(1).eq.ivtype(2*i-1)) then
            nosneutrino=.false.
	  endif
        enddo
***** no u-channel:  
        nsferu=0
***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kz
        morga(1)=1.d0/mass(kz)**2
        if(nosneutrino) then         
          kgbs(2)=kgamma 
          morga(2)=0.d0
          ngbs=2
        else
          kgbs(2)=0 
          morga(2)=0.d0
          ngbs=1
        endif
***** in the propagator, morga(k) is a factor that we have
***** defined for the term that is inversely proportional 
***** to the mass squared
***** morga(k) is mass(k)**-2 for a massive gauge boson in the 
***** exchange and 0 for the photon (or no exchange)
*****
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,khb1,khb2)=g4p(f~_u,f~_u,H^+_out,H^+_in)
***** =g4p(f~_u,f~_u,khc,khc)
***** g4p interpretes g4p(1,2,3,4) as g4p(out,in,out,in),
***** and it takes H^+ (not H^-) to give the direction of khc
***** so we have got it correct in this case
*****
        goto 777
      endif 
*****
*****
******************* the seventh case ************************
***** slepton_i^- + slepton_i^+ -> H^+ + H^-
***** down-type-squark_i + anti-down-type-squark_j 
*****   -> H^+ + H^-   with i,j = 1, 2, ... 6
***** we take the first higgs boson 
***** (kp3, with four-momentum k1 in the form code) 
***** to be H^+ and the other one is then H^-
      if(icase.eq.7) then
***** symmetry factor, for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** no t-channel:
        nsfert=0 
***** in the u-channel:
        nsferu=nsfertc
        do i=1,nsferu
          ksferu(i)=ksfertc(i)
        enddo
***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kz
        morga(1)=1.d0/mass(kz)**2
        kgbs(2)=kgamma 
        morga(2)=0.d0
        ngbs=2
***** in the propagator, morga(k) is a factor that we have
***** defined for the term that is inversely proportional 
***** to the mass squared
***** morga(k) is mass(k)**-2 for a massive gauge boson in the 
***** exchange and 0 for the photon (or no exchange)
*****
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
***** for the point interaction
        npoint=1
***** for the point interaction, a comment similar to the one 
***** for the sixth case applies
        goto 777
      endif 
*****
      write(*,*) 'DS: dsashbhb called with wrong icase : ',icase   
      write(*,*) 'DS: initial and final states : ',  
     &             pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  

 777  continue        


c....specify the coefficients used in the form-code
c....first for the t-channel(sfermion exchange)
      tsf=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
       do k=1,nsfert
        tsf=tsf+gl(khb1,ksfert(k),kp1)
     &    *(-dsasdepro(t,ksfert(k)))
     &    *gl(khb2,kp2,ksfert(k))
       enddo
      endif

c....then for the u-channel(sfermion exchange)
      usf=dcmplx(0.d0,0.d0)
      if(nsferu.gt.0) then
       do k=1,nsferu
        usf=usf+gl(khb2,ksferu(k),kp1)
     &    *(-dsasdepro(u,ksferu(k)))
     &    *gl(khb1,kp2,ksferu(k))
       enddo
      endif

c....then the s-channel with gauge-boson exchange
      sgba=dcmplx(0.d0,0.d0)
      sgbb=dcmplx(0.d0,0.d0)
      if(ngbs.gt.0) then
       do k=1,ngbs
        sgba=sgba+gl(kgbs(k),kp2,kp1)*dsasdepro(s,kgbs(k))
     &     *gl(kgbs(k),khb1,khb2)
        sgbb=sgbb+morga(k)*
     &     gl(kgbs(k),kp2,kp1)*(-dsasdepro(s,kgbs(k)))
     &     *gl(kgbs(k),khb1,khb2)
       enddo
      endif   

c....then the s-channel with higgs exchange
      shb=dcmplx(0.d0,0.d0)
      if(nhbs.gt.0) then
       do k=1,nhbs
        shb=shb+gl(khbs(k),kp2,kp1)
     &    *(-dsasdepro(s,khbs(k)))
     &    *mass(kw)*gl(khb2,khb1,khbs(k))
       enddo
      endif

c....then the point interaction
      point=dcmplx(0.d0,0.d0)
      if(npoint.gt.0) then
       point=g4p(kp2,kp1,khb1,khb2)
      endif     

***** now set the color factors to 1 for two sleptons in
***** the initial state and to 3 for two squarks in
***** the initial state  

      col=1.d0

      if(abs(itype(1)-ivtype(ku)).le.1) then
        col=3.d0     
      endif

***** After the template file follows the form expression of the
***** amplitude: called amplitude
***** NOTE: THE NAME HAS BEEN CHANGED  
***** from dsashbhbas to amplitude 
***** after running the PERL script
***** NOTE: this is NOT THE AMPLITUDE SQUARED 
***** JUST THE AMPLITUDE

***** Then we calculate the amplitude SQUARED 
***** and multiply it by the color factor
***** This is done by the following lines
***** which have to be PUT IN BY HAND
***** after running the PERL script

*      dsashbhbas=col*dsabsq(amplitude)
**    dsabsq is a function that calculates the 
**    absolute squared of a complex number
c....Template file for dsashbhb ends here



      amplitude =
     &  + msferi2*mhb12*sgbb
     &  - msferi2*mhb22*sgbb
     &  + msferi2*sgba
     &  - msferj2*mhb12*sgbb
     &  + msferj2*mhb22*sgbb
     &  + msferj2*sgba
     &  + mhb12*sgba
     &  + mhb22*sgba
     &  - s*sgba
     &  - 2.d0*t*sgba
     &  + tsf
     &  + usf
     &  + shb
     &  + point

      dsashbhbas=col*dsabsq(amplitude)
**    dsabsq is a function that calculates the 
**    absolute squared of a complex number

      if (dsashbhbas.lt.0.0d0) then
        if(aszeroprint) then
          write(*,*) ' '
          write(*,*) 'DS: negative cross section in dsashbhb:'
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
          write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
          write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
          write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
          write(*,*) 'DS: dsashbhbas = ',dsashbhbas
          write(*,*) 'DS: mass1=',mass1
          write(*,*) 'DS: mass2=',mass2
          write(*,*) 'DS: mass3=',mass3
          write(*,*) 'DS: mass4=',mass4
          write(*,*) 'DS: s=',s
          write(*,*) 'DS: t=',t
          write(*,*) 'DS: u=',u
          write(*,*) 'DS: p12=',p12
          write(*,*) 'DS: costheta=',costheta
        endif
        dsashbhbas=0.0d0
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsashbhbas)*k34/(8.0d0*pi*gg1*gg2*s34*dsqrt(s))
      return
      end



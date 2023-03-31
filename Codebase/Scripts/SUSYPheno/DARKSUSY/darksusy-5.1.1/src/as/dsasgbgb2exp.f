c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.35, May 23, 2002, edsjo@physto.se)
c....Template file for dsasgbgb2exp begins here

**************************************************************
*** SUBROUTINE dsasgbgb2exp                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + anti-sfermion(j)                         ***
*** -> gluon+gluon, photon+photon, photon+gluon            ***
***                                                        ***
*** The first mentioned particle (kp1) will be taken as    ***
*** a sfermion and the second particle (kp2) as an         ***
*** anti-sfermion -- not the opposite.                     ***
***                                                        ***
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 02-05-21                                         ***
*** Rewritten: 02-07-03 (to have exactly two massless gb)  ***  
*** explicite pol. vectors introduced: 02-07-08            ***
*** sum over pol. moved from fortran to form: 02-07-09     *** 
*** two colour factors(c.f.) made complex: 02-07-10        ***
***(+these c.f. changed for g+g as ggg vertex code changed)*** 
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      subroutine dsasgbgb2exp(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      real*8      dsasgbgb2expas
      real*8      mi2,mj2
      real*8      f(2)
      complex*16  g4p
      complex*16  tsf,usf,sgba,point
      real*8      tsfctsf,usfcusf,sgbacsgba,pointcpoint,
     &            tsfcusf2r,tsfcsgba2r,tsfcpoint2r,
     &            usfcsgba2r,usfcpoint2r,
     &            sgbacpoint2r
      real*8      ctt,cuu,csgbsgb,cpp,
     &            ctu,ctp,cup,csgbp
      complex*16  ctsgb,cusgb
      real*8      ek11p1,ek11p2,ek21p1,ek21p2
      integer     kgb1,kgb2
      integer     k,i
      integer     ksfert(6),ksferu(6),kgbs(2)
      integer     nsfert,nsferu,ngbs,npoint
      real*8 pi
      parameter (pi=3.141592653589793238d0)


c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      

***** set particles in the final state:
      kgb1=kp3
      kgb2=kp4

***** masses in final state:  
      mass3=mass(kgb1)  
      mass4=mass(kgb2)  
***** define the kinematic variables  
      call dsaskinset2


c....mass symbols used in form-code
      mi2=mass(kp1)**2
      mj2=mass(kp2)**2
      s=Svar  
      t=Tvar  
      u=Uvar  


***** now set the color factors to 1.d0 for two sleptons in
***** the initial state and to 3.d0 for two squarks in
***** the initial state
***** if one of the gauge bosons in the final state is a gluon, then the
***** color factor is changed later on

      ctt=1.d0
      cuu=1.d0
      csgbsgb=1.d0
      cpp=1.d0
      ctu=1.d0
      ctsgb=dcmplx(1.d0,0.d0)
      ctp=1.d0
      cusgb=dcmplx(1.d0,0.d0)
      cup=1.d0
      csgbp=1.d0

      if(abs(itype(1)-ivtype(ku)).le.1) then
        ctt=3.d0
        cuu=3.d0
        csgbsgb=3.d0
        cpp=3.d0
        ctu=3.d0
        ctsgb=dcmplx(3.d0,0.d0)
        ctp=3.d0
        cusgb=dcmplx(3.d0,0.d0)
        cup=3.d0
        csgbp=3.d0
      endif  

***** now look at the different initial/final state cases  

***** sfermion(i) + anti-sfermion(j)
***** -> photon + photon
***** (the third case in the original dsasgbgb code)
***** no sneutrinos are allowed in the initial state
      if(icase.eq.1) then
*****
***** set the symmetry factor for identical final state particles
        s34=2.d0
***** now set particles in the intermediate states
***** in the t-channel:  
        nsfert=nsfertn
        do i=1,nsfert
        ksfert(i)=ksfertn(i)
        enddo
***** in the u-channel:  
        nsferu=nsferun
        do i=1,nsferu
        ksferu(i)=ksferun(i)
        enddo
***** no s-channel with gauge-boson exchange:
        kgbs(1)=0
        f(1)=0.d0
        kgbs(2)=0
        f(2)=0.d0
        ngbs=0
***** for the point interaction  
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~,f~,gamma,gamma)
***** The vertex factor we want is constructed in g4p in
***** the following way:
***** g4p(gamma,gamma,f~,f~) is defined in dsvertx3.f,
***** and it is identical to the one we want, which is
***** obtained by the particle interchange (1<->3)+(2<->4)
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles
*****
        goto 567
      endif  

***** squark(i) + anti-squark(j) 
***** -> gluon + gluon
***** (the eighth case in the original dsasgbgb code)
      if(icase.eq.2) then
***** set the color factors
        ctt=16.d0/3.d0
        cuu=16.d0/3.d0
        csgbsgb=12.d0
        cpp=28.d0/3.d0
        ctu=-2.d0/3.d0
        ctsgb=dcmplx(0.d0,-6.d0)
*ctsgb changed 02-07-10: from 6 to (-6i)
*it is the color factor of M_t*(M_s)^\dagger
*color factor of M_s*(M_t)^\dagger is (6i)
*the imaginary unit from the ggg-coupling was included in the color
*factor before, which is why it was real, but now it is included in 
*the vertex factor code
        ctp=14.d0/3.d0
        cusgb=dcmplx(0.d0,6.d0)
*cusgb changed 02-07-10: from (-6) to (6i)
*it is the color factor of M_u*(M_s)^\dagger
*color factor of M_s*(M_u)^\dagger is (-6i)
*the imaginary unit from the ggg-coupling was included in the color
*factor before, which is why it was real, but now it is included in 
*the vertex factor code
        cup=14.d0/3.d0
        csgbp=0.d0
*****
***** set the symmetry factor for identical final state particles
        s34=2.d0  
***** now set particles in the intermediate states 
***** in the t-channel:  
        nsfert=nsfertn
        do i=1,nsfert
        ksfert(i)=ksfertn(i)
        enddo
***** in the u-channel:  
        nsferu=nsferun
        do i=1,nsferu
        ksferu(i)=ksferun(i)
        enddo
***** in the s-channel with gauge-boson exchange:
        kgbs(1)=kgluon
        f(1)=1.d0
        kgbs(2)=0
        f(2)=0.d0
        ngbs=1
***** for the point interaction
        npoint=1
*****
        goto 567
      endif 

*****
***** squark(i) + anti-squark(j)
***** -> gluon + photon
***** (the nineth case in the original dsasgbgb code)
      if(icase.eq.3) then
***** set the color factors
        ctt=4.d0
        cuu=4.d0
        csgbsgb=0.d0
        cpp=4.d0
        ctu=4.d0
        ctsgb=dcmplx(0.d0,0.d0)
        ctp=4.d0
        cusgb=dcmplx(0.d0,0.d0)
        cup=4.d0
        csgbp=0.d0
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0
***** now set particles in the intermediate states
***** in the t-channel:  
        nsfert=nsfertn
        do i=1,nsfert
        ksfert(i)=ksfertn(i)
        enddo
***** in the u-channel:  
        nsferu=nsferun
        do i=1,nsferu
        ksferu(i)=ksferun(i)
        enddo
***** no s-channel with gauge-boson exchange:
        kgbs(1)=0
        f(1)=0.d0
        kgbs(2)=0
        f(2)=0.d0
        ngbs=0
*****
***** for the point interaction
        npoint=1
        goto 567
      endif  

      write(*,*) 'DS: dsasgbgb2exp called with wrong icase : ',icase   
      write(*,*) 'DS: initial or final state'  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  

 567  continue     


c....specify the coefficients used in the form-code
c....first for the t-channel(sfermion exchange)
      tsf=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
       do k=1,nsfert
        tsf=tsf+gl(kgb1,ksfert(k),kp1)
     &    *(-dsasdepro(t,ksfert(k)))
     &    *gl(kgb2,kp2,ksfert(k))
       enddo
      endif

c....then the u-channel(sfermion exchange)
      usf=dcmplx(0.d0,0.d0)
      if(nsferu.gt.0) then
       do k=1,nsferu
        usf=usf+gl(kgb2,ksferu(k),kp1)
     &    *(-dsasdepro(u,ksferu(k)))
     &    *gl(kgb1,kp2,ksferu(k))
       enddo  
      endif

c....then the s-channel with gauge-boson exchange
      sgba=dcmplx(0.d0,0.d0)
      if(ngbs.gt.0) then
       do k=1,ngbs
        sgba=sgba+gl(kgbs(k),kp2,kp1)
     &     *dsasdepro(s,kgbs(k))
     &     *f(k)*gl(kgb1,kgb2,kgbs(k))
       enddo
      endif

c....then the point interaction
      point=dcmplx(0.d0,0.d0)
      if(npoint.gt.0) then
       point=g4p(kp2,kp1,kgb1,kgb2)
      endif    

c....then all the real combinations of channels
c....used in form code
c....the color factors ctt,utt etc. are included here
c....note:the two color factors ctsgb and cusgb can be complex 
      tsfctsf=ctt*dble(tsf*conjg(tsf))
      usfcusf=cuu*dble(usf*conjg(usf))
      sgbacsgba=csgbsgb*dble(sgba*conjg(sgba))
      pointcpoint=cpp*dble(point*conjg(point))
      tsfcusf2r=ctu*2*dble(tsf*conjg(usf))
      tsfcsgba2r=2*dble(ctsgb*tsf*conjg(sgba))
      tsfcpoint2r=ctp*2*dble(tsf*conjg(point))
      usfcsgba2r=2*dble(cusgb*usf*conjg(sgba))
      usfcpoint2r=cup*2*dble(usf*conjg(point))
      sgbacpoint2r=csgbp*2*dble(sgba*conjg(point))

*Now set the 4-scalar product of the the 1st
*pol.vector of the first gb with p1 and p2
*for the 2nd pol. vector these numbers vanish
*and this was already specified in form
      ek11p1=p12*dsqrt(1.d0-costheta**2)
      ek11p2=-p12*dsqrt(1.d0-costheta**2)
*similarly for the pol. vectors of the other gb
      ek21p1=p12*dsqrt(1.d0-costheta**2)
      ek21p2=-p12*dsqrt(1.d0-costheta**2)
*p12=p in dsasdwdcossfsf.f
*costheta=costhe in dsasdwdcossfsf.f


***** After the template file follows the form expression of the 
***** amplitude squared: dsasgbgb2expas

c....Template file for dsasgbgb2exp ends here




      dsasgbgb2expas =
     &  + 4.d0*mi2*mj2*sgbacsgba
     &  - 4.d0*mi2*s*sgbacsgba
     &  - 8.d0*mi2*t*sgbacsgba
     &  + 2.d0*mi2*ek11p1*ek21p1*tsfcsgba2r
     &  + 2.d0*mi2*ek11p1*ek21p1*usfcsgba2r
     &  + 4.d0*mi2*ek11p1*ek21p2*sgbacsgba
     &  - 2.d0*mi2*ek11p1*ek21p2*tsfcsgba2r
     &  - 4.d0*mi2*ek11p2*ek21p1*sgbacsgba
     &  - 2.d0*mi2*ek11p2*ek21p1*usfcsgba2r
     &  - 2.d0*mi2*sgbacpoint2r
     &  + 2.d0*mi2**2*sgbacsgba
     &  - 4.d0*mj2*s*sgbacsgba
     &  - 8.d0*mj2*t*sgbacsgba
     &  + 2.d0*mj2*ek11p1*ek21p1*tsfcsgba2r
     &  + 2.d0*mj2*ek11p1*ek21p1*usfcsgba2r
     &  + 4.d0*mj2*ek11p1*ek21p2*sgbacsgba
     &  - 2.d0*mj2*ek11p1*ek21p2*tsfcsgba2r
     &  - 4.d0*mj2*ek11p2*ek21p1*sgbacsgba
     &  - 2.d0*mj2*ek11p2*ek21p1*usfcsgba2r
     &  - 2.d0*mj2*sgbacpoint2r
     &  + 2.d0*mj2**2*sgbacsgba
     &  + 8.d0*s*t*sgbacsgba
     &  - 2.d0*s*ek11p1*ek21p1*tsfcsgba2r
     &  - 2.d0*s*ek11p1*ek21p1*usfcsgba2r
     &  - 4.d0*s*ek11p1*ek21p2*sgbacsgba
     &  + 2.d0*s*ek11p1*ek21p2*tsfcsgba2r
     &  + 4.d0*s*ek11p2*ek21p1*sgbacsgba
     &  + 2.d0*s*ek11p2*ek21p1*usfcsgba2r
     &  + 2.d0*s*sgbacpoint2r
     &  + 2.d0*s**2*sgbacsgba
     &  - 4.d0*t*ek11p1*ek21p1*tsfcsgba2r
     &  - 4.d0*t*ek11p1*ek21p1*usfcsgba2r
     &  - 8.d0*t*ek11p1*ek21p2*sgbacsgba
     &  + 4.d0*t*ek11p1*ek21p2*tsfcsgba2r
     &  + 8.d0*t*ek11p2*ek21p1*sgbacsgba
     &  + 4.d0*t*ek11p2*ek21p1*usfcsgba2r
     &  + 4.d0*t*sgbacpoint2r
     &  + 8.d0*t**2*sgbacsgba
     &  - 8.d0*ek11p1*ek11p2*ek21p1*ek21p2*sgbacsgba
     &  + 4.d0*ek11p1*ek11p2*ek21p1*ek21p2*tsfcusf2r
     &  + 4.d0*ek11p1*ek11p2*ek21p1*ek21p2*tsfcsgba2r
     &  - 4.d0*ek11p1*ek11p2*ek21p1*ek21p2*usfcsgba2r
     &  - 8.d0*ek11p1*ek11p2*ek21p1**2*usfcusf
     &  - 4.d0*ek11p1*ek11p2*ek21p1**2*tsfcusf2r
     &  - 4.d0*ek11p1*ek11p2*ek21p1**2*tsfcsgba2r
     &  - 4.d0*ek11p1*ek11p2*ek21p1**2*usfcsgba2r
     &  - 2.d0*ek11p1*ek21p1*tsfcpoint2r
     &  - 2.d0*ek11p1*ek21p1*usfcpoint2r
     &  + 2.d0*ek11p1*ek21p2*tsfcpoint2r
     &  - 2.d0*ek11p1*ek21p2*sgbacpoint2r
     &  - 8.d0*ek11p1**2*ek21p1*ek21p2*tsfctsf
     &  - 4.d0*ek11p1**2*ek21p1*ek21p2*tsfcusf2r
     &  + 4.d0*ek11p1**2*ek21p1*ek21p2*tsfcsgba2r
     &  + 4.d0*ek11p1**2*ek21p1*ek21p2*usfcsgba2r
     &  + 4.d0*ek11p1**2*ek21p1**2*tsfctsf
     &  + 4.d0*ek11p1**2*ek21p1**2*usfcusf
     &  + 4.d0*ek11p1**2*ek21p1**2*tsfcusf2r
     &  + 4.d0*ek11p1**2*ek21p2**2*tsfctsf
     &  + 4.d0*ek11p1**2*ek21p2**2*sgbacsgba
     &  - 4.d0*ek11p1**2*ek21p2**2*tsfcsgba2r
     &  + 2.d0*ek11p2*ek21p1*usfcpoint2r
     &  + 2.d0*ek11p2*ek21p1*sgbacpoint2r
     &  + 4.d0*ek11p2**2*ek21p1**2*usfcusf
     &  + 4.d0*ek11p2**2*ek21p1**2*sgbacsgba
     &  + 4.d0*ek11p2**2*ek21p1**2*usfcsgba2r
     &  + 2.d0*pointcpoint


      if (dsasgbgb2expas.lt.0.0d0) then  
        if(aszeroprint) then
          write(*,*) ' '
          write(*,*) 'DS: negative cross section in dsasgbgb2exp:'
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
          write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
          write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
          write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
          write(*,*) 'DS: dsasgbgb2expas = ',dsasgbgb2expas
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
        dsasgbgb2expas=0.0d0
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsasgbgb2expas)*k34/(8.0d0*pi*gg1*gg2*s34*dsqrt(s))
      return
      end

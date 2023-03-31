c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.35, May 23, 2002, edsjo@physto.se)
c....Template file for dsasgbgb1exp begins here

**************************************************************
*** SUBROUTINE dsasgbgb1exp                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + anti-sfermion(j)                         ***
*** -> massive gauge-boson + massless gauge-boson          ***
***                                                        ***
***                                                        ***
*** The first mentioned particle (kp1) will be taken as    ***
*** a sfermion and the second particle (kp2) as an         ***
*** anti-sfermion -- not the opposite.                     ***
***                                                        ***
*** When kp1 and kp2 have different                        ***
*** weak isospin (T^3=+,-1/2), then kp1 must be an         ***
*** up-type-sfermion and kp2 a down-type-anti-sfermion.    ***
***                                                        ***
*** NOTE: for the gauge bosons, the MASSIVE must be        ***
*** mentioned first (kp3) and then the MASSLESS one (kp4)  ***
*** -- not the opposite.                                   ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-23  rewritten:02-03-12                     ***
*** QCD included: 02-03-20                                 ***
*** rewritten: 02-07-04 (to have exactly one massless gb)  ***
*** explicite pol. vectors introduced: 02-07-05            ***
*** sum over massless pol. moved from                      ***
*** fortran to form: 02-07-09                              ***
***                                                        ***
*** addition: Jan 2007 by Mia Schelke                      ***
*** squark mixing allowed for W+ gamma,gluon               ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      subroutine dsasgbgb1exp(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      real*8      dsasgbgb1expas
      real*8      mi2,mj2,gb
      real*8      d,smorga(2),f(2)
*symbol d was called mog1 in previous code
*but renamed to give shorter form output lines 
      complex*16  g4p
      complex*16  tsf,usf,sgba,sgbb,point
      real*8      tsfctsf,usfcusf,sgbacsgba,
     &            sgbbcsgbb,pointcpoint,
     &            tsfcusf2r,tsfcsgba2r,
     &            tsfcsgbb2r,tsfcpoint2r,
     &            usfcsgba2r,usfcsgbb2r,usfcpoint2r,
     &            sgbacsgbb2r,sgbacpoint2r,
     &            sgbbcpoint2r
      real*8      ctt,cuu,csgbsgb,cpp,
     &            ctu,ctsgb,ctp,
     &            cusgb,cup,csgbp
      real*8      eps1p1,eps1p2
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

***** check initial state is ok and set type  
c      if(abs(itype(2)-itype(1)).gt.1) then
c        write(*,*) 'dsasgbgb1exp called with wrong particles'   
c        write(*,*) 'in the initial state :'  
c        write(*,*) pname(kp1),pname(kp2)  
c        stop  
c      endif    


***** set and check particles in the final state:
      kgb1=kp3
      kgb2=kp4
c      if(kgb1.ne.kz.and.kgb1.ne.kw) then
c        write(*,*) 'dsasgbgb1exp called with kp3 ='   
c        write(*,*) pname(kgb1)  
c        write(*,*) 'instead of a MASSIVE gauge boson'
c        stop
c      endif 
c      if(kgb2.ne.kgamma.and.kgb2.ne.kgluon) then
c         write(*,*) 'dsasgbgb1exp called with kp4 ='   
c         write(*,*) pname(kgb2)  
c         write(*,*) 'instead of a MASSLESS gauge boson'
c        stop
c      endif 

***** masses in final state:  
      mass3=mass(kgb1)  
      mass4=mass(kgb2)  
***** define the kinematic variables  
      call dsaskinset2   


c....mass symbols used in form-code
      mi2=mass(kp1)**2
      mj2=mass(kp2)**2
      gb=mass(kgb1)**2
      s=Svar  
      t=Tvar  
      u=Uvar  

***** now set the color factors to 1.d0 for two sleptons in
***** the initial state and to 3.d0 for two squarks in 
***** the initial state
***** if one of the gauge bosons in the final state is a gluon, 
***** then the color factor is changed later on 

      ctt=1.d0
      cuu=1.d0
      csgbsgb=1.d0
      cpp=1.d0
      ctu=1.d0
      ctsgb=1.d0
      ctp=1.d0
      cusgb=1.d0
      cup=1.d0
      csgbp=1.d0

      if(abs(itype(1)-ivtype(ku)).le.1) then
        ctt=3.d0
        cuu=3.d0
        csgbsgb=3.d0
        cpp=3.d0
        ctu=3.d0
        ctsgb=3.d0
        ctp=3.d0
        cusgb=3.d0
        cup=3.d0
        csgbp=3.d0
      endif

*****
*****
***** the first case    
***** up-type-sfermion(i) + anti-down-type-sfermion(j) 
***** -> W^+ + photon
      if(icase.eq.1) then
        d=1.d0/mass(kgb1)**2
***** in the completeness relation for the polarization vectors 
***** d is a factor that we have defined for the 
***** terms that are inversely proportional to the mass squared
***** d is mass(kgb1)**-2 for the massive gauge boson in the 
***** final state 
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** in the t-channel (2 sfermions):  
        nsfert=nsfertc
        do i=1,nsfert
          ksfert(i)=ksfertc(i)
        enddo
***** in the u-channel (2 sfermions or none):  
        if(itype(1).eq.ivtype(ku)) then
          nsferu=nsferuc
          do i=1,nsferu
            ksferu(i)=ksferuc(i)
          enddo
        else
          nsferu=0
        endif
***** in the s-channel with gauge-boson exchange:
***** We will write the trilinear gauge boson coupling on the form:
***** f(k)*gl(kgb1,kgb2,kgbs(k))*mom(kgb1,kgb2,kgbs(k)).
***** Here mom(kgb1,kgb2,kgbs(k)) is 'the momentum part', 
***** that can be found in the form code, and it depends on
***** the particle momenta in an even permutation of the order 
***** in which they are mentioned.
***** The 'structure constant',f(k), will be +,-1 depending on
***** the 'darkSUSY' interpretation of gl(kgb1,kgb2,kgbs(k).
***** In the present case we have:
***** gl(kgb1,kgb2,kgbs(k)=gl(W^+_out,photon,W^+_in)
***** =gl(W^-_in,photon,W^-_out)
***** The 'problem' is that when we write gl(kw,kgamma,kw), then 
***** 'darkSUSY' will instead interpret this as 
***** gl(W^-_out,photon,W^-_in), i.e. an odd permutation of the 
***** gauge bosons compared to the correct interpretation.
***** The difference between even and odd permutations is 
***** a minus sign (that would cancel the minus sign you 
***** would get if you also changed from even to odd permutation
***** in the momentum part).
***** Therefore, we must set the 'structure constant', f(1), 
***** equal to (-1), to get the correct result.
        kgbs(1)=kw
        smorga(1)=1.d0/mass(kw)**2
        f(1)=-1.d0
        kgbs(2)=0
        smorga(2)=0.d0
        f(2)=0.d0
        ngbs=1
***** in the propagator, smorga(k) is a factor that we have
***** defined for the term that is inversely proportional to 
***** the mass squared 
***** smorga(k) is mass(k)**-2 for a massive gauge boson in the 
***** exchange and 0 for the photon (or no exchange)
*****
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~_d,f~_u,W^+_out,gamma)
***** g4p interpretates g4p(1,2,3,4) as g4p(out,in,out,in),
***** and it takes W^+ (i.e. not W^- as above!) to give the direction of kw,
***** so we have got it correct in this case
***** but this is only because we did take the order of the 
***** gauge bosons to be (kgb1,kgb2)
***** (If we had written g4p(f~_d,f~_u,gamma,W^-_in)
***** =g4p(f~_d,f~_u,kgamma,kw), then g4p would have interpreted
***** it as g4p(f~_d_out,f~_u_in,gamma_out,W^+_in), which violates 
***** charge conservation)
***** The vertex factor we want, is constructed in g4p in
***** the following way:
***** g4p(kgamma,kw,f~_u,f~_d) is defined in dsvertx3.f,
***** the hermitian conjugated of this is
***** g4p(kw,kgamma,f~_d,f~_u), i.e. (1<->2,3<->4)
***** which is identical to what we want (1<->3)+(2<->4),  
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles   
*****
        goto 555
      endif   

***** the second case 
***** (= the fourth case in the original dsasgbgb file)
***** sfermion(i) + anti-sfermion(j) 
***** -> Z + photon 
***** no sneutrinos are allowed in the initial state
      if(icase.eq.2) then
        d=1.d0/mass(kgb1)**2
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
        smorga(1)=0.d0
        f(1)=0.d0
        kgbs(2)=0
        smorga(2)=0.d0
        f(2)=0.d0
        ngbs=0
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p either like 
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~,f~,gamma,Z)
***** or g4p(f~,f~,Z,gamma)
***** The vertex factors we want are constructed in g4p in
***** the following way:
***** g4p(gamma,Z,f~,f~) and g4p(Z,gamma,f~,f~)
***** are defined in dsvertx3.f (and are identical),
***** and they are identical to those we have, that are 
***** obtained by the particle interchange (1<->3)+(2<->4)
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles  
*****
        goto 555
      endif 

*****
***** the third case
***** (= the tenth case in the original dsasgbgb file)
***** squark(i) + anti-squark(j) 
***** -> Z + gluon 
      if(icase.eq.3) then
        d=1.d0/mass(kgb1)**2
***** set the color factors
        ctt=4.d0
        cuu=4.d0
        csgbsgb=0.d0
        cpp=4.d0
        ctu=4.d0
        ctsgb=0.d0
        ctp=4.d0
        cusgb=0.d0
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
        smorga(1)=0.d0
        f(1)=0.d0
        kgbs(2)=0 
        smorga(2)=0.d0
        f(2)=0.d0
        ngbs=0
*****
***** for the point interaction
        npoint=1
        goto 555
      endif
***** 
***** the fourth case
***** (= the eleventh case in the original dsasgbgb file)     
***** up-type-squark(i) + anti-down-type-squark(j) 
***** -> W^+ + gluon
      if(icase.eq.4) then
        d=1.d0/mass(kgb1)**2
***** in the completeness relation for the polarization vectors 
***** d is a factor that we have defined for the 
***** terms that are inversely proportional to the mass squared
***** d is mass(kgb1)**-2 for the massive gauge boson in the 
***** final state 
*****
***** set the color factors
        ctt=4.d0
        cuu=4.d0
        csgbsgb=0.d0
        cpp=4.d0
        ctu=4.d0
        ctsgb=0.d0
        ctp=4.d0
        cusgb=0.d0
        cup=4.d0
        csgbp=0.d0     
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** in the t-channel (2 sfermions):
        nsfert=nsfertc
        do i=1,nsfert
          ksfert(i)=ksfertc(i)
        enddo
***** in the u-channel (2 sfermions or none):  
        nsferu=nsferuc
        do i=1,nsferu
          ksferu(i)=ksferuc(i)
        enddo
***** no s-channel with gauge-boson exchange:
        kgbs(1)=0
        smorga(1)=0.d0
        f(1)=0.d0
        kgbs(2)=0
        smorga(2)=0.d0
        f(2)=0.d0
        ngbs=0
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~_d,f~_u,W^+_out,gluon)
***** g4p interpretates g4p(1,2,3,4) as g4p(out,in,out,in),
***** and it takes W^+ (i.e. not W^-) to give the direction of kw,
***** so we have got it correct in this case
        goto 555
      endif   

      write(*,*) 'DS: dsasgbgb1exp called with wrong icase : ',icase   
      write(*,*) 'DS: initial and  final states : '  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  


 555  continue


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
      sgbb=dcmplx(0.d0,0.d0)
      if(ngbs.gt.0) then
       do k=1,ngbs
        sgba=sgba+gl(kgbs(k),kp2,kp1)
     &     *dsasdepro(s,kgbs(k))
     &     *f(k)*gl(kgb1,kgb2,kgbs(k))
        sgbb=sgbb+smorga(k)*
     &     gl(kgbs(k),kp2,kp1)
     &     *(-dsasdepro(s,kgbs(k)))
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

      tsfctsf=ctt*dble(tsf*conjg(tsf))
      usfcusf=cuu*dble(usf*conjg(usf))
      sgbacsgba=csgbsgb*dble(sgba*conjg(sgba))
      sgbbcsgbb=csgbsgb*dble(sgbb*conjg(sgbb))
      pointcpoint=cpp*dble(point*conjg(point))
      tsfcusf2r=ctu*2*dble(tsf*conjg(usf))
      tsfcsgba2r=ctsgb*2*dble(tsf*conjg(sgba))
      tsfcsgbb2r=ctsgb*2*dble(tsf*conjg(sgbb))
      tsfcpoint2r=ctp*2*dble(tsf*conjg(point))
      usfcsgba2r=cusgb*2*dble(usf*conjg(sgba))
      usfcsgbb2r=cusgb*2*dble(usf*conjg(sgbb))
      usfcpoint2r=cup*2*dble(usf*conjg(point))
      sgbacsgbb2r=csgbsgb*2*dble(sgba*conjg(sgbb))
      sgbacpoint2r=csgbp*2*dble(sgba*conjg(point))
      sgbbcpoint2r=csgbp*2*dble(sgbb*conjg(point))      


*Now set the 4-scalar product of the the 1st
*pol.vector of the massless gb with p1 and p2
*for the 2nd pol. vector these numbers vanish
*and this was already specified in form
      eps1p1=p12*dsqrt(1.d0-costheta**2)
      eps1p2=-p12*dsqrt(1.d0-costheta**2)
*p12=p in dsasdwdcossfsf.f
*costheta=costhe in dsasdwdcossfsf.f 	
    

***** After the template file follows the form expression of the 
***** amplitude squared: dsasgbgb1expas

c....Template file for dsasgbgb1exp ends here



      dsasgbgb1expas =
     &  - 4.d0*mi2*mj2*gb*eps1p1*eps1p2*sgbbcsgbb
     &  - 2.d0*mi2*mj2*gb*eps1p1**2*sgbbcsgbb
     &  - 2.d0*mi2*mj2*gb*eps1p2**2*sgbbcsgbb
     &  - 4.d0*mi2*mj2*gb**2*sgbbcsgbb
     &  - 4.d0*mi2*mj2*d*eps1p1*eps1p2*sgbacsgba
     &  + 2.d0*mi2*mj2*d*eps1p1*eps1p2*tsfcsgba2r
     &  - 2.d0*mi2*mj2*d*eps1p1*eps1p2*usfcsgba2r
     &  - 2.d0*mi2*mj2*d*eps1p1**2*sgbacsgba
     &  - 2.d0*mi2*mj2*d*eps1p1**2*usfcsgba2r
     &  - 2.d0*mi2*mj2*d*eps1p2**2*sgbacsgba
     &  + 2.d0*mi2*mj2*d*eps1p2**2*tsfcsgba2r
     &  - 2.d0*mi2*mj2*eps1p1*eps1p2*tsfcsgbb2r
     &  + 2.d0*mi2*mj2*eps1p1*eps1p2*usfcsgbb2r
     &  + 4.d0*mi2*mj2*eps1p1*eps1p2*sgbacsgbb2r
     &  + 2.d0*mi2*mj2*eps1p1**2*usfcsgbb2r
     &  + 2.d0*mi2*mj2*eps1p1**2*sgbacsgbb2r
     &  - 2.d0*mi2*mj2*eps1p2**2*tsfcsgbb2r
     &  + 2.d0*mi2*mj2*eps1p2**2*sgbacsgbb2r
     &  + 4.d0*mi2*mj2*sgbacsgba
     &  - 2.d0*mi2*gb*s*sgbacsgbb2r
     &  - 4.d0*mi2*gb*t*sgbacsgbb2r
     &  - 4.d0*mi2*gb*d*eps1p1*eps1p2*tsfcusf2r
     &  + 2.d0*mi2*gb*d*eps1p1*eps1p2*usfcsgba2r
     &  + 8.d0*mi2*gb*d*eps1p1**2*usfcusf
     &  + 2.d0*mi2*gb*d*eps1p1**2*usfcsgba2r
     &  - 2.d0*mi2*gb*eps1p1*eps1p2*tsfcsgbb2r
     &  - 4.d0*mi2*gb*eps1p1*eps1p2*usfcsgbb2r
     &  + mi2*gb*eps1p1**2*sgbacsgbb2r
     &  + 2.d0*mi2*gb*eps1p2**2*tsfcsgbb2r
     &  - mi2*gb*eps1p2**2*sgbacsgbb2r
     &  + 4.d0*mi2*gb*sgbacsgba
     &  - 2.d0*mi2*gb*sgbbcpoint2r
     &  + 2.d0*mi2*gb**2*sgbacsgbb2r
     &  + 4.d0*mi2*s*d*eps1p1*eps1p2*tsfcusf2r
     &  + 2.d0*mi2*s*d*eps1p1*eps1p2*tsfcsgba2r
     &  - 8.d0*mi2*s*d*eps1p1**2*usfcusf
     &  - 2.d0*mi2*s*d*eps1p1**2*sgbacsgba
     &  - 4.d0*mi2*s*d*eps1p1**2*usfcsgba2r
     &  + 2.d0*mi2*s*d*eps1p2**2*sgbacsgba
     &  - 2.d0*mi2*s*d*eps1p2**2*tsfcsgba2r
     &  + 2.d0*mi2*s*eps1p1*eps1p2*usfcsgbb2r
     &  + 2.d0*mi2*s*eps1p1**2*usfcsgbb2r
     &  + mi2*s*eps1p1**2*sgbacsgbb2r
     &  - mi2*s*eps1p2**2*sgbacsgbb2r
     &  - 4.d0*mi2*s*sgbacsgba
     &  + 8.d0*mi2*t*d*eps1p1*eps1p2*tsfcusf2r
     &  + 2.d0*mi2*t*d*eps1p1*eps1p2*tsfcsgba2r
     &  - 2.d0*mi2*t*d*eps1p1*eps1p2*usfcsgba2r
     &  - 8.d0*mi2*t*d*eps1p1**2*usfcusf
     &  - 2.d0*mi2*t*d*eps1p1**2*usfcsgba2r
     &  - 8.d0*mi2*t*d*eps1p2**2*tsfctsf
     &  + 2.d0*mi2*t*d*eps1p2**2*tsfcsgba2r
     &  - 2.d0*mi2*t*eps1p1*eps1p2*tsfcsgbb2r
     &  + 2.d0*mi2*t*eps1p1*eps1p2*usfcsgbb2r
     &  + 2.d0*mi2*t*eps1p1**2*usfcsgbb2r
     &  - 2.d0*mi2*t*eps1p2**2*tsfcsgbb2r
     &  - 8.d0*mi2*t*sgbacsgba
     &  - 2.d0*mi2*d*eps1p1*eps1p2*tsfcpoint2r
     &  + 2.d0*mi2*d*eps1p1*eps1p2*usfcpoint2r
     &  + 2.d0*mi2*d*eps1p1*eps1p2*sgbacpoint2r
     &  + 2.d0*mi2*d*eps1p1**2*usfcpoint2r
     &  + mi2*d*eps1p1**2*sgbacpoint2r
     &  - 2.d0*mi2*d*eps1p2**2*tsfcpoint2r
     &  + mi2*d*eps1p2**2*sgbacpoint2r
     &  - 16.d0*mi2*eps1p1*eps1p2*sgbacsgba
     &  + 8.d0*mi2*eps1p1*eps1p2*tsfcusf2r
     &  + 4.d0*mi2*eps1p1*eps1p2*tsfcsgba2r
     &  - 12.d0*mi2*eps1p1*eps1p2*usfcsgba2r
     &  - 2.d0*mi2*eps1p1*eps1p2*sgbbcpoint2r
     &  - 8.d0*mi2*eps1p1**2*usfcusf
     &  - 2.d0*mi2*eps1p1**2*sgbacsgba
     &  - 4.d0*mi2*eps1p1**2*usfcsgba2r
     &  - mi2*eps1p1**2*sgbbcpoint2r
     &  - 8.d0*mi2*eps1p2**2*tsfctsf
     &  - 14.d0*mi2*eps1p2**2*sgbacsgba
     &  + 12.d0*mi2*eps1p2**2*tsfcsgba2r
     &  - mi2*eps1p2**2*sgbbcpoint2r
     &  - 2.d0*mi2*sgbacpoint2r
     &  + 2.d0*mi2**2*gb*eps1p1*eps1p2*sgbbcsgbb
     &  + mi2**2*gb*eps1p1**2*sgbbcsgbb
     &  + mi2**2*gb*eps1p2**2*sgbbcsgbb
     &  + 2.d0*mi2**2*gb*sgbacsgbb2r
     &  + 2.d0*mi2**2*gb**2*sgbbcsgbb
     &  + 2.d0*mi2**2*d*eps1p1*eps1p2*sgbacsgba
     &  - 4.d0*mi2**2*d*eps1p1*eps1p2*tsfcusf2r
     &  - 2.d0*mi2**2*d*eps1p1*eps1p2*tsfcsgba2r
     &  + 2.d0*mi2**2*d*eps1p1*eps1p2*usfcsgba2r
     &  + 4.d0*mi2**2*d*eps1p1**2*usfcusf
     &  + mi2**2*d*eps1p1**2*sgbacsgba
     &  + 2.d0*mi2**2*d*eps1p1**2*usfcsgba2r
     &  + 4.d0*mi2**2*d*eps1p2**2*tsfctsf
     &  + mi2**2*d*eps1p2**2*sgbacsgba
     &  - 2.d0*mi2**2*d*eps1p2**2*tsfcsgba2r
     &  + 2.d0*mi2**2*eps1p1*eps1p2*tsfcsgbb2r
     &  - 2.d0*mi2**2*eps1p1*eps1p2*usfcsgbb2r
     &  - 2.d0*mi2**2*eps1p1*eps1p2*sgbacsgbb2r
     &  - 2.d0*mi2**2*eps1p1**2*usfcsgbb2r
     &  - mi2**2*eps1p1**2*sgbacsgbb2r
     &  + 2.d0*mi2**2*eps1p2**2*tsfcsgbb2r
     &  - mi2**2*eps1p2**2*sgbacsgbb2r
     &  + 2.d0*mi2**2*sgbacsgba
     &  + 2.d0*mj2*gb*s*sgbacsgbb2r
     &  + 4.d0*mj2*gb*t*sgbacsgbb2r
     &  - 2.d0*mj2*gb*d*eps1p1*eps1p2*usfcsgba2r
     &  - 2.d0*mj2*gb*d*eps1p1**2*usfcsgba2r
     &  + 2.d0*mj2*gb*eps1p1*eps1p2*tsfcsgbb2r
     &  + 4.d0*mj2*gb*eps1p1*eps1p2*usfcsgbb2r
     &  - mj2*gb*eps1p1**2*sgbacsgbb2r
     &  - 2.d0*mj2*gb*eps1p2**2*tsfcsgbb2r
     &  + mj2*gb*eps1p2**2*sgbacsgbb2r
     &  + 4.d0*mj2*gb*sgbacsgba
     &  + 2.d0*mj2*gb*sgbbcpoint2r
     &  - 2.d0*mj2*gb**2*sgbacsgbb2r
     &  + 2.d0*mj2*s*d*eps1p1*eps1p2*usfcsgba2r
     &  + 2.d0*mj2*s*d*eps1p1**2*sgbacsgba
     &  + 2.d0*mj2*s*d*eps1p1**2*usfcsgba2r
     &  - 2.d0*mj2*s*d*eps1p2**2*sgbacsgba
     &  - 2.d0*mj2*s*eps1p1*eps1p2*usfcsgbb2r
     &  - 2.d0*mj2*s*eps1p1**2*usfcsgbb2r
     &  - mj2*s*eps1p1**2*sgbacsgbb2r
     &  + mj2*s*eps1p2**2*sgbacsgbb2r
     &  - 4.d0*mj2*s*sgbacsgba
     &  - 2.d0*mj2*t*d*eps1p1*eps1p2*tsfcsgba2r
     &  + 2.d0*mj2*t*d*eps1p1*eps1p2*usfcsgba2r
     &  + 2.d0*mj2*t*d*eps1p1**2*usfcsgba2r
     &  - 2.d0*mj2*t*d*eps1p2**2*tsfcsgba2r
     &  + 2.d0*mj2*t*eps1p1*eps1p2*tsfcsgbb2r
     &  - 2.d0*mj2*t*eps1p1*eps1p2*usfcsgbb2r
     &  - 2.d0*mj2*t*eps1p1**2*usfcsgbb2r
     &  + 2.d0*mj2*t*eps1p2**2*tsfcsgbb2r
     &  - 8.d0*mj2*t*sgbacsgba
     &  - 2.d0*mj2*d*eps1p1*eps1p2*sgbacpoint2r
     &  - mj2*d*eps1p1**2*sgbacpoint2r
     &  - mj2*d*eps1p2**2*sgbacpoint2r
     &  - 16.d0*mj2*eps1p1*eps1p2*sgbacsgba
     &  + 8.d0*mj2*eps1p1*eps1p2*tsfcusf2r
     &  + 6.d0*mj2*eps1p1*eps1p2*tsfcsgba2r
     &  - 10.d0*mj2*eps1p1*eps1p2*usfcsgba2r
     &  + 2.d0*mj2*eps1p1*eps1p2*sgbbcpoint2r
     &  - 16.d0*mj2*eps1p1**2*usfcusf
     &  - 14.d0*mj2*eps1p1**2*sgbacsgba
     &  - 14.d0*mj2*eps1p1**2*usfcsgba2r
     &  + mj2*eps1p1**2*sgbbcpoint2r
     &  - 2.d0*mj2*eps1p2**2*sgbacsgba
     &  + 2.d0*mj2*eps1p2**2*tsfcsgba2r
     &  + mj2*eps1p2**2*sgbbcpoint2r
     &  - 2.d0*mj2*sgbacpoint2r
     &  + 2.d0*mj2**2*gb*eps1p1*eps1p2*sgbbcsgbb
     &  + mj2**2*gb*eps1p1**2*sgbbcsgbb
     &  + mj2**2*gb*eps1p2**2*sgbbcsgbb
     &  - 2.d0*mj2**2*gb*sgbacsgbb2r
     &  + 2.d0*mj2**2*gb**2*sgbbcsgbb
     &  + 2.d0*mj2**2*d*eps1p1*eps1p2*sgbacsgba
     &  + mj2**2*d*eps1p1**2*sgbacsgba
     &  + mj2**2*d*eps1p2**2*sgbacsgba
     &  - 2.d0*mj2**2*eps1p1*eps1p2*sgbacsgbb2r
     &  - mj2**2*eps1p1**2*sgbacsgbb2r
     &  - mj2**2*eps1p2**2*sgbacsgbb2r
     &  + 2.d0*mj2**2*sgbacsgba
     &  + 2.d0*gb*s*d*eps1p1*eps1p2*usfcsgba2r
     &  - 8.d0*gb*s*d*eps1p1**2*usfcusf
     &  - 2.d0*gb*s*d*eps1p1**2*usfcsgba2r
     &  - 4.d0*gb*s*sgbacsgba
     &  + 4.d0*gb*t*d*eps1p1*eps1p2*tsfcusf2r
     &  - 8.d0*gb*t*d*eps1p1**2*usfcusf
     &  - 8.d0*gb*t*sgbacsgba
     &  + 2.d0*gb*d*eps1p1*eps1p2*usfcpoint2r
     &  + 2.d0*gb*d*eps1p1**2*usfcpoint2r
     &  - 2.d0*gb*eps1p1*eps1p2*sgbacsgba
     &  - 2.d0*gb*eps1p1*eps1p2*tsfcsgba2r
     &  - 4.d0*gb*eps1p1*eps1p2*usfcsgba2r
     &  - 4.d0*gb*eps1p1**2*usfcusf
     &  + gb*eps1p1**2*sgbacsgba
     &  + 4.d0*gb*eps1p2**2*tsfctsf
     &  + gb*eps1p2**2*sgbacsgba
     &  - 2.d0*gb*eps1p2**2*tsfcsgba2r
     &  - 2.d0*gb*sgbacpoint2r
     &  + 4.d0*gb**2*d*eps1p1**2*usfcusf
     &  + 2.d0*gb**2*sgbacsgba
     &  - 4.d0*s*t*d*eps1p1*eps1p2*tsfcusf2r
     &  - 2.d0*s*t*d*eps1p1*eps1p2*tsfcsgba2r
     &  - 2.d0*s*t*d*eps1p1*eps1p2*usfcsgba2r
     &  + 8.d0*s*t*d*eps1p1**2*usfcusf
     &  + 2.d0*s*t*d*eps1p1**2*usfcsgba2r
     &  + 2.d0*s*t*d*eps1p2**2*tsfcsgba2r
     &  + 8.d0*s*t*sgbacsgba
     &  - 2.d0*s*d*eps1p1*eps1p2*usfcpoint2r
     &  - 2.d0*s*d*eps1p1**2*usfcpoint2r
     &  - s*d*eps1p1**2*sgbacpoint2r
     &  + s*d*eps1p2**2*sgbacpoint2r
     &  + 12.d0*s*eps1p1*eps1p2*sgbacsgba
     &  - 4.d0*s*eps1p1*eps1p2*tsfcusf2r
     &  - 2.d0*s*eps1p1*eps1p2*tsfcsgba2r
     &  + 8.d0*s*eps1p1*eps1p2*usfcsgba2r
     &  + 8.d0*s*eps1p1**2*usfcusf
     &  + 2.d0*s*eps1p1**2*sgbacsgba
     &  + 4.d0*s*eps1p1**2*usfcsgba2r
     &  + 2.d0*s*eps1p2**2*sgbacsgba
     &  - 2.d0*s*eps1p2**2*tsfcsgba2r
     &  + 2.d0*s*sgbacpoint2r
     &  - 2.d0*s**2*d*eps1p1*eps1p2*sgbacsgba
     &  - 2.d0*s**2*d*eps1p1*eps1p2*usfcsgba2r
     &  + 4.d0*s**2*d*eps1p1**2*usfcusf
     &  + s**2*d*eps1p1**2*sgbacsgba
     &  + 2.d0*s**2*d*eps1p1**2*usfcsgba2r
     &  + s**2*d*eps1p2**2*sgbacsgba
     &  + 2.d0*s**2*sgbacsgba
     &  + 2.d0*t*d*eps1p1*eps1p2*tsfcpoint2r
     &  - 2.d0*t*d*eps1p1*eps1p2*usfcpoint2r
     &  - 2.d0*t*d*eps1p1**2*usfcpoint2r
     &  + 2.d0*t*d*eps1p2**2*tsfcpoint2r
     &  + 6.d0*t*eps1p1*eps1p2*tsfcsgba2r
     &  + 6.d0*t*eps1p1*eps1p2*usfcsgba2r
     &  + 8.d0*t*eps1p1**2*usfcusf
     &  + 2.d0*t*eps1p1**2*usfcsgba2r
     &  - 8.d0*t*eps1p2**2*tsfctsf
     &  + 2.d0*t*eps1p2**2*tsfcsgba2r
     &  + 4.d0*t*sgbacpoint2r
     &  - 4.d0*t**2*d*eps1p1*eps1p2*tsfcusf2r
     &  + 4.d0*t**2*d*eps1p1**2*usfcusf
     &  + 4.d0*t**2*d*eps1p2**2*tsfctsf
     &  + 8.d0*t**2*sgbacsgba
     &  + 2.d0*d*eps1p1*eps1p2*pointcpoint
     &  + d*eps1p1**2*pointcpoint
     &  + d*eps1p2**2*pointcpoint
     &  + 2.d0*eps1p1*eps1p2*tsfcpoint2r
     &  + 2.d0*eps1p1*eps1p2*usfcpoint2r
     &  - 2.d0*eps1p1**2*usfcpoint2r
     &  - eps1p1**2*sgbacpoint2r
     &  - 2.d0*eps1p2**2*tsfcpoint2r
     &  + eps1p2**2*sgbacpoint2r
     &  + 2.d0*pointcpoint


      if (dsasgbgb1expas.lt.0.0d0) then     
        if(aszeroprint) then
          write(*,*) ' '
          write(*,*) 'DS: negative cross section in dsasgbgb1exp:'
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
          write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
          write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
          write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
          write(*,*) 'DS: dsasgbgb1expas = ',dsasgbgb1expas
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
        dsasgbgb1expas=0.0d0
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsasgbgb1expas)*k34/(8.0d0*pi*gg1*gg2*s34*dsqrt(s))
      return
      end

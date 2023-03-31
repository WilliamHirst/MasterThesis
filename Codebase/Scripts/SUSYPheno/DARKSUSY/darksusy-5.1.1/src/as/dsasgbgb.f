c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.35, May 23, 2002, edsjo@physto.se)
c....Template file for dsasgbgb begins here

**************************************************************
*** SUBROUTINE dsasgbgb                                    ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + anti-sfermion(j)                         ***
*** -> MASSIVE gauge-boson + MASSIVE gauge-boson           ***
***                                                        ***
*** for one massive and one massless gb use dsasgbgb1exp   ***
*** for two massless gb use code dsasgbgb2exp              ***
***                                                        ***
*** The first mentioned particle (kp1) will be taken as    ***
*** a sfermion and the second particle (kp2) as an         ***
*** anti-sfermion -- not the opposite.                     ***
***                                                        ***
*** When kp1 and kp2 have different                        ***
*** weak isospin (T^3=+,-1/2), then kp1 must be an         ***
*** up-type-sfermion and kp2 a down-type-anti-sfermion.    ***
***                                                        ***
*** When one gauge boson have electric charge while the    ***
*** other is neutral, then the charged one must be         ***
*** mentioned first (kp3) and then the neutral one (kp4)   ***
*** -- not the opposite.                                   ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-23  rewritten:02-03-12                     ***
*** QCD included: 02-03-20                                 ***
*** Ghost term excluded: 02-05-22                          ***
*** rewritten: 02-07-04 (now only massive gb)              ***
*** added flavour changing charged exchange for W^-W^+:    ***
*** added by Mia Schelke 2005-06-14 + Jan 2007             ***
*** added flavour changing charged exchange for W^+Z:      ***
*** added by Mia Schelke Jan 2007                          ***
*** terms rearranged by Paolo Gondolo, 2005-06             ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, plus new factorization of first terms          ***
*** 08-05-30                                               ***
**************************************************************

      subroutine dsasgbgb(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      real*8      dsasgbgbas
      real*8      mi2,mj2,m12,m22
      real*8      mog1,mog2,smorga(2),f(2)
      real*8      cmog0,cmog1,cmog2,cmog12
      real*8      c1,c2,c3,c4,c6,c8,c11,c12,c13,c15,csum,c04,c03,d1
      complex*16  g4p
      complex*16  tsf,usf,shb,sgba,sgbb,point
      real*8      tsfctsf,usfcusf,shbcshb,sgbacsgba,
     &            sgbbcsgbb,pointcpoint,
     &            tsfcusf2r,tsfcshb2r,tsfcsgba2r,
     &            tsfcsgbb2r,tsfcpoint2r,
     &            usfcshb2r,usfcsgba2r,usfcsgbb2r,usfcpoint2r,
     &            shbcsgba2r,shbcsgbb2r,shbcpoint2r,
     &            sgbacsgbb2r,sgbacpoint2r,
     &            sgbbcpoint2r
*,ghost
      real*8      ctt,cuu,cshbshb,csgbsgb,cpp,
     &            ctu,ctshb,ctsgb,ctp,cushb,
     &            cusgb,cup,cshbsgb,cshbp,csgbp
      integer     kgb1,kgb2
      integer     k,i
      integer     ksfert(6),ksferu(6),khbs(2),kgbs(2)
      integer     nsfert,nsferu,nhbs,ngbs,npoint
      real*8 pi
      parameter (pi=3.141592653589793238d0)


c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      

***** check initial state is ok:  
c      if(abs(itype(1)-itype(2)).gt.1) then
c        write(*,*) 'DS: dsasgbgbnew called with wrong particles'   
c        write(*,*) 'DS: in the initial state : ',pname(kp1),pname(kp2)  
c        stop  
c      endif    

***** set and check particles in the final state:
      kgb1=kp3
      kgb2=kp4
c      if(kgb1.ne.kz.and.kgb1.ne.kw) then
c        write(*,*) 'DS: dsasgbgb called with kp3 = ',pname(kgb1)  
c        write(*,*) 'DS: instead of a massive gauge boson'
c        stop
c      endif 
c      if(kgb2.ne.kz.and.kgb2.ne.kw) then
c         write(*,*) 'DS: dsasgbgb called with kp4 = ',pname(kgb2)  
c         write(*,*) 'DS: instead of a massive gauge boson'
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
      m12=mass(kgb1)**2
      m22=mass(kgb2)**2
      s=Svar  
      t=Tvar  
      u=Uvar

***** now set the color factors to 1.d0 for two sleptons in
***** the initial state and to 3.d0 for two squarks in 
***** the initial state
***** if the gauge boson in the final state is a gluon, then the 
***** color factor is changed -- see case eight -> eleven below 

      ctt=1.d0
      cuu=1.d0
      cshbshb=1.d0
      csgbsgb=1.d0
      cpp=1.d0
      ctu=1.d0
      ctshb=1.d0
      ctsgb=1.d0
      ctp=1.d0
      cushb=1.d0
      cusgb=1.d0
      cup=1.d0
      cshbsgb=1.d0
      cshbp=1.d0
      csgbp=1.d0

      if(itype(1).eq.ivtype(ku).or.itype(1).eq.ivtype(kd)) then
        ctt=3.d0
        cuu=3.d0
        cshbshb=3.d0
        csgbsgb=3.d0
        cpp=3.d0
        ctu=3.d0
        ctshb=3.d0
        ctsgb=3.d0
        ctp=3.d0
        cushb=3.d0
        cusgb=3.d0
        cup=3.d0
        cshbsgb=3.d0
        cshbp=3.d0
        csgbp=3.d0
      endif

***** the ghost term is set to zero
***** it will be reintroduced for case eight below
***** 02-05-22 no, see code dsasgbgb2exp. instead
*      ghost=0.d0

***** initially setting the arrays of exchange sfermions 
***** to be empty
      do i=1,6
	ksfert(i)=0
        ksferu(i)=0
      enddo
      nsfert=0
      nsferu=0

*****
*****
******************* the first case ***************************
***** (= the second case in the original dsasgbgb code)
***** up-type-squark(i) + anti-down-type-squark(j) -> W^+ + Z
*****   with i,j = 1, 2, ... 6
***** sneutrino_i + slepton_i^+ -> W^+ + Z
*****
      if(icase.eq.1) then
        mog1=1.d0/mass(kgb1)**2
        mog2=1.d0/mass(kgb2)**2
***** in the completeness relation for the polarization vectors 
***** mog1 and mog2 are factors that we have defined for the 
***** terms that are inversely proportional to the mass squared
***** mogi is mass(kgbi)**-2 for a massive gauge boson in the 
***** final state and 0 for the photon
*****
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** in the t- and u-channel (1 up to 6 sfermions):  
        nsfert=nsfertc
        do i=1,nsfertc
          ksfert(i)=ksfertc(i)
        enddo
        nsferu=nsferuc
        do i=1,nsferuc
          ksferu(i)=ksferuc(i)
        enddo
***** no s-channel with higgs exchange:
        khbs(1)=0
        khbs(2)=0
        nhbs=0
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
***** gl(kgb1,kgb2,kgbs(k)=gl(W^+_out,Z,W^+_in)
***** =gl(W^-_in,Z,W^-_out)
***** The 'problem' is that when we write gl(kw,kz,kw), then 
***** 'darkSUSY' will instead interpret this as 
***** gl(W^-_out,Z,W^-_in), i.e. an odd permutation of the 
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
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~_d,f~_u,W^+_out,Z)
***** g4p interpretates g4p(1,2,3,4) as g4p(out,in,out,in),
***** and it takes W^+ (i.e. not W^- as above!) to give the direction of kw,
***** so we have got it correct in this case
***** but this is only because we did take the order of the 
***** gauge bosons to be (kgb1,kgb2)
***** (If we had written g4p(f~_d,f~_u,Z,W^-_in)
***** =g4p(f~_d,f~_u,kz,kw), then g4p would have interpreted
***** it as g4p(f~_d_out,f~_u_in,Z_out,W^+_in), which violates 
***** charge conservation)
***** The vertex factor we want, is constructed in g4p in
***** the following way:
***** g4p(kz,kw,f~_u,f~_d) is defined in dsvertx3.f,
***** the hermitian conjugated of this is
***** g4p(kw,kz,f~_d,f~_u), i.e. (1<->2,3<->4)
***** which is identical to what we want (1<->3)+(2<->4),  
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles   
*****
*****
        goto 500
      endif 
*****
*****
******************* the second case *************************
***** (= the fifth case in the original dsasgbgb code)
***** up-type-sfermion(i) + anti-up-type-sfermion(j) -> W^- + W^+ 
*****   with i,j = 1, 2, ... 6
***** sneutrino_i + anti-sneutrino_i -> W^- + W^+
***** we take the first gauge boson 
***** (kp3, with four-momentum k1 in the form code) 
***** to be W^- and the other one is then W^+
      if(icase.eq.2) then
        mog1=1.d0/mass(kgb1)**2
        mog2=1.d0/mass(kgb2)**2
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** no t-channel:  
        ksfert(1)=0
        ksfert(2)=0
        nsfert=0
***** in the u-channel (sfermions of down-type):  
***** for squarks in the initial state we take all
***** down-type sq in the exchange, while for 
***** sleptons we only take those sl-down of the 
***** family of the initial state 
        nsferu=nsferuc
        do i=1,nsferu
          ksferu(i)=ksferuc(i)
        enddo  
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
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
***** gl(kgb1,kgb2,kgbs(k)=gl(W^-_out,W^-_in,Z/photon), which is  
***** what 'darkSUSY' understand by gl(kw,kw,kz/kgamma), and thus
***** we just set the 'structure constants', f(k), to 1.
        kgbs(1)=kz
        smorga(1)=1.d0/mass(kz)**2
        f(1)=1.d0
        if(itype(1).eq.ivtype(ku)) then
          kgbs(2)=kgamma 
          smorga(2)=0.d0
          f(2)=1.d0
          ngbs=2
        else
          kgbs(2)=0 
          smorga(2)=0.d0
          f(2)=0.d0
          ngbs=1
        endif        
***** in the propagator, smorga(k) is a factor that we have
***** defined for the term that is inversely proportional to 
***** the mass squared
***** smorga(k) is mass(k)**-2 for a massive gauge boson in the 
***** exchange and 0 for the photon (or no exchange)
*****
***** for the point interaction
        npoint=1
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~_u,f~_u,W^-_out,W^-_in)
***** g4p interpretates g4p(1,2,3,4) as g4p(out,in,out,in),
***** and it interpretes the direction of the W boson to
***** be that of the W^+ (and not W^- !), i.e. when we write
***** g4p(f~_u,f~_u,kw,kw) it takes it to be
***** g4p(f~_u,f~_u,W^+_out,W^+_in) 
***** so we should either have interchanged kgb1 and kgb2 in g4p
*****(but this would have created problems for some of the other cases)
***** or have considered the hermitian conjugated vertex,
***** then the interpretation g4p makes would be correct
***** We don't have to do so, however, since the two vertices
***** (the one we have and the one g4p thinks we have)
***** are identical to eachother (the Lagrangian term is its 
***** own hermitian conjugate; also note that there is no 
***** momentum dependence in the vertex factor, just the 
***** symmetric metric g_{\mu\nu})      
***** 
        goto 500
      endif 
*****
*****
******************* the third case **************************
***** (= the sixth case in the original dsasgbgb code) 
***** down-type-sfermion(i) + anti-down-type-sfermion(j) -> W^- + W^+ 
*****   with i,j = 1, 2, ... 6
***** slepton_i^- + slepton_i^+ -> W^- + W^+
***** we take the first gauge boson 
***** (kp3, with four-momentum k1 in the form code) 
***** to be W^- and the other one is then W^+
      if(icase.eq.3) then 
        mog1=1.d0/mass(kgb1)**2
        mog2=1.d0/mass(kgb2)**2
*****
***** set the symmetry factor for non-identical final state particles
        s34=1.d0  
***** now set particles in the intermediate states 
***** in the t-channel(up-type sfermions): 
***** for squarks in the initial state we take all
***** up-type sq in the exchange, while for 
***** sleptons we only take the sneutrino of the 
***** family of the initial state
        nsfert=nsfertc
        do i=1,nsfert
          ksfert(i)=ksfertc(i)
        enddo
***** no u-channel:  
        nsferu=0
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
***** in the s-channel with gauge-boson exchange:
***** for comments see the previous case 
        kgbs(1)=kz
        smorga(1)=1.d0/mass(kz)**2
        f(1)=1.d0
        kgbs(2)=kgamma 
        smorga(2)=0.d0
        f(2)=1.d0
        ngbs=2
*****
***** for the point interaction
        npoint=1
***** for g4p, a comment similar to the one for
***** the previous case applies    
*****
        goto 500
      endif 
*****
*****
******************* the fourth case *************************
***** (= the seventh case in the original dsasgbgb code)
***** squark(i) + anti-squark(j) -> Z + Z
*****   with i,j = 1, 2, ... 6
***** slepton_i + slepton_i^* -> Z + Z
      if(icase.eq.4) then 
        mog1=1.d0/mass(kgb1)**2
        mog2=1.d0/mass(kgb2)**2
*****
***** set the symmetry factor for identical final state particles
        s34=2.d0  
***** now set particles in the intermediate states 
***** in the t- and u-channels (1 up to 6 sfermions): 
        nsferu=nsferun
        do i=1,nsferu
          ksferu(i)=ksferun(i)
        enddo
        nsfert=nsfertn
        do i=1,nsfert
          ksfert(i)=ksfertn(i)
        enddo
***** in the s-channel with higgs exchange:
        khbs(1)=kh1
        khbs(2)=kh2
        nhbs=2
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
***** Below, we will call the function g4p with the following argument
***** g4p(kp2,kp1,kgb1,kgb2)=g4p(f~,f~,Z,Z)
***** The vertex factor we want is constructed in g4p in
***** the following way:
***** g4p(Z,Z,f~,f~) is defined in dsvertx3.f,
***** and it is identical to the one we want, which is
***** obtained by the particle interchange (1<->3)+(2<->4)
***** i.e. interchanging the order in which the two outgoing particles
***** are mentioned, and likewise for the two incoming particles   
        goto 500
      endif
 
      write(*,*) 'DS: dsasgbgb called with wrong icase : ',icase   
      write(*,*) 'DS: initial and final states : ',  
     &             pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  


 500  continue


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

c....then the s-channel with higgs exchange
      shb=dcmplx(0.d0,0.d0)
      if(nhbs.gt.0) then
       do k=1,nhbs
        shb=shb+gl(khbs(k),kp2,kp1)
     &    *(-dsasdepro(s,khbs(k)))
     &    *mass(kw)*gl(khbs(k),kgb1,kgb2)
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
      shbcshb=cshbshb*dble(shb*conjg(shb))
      sgbacsgba=csgbsgb*dble(sgba*conjg(sgba))
      sgbbcsgbb=csgbsgb*dble(sgbb*conjg(sgbb))
      pointcpoint=cpp*dble(point*conjg(point))
      tsfcusf2r=ctu*2.d0*dble(tsf*conjg(usf))
      tsfcshb2r=ctshb*2.d0*dble(tsf*conjg(shb))
      tsfcsgba2r=ctsgb*2.d0*dble(tsf*conjg(sgba))
      tsfcsgbb2r=ctsgb*2.d0*dble(tsf*conjg(sgbb))
      tsfcpoint2r=ctp*2.d0*dble(tsf*conjg(point))
      usfcshb2r=cushb*2.d0*dble(usf*conjg(shb))
      usfcsgba2r=cusgb*2.d0*dble(usf*conjg(sgba))
      usfcsgbb2r=cusgb*2.d0*dble(usf*conjg(sgbb))
      usfcpoint2r=cup*2.d0*dble(usf*conjg(point))
      shbcsgba2r=cshbsgb*2.d0*dble(shb*conjg(sgba))
      shbcsgbb2r=cshbsgb*2.d0*dble(shb*conjg(sgbb))
      shbcpoint2r=cshbp*2.d0*dble(shb*conjg(point))
      sgbacsgbb2r=csgbsgb*2.d0*dble(sgba*conjg(sgbb))
      sgbacpoint2r=csgbp*2.d0*dble(sgba*conjg(point))
      sgbbcpoint2r=csgbp*2.d0*dble(sgbb*conjg(point))      


***** After the template file follows the form expression of the 
***** amplitude squared: dsasgbgbas

c....Template file for dsasgbgb ends here

      c1 = ( 
     &  + ((s+t)-((mi2+mj2)/2.d0+(m12+m22)))**4
     &  + (mi2-mj2)**2*(-1.d0/2.d0*(s+t)**2
     &                  +(s+t)*((m12+m22)+(mi2+mj2)/2.d0)
     &                  -(m12+m22)*((m12+m22)/2.d0+(mi2+mj2)/2.d0))
     &  - (mi2+mj2)**4/16.d0+(mi2*mj2)**2
     &  ) *usfcusf

      c2 = (
     &  + 1.d0/4.d0*s**4
     &  + s**2*(t-1.d0/2.d0*(mi2+mj2))*(t+s-1.d0/2.d0*(mi2+mj2))
     &  + s**2*((m12+m22)*(-t-1.d0/2.d0*s+1.d0/4.d0*(m12+m22))
     &           +mi2*m12+mj2*m22) 
     &  +(mi2-mj2)*(m12-m22)*s*t 
     &  + 1.d0/4.d0*(mi2-mj2)**2*(m12-m22)**2
     &  - 1.d0/2.d0*(mi2-mj2)*(m12-m22)*((mi2+mj2)+(m12+m22))*s
     &  ) *sgbacsgba

      c3 = (
     &  + 1.d0/2.d0*s*(s+2.d0*t)*(s+t)**2
     &  - 1.d0/2.d0*s*(2.d0*s+3.d0*t)*(s+t)*(mi2+mj2)
     &  - 1.d0/2.d0*s*(3.d0*s+5.d0*t)*(s+t)*(m12+m22)
     &  + 1.d0/2.d0*s*(s+t)*(mi2+mj2)**2
     &  + 1.d0/2.d0*s*(3.d0*s+4.d0*t)*(m12+m22)**2
     &  + s*mi2*mj2*(1.d0/2.d0*s+t)
     &  + 3.d0/2.d0*s*(s+t)*(mi2+mj2)*(m12+m22)
     &  + s*(s+2.d0*t)*(mi2*m12+mj2*m22)
     &  + 1.d0/2.d0*(mi2-mj2)*(m12-m22)*t*(t-((mi2+mj2)+2.d0*(m12+m22)))
     &  + 1.d0/2.d0*(mi2-mj2)*(m12-m22)
     &             *(mi2*mj2+(mi2+mj2)*(m12+m22)+(m12+m22)**2)
     &  - 1.d0/2.d0*mi2*mj2*((mi2+mj2)+3.d0*(m12+m22))*s
     &  - 2.d0*(mi2*m12+mj2*m22)*(m12+m22)*s
     &  - (mi2**2*m12+mj2**2*m22)*s
     &  - 1.d0/2.d0*(m12+m22)**3*s
     &  ) *usfcsgba2r

      c4 = (
     &  - 1.d0/2.d0*s*(s+t)**2
     &  + (m12+m22)*(1.d0/2.d0*(s+t)**2+s*(s+t))
     &  + (mi2+mj2)*(1.d0/2.d0*s*(s+t))
     &  - 1.d0/2.d0*s*(m12+m22)*((mi2+mj2)+(m12+m22))
     &  - (s+t)*(m12+m22)*(1.d0/2.d0*(mi2+mj2)+(m12+m22))
     &  + 1.d0/2.d0*(m12+m22)**2*((m12+m22)+(mi2+mj2))
     &  + 1.d0/2.d0*mi2*mj2*(m22+m12-s)
     &  ) *(usfcshb2r+usfcpoint2r)

      c6 = (
     &  s**2* (
     &  - 1.d0/2.d0*(t-m12-m22)
     &  - 1.d0/4.d0*(s-mi2-mj2))
     &  + 1.d0/2.d0*(m12+m22)*s*t
     &  + 1.d0/4.d0*(mi2-mj2)*(m12**2-m22**2)
     &  - 1.d0/2.d0*(mi2*m12+mj2*m22)*s
     &  - 1.d0/4.d0*(m12+m22)**2*s
     &  ) *(shbcsgba2r+sgbacpoint2r)

      c8 = (
     &  + 1.d0/4.d0*s**2
     &  + 1.d0/4.d0*(m12+m22)*((m12+m22)-2.d0*s)
     &  )*(shbcshb+pointcpoint+shbcpoint2r)

      c11 = (
     &  + t**2*(s+t)**2 
     &  - 2.d0*t**2*(s+t)*((mi2+mj2)+(m12+m22))
     &  - s*t*((s+t)*(mi2+mj2)-(mi2**2+mj2**2))
     &  + (2.d0*mi2*mj2+(mi2+mj2)*(m12+m22))*(t**2+2.d0*s*t+mi2*mj2
     &                                        -t*((mi2+mj2)+(m12+m22)))
     &  + t**2*((mi2+mj2)+(m12+m22))**2
     &  + mi2*mj2*(s**2-s*((mi2+mj2)+2.d0*(m12+m22))
     &             +((m12+m22)**2-mi2*mj2))
     &  ) *tsfcusf2r

      c12 = (
     &  + s*t**2*(t+1.d0/2.d0*s-3.d0/2.d0*(mi2+mj2)-1.d0/2.d0*(m12+m22))
     &  + s**2*(-1.d0/2.d0*t*(mi2+mj2)+1.d0/2.d0*mi2*mj2) 
     &  + 1.d0/2.d0*t**2*((mi2-mj2)*(m12-m22))
     &  + 1.d0/2.d0*s*t*((mi2+mj2)**2+2.d0*mi2*mj2+(mi2+mj2)*(m12+m22))
     &  - 1.d0/2.d0*s*mi2*mj2*((mi2+mj2)+(m12+m22))
     &  - 1.d0/2.d0*t*(m12-m22)*(mi2**2-mj2**2)
     &  + 1.d0/2.d0*mi2*mj2*(mi2-mj2)*(m12-m22)
     &  ) *tsfcsgba2r

      c13 = (
     &  + 1.d0/2.d0*mi2*mj2*((m12+m22)-s)
     &  - 1.d0/2.d0*(mi2+mj2)*(m12+m22)*t
     &  + 1.d0/2.d0*((m12+m22)-s)*t**2
     &  + 1.d0/2.d0*(mi2+mj2)*s*t
     &  ) *(tsfcshb2r+tsfcpoint2r)

      c15 = (
     &  + mi2**2*mj2**2
     &  - 2.d0*mi2*mj2*(mi2+mj2-t)*t
     &  + (mi2+mj2)**2*t**2
     &  + t**4
     &  - 2.d0*(mi2+mj2)*t**3
     &  ) *tsfctsf

      csum=c1+c2+c3+c4+c6+c8+c11+c12+c13+c15


      cmog12 = mog1*mog2*csum
      
      cmog1 = mog1*(
     &  + mi2*mj2*m12*m22*s*sgbbcsgbb
     &  - 1.d0/2.d0*mi2*mj2*m12*m22*tsfcsgbb2r
     &  + 5.d0/2.d0*mi2*mj2*m12*m22*usfcsgbb2r
     &  + 3.d0/2.d0*mi2*mj2*m12*m22*sgbacsgbb2r
     &  + mi2*mj2*m12*m22**2*sgbbcsgbb
     &  + 1.d0/2.d0*mi2*mj2*m12*s*sgbacsgbb2r
     &  - 4.d0*mi2*mj2*m12*usfcusf
     &  + mi2*mj2*m12*sgbacsgba
     &  - 2.d0*mi2*mj2*m12*tsfcusf2r
     &  + 4.d0*mi2*mj2*m12*usfcsgba2r
     &  - 1.d0/2.d0*mi2*mj2*m12**2*m22*sgbbcsgbb
     &  - 1.d0/2.d0*mi2*mj2*m12**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m22*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mi2*mj2*m22*s*usfcsgbb2r
     &  - mi2*mj2*m22*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m22*s**2*sgbbcsgbb
     &  - 2.d0*mi2*mj2*m22*t*tsfcsgbb2r
     &  - 2.d0*mi2*mj2*m22*t*usfcsgbb2r
     &  - 4.d0*mi2*mj2*m22*usfcusf
     &  + 3.d0*mi2*mj2*m22*sgbacsgba
     &  - 2.d0*mi2*mj2*m22*tsfcusf2r
     &  - mi2*mj2*m22*tsfcsgba2r
     &  + mi2*mj2*m22*usfcsgba2r
     &  + mi2*mj2*m22**2*s*sgbbcsgbb
     &  + 1.d0/2.d0*mi2*mj2*m22**2*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2*mj2*m22**2*usfcsgbb2r
     &  + mi2*mj2*m22**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m22**3*sgbbcsgbb
     &  + 4.d0*mi2*mj2*s*usfcusf
     &  - 3.d0*mi2*mj2*s*sgbacsgba
     &  + 2.d0*mi2*mj2*s*tsfcusf2r
     &  + mi2*mj2*s*tsfcsgba2r
     &  - mi2*mj2*s*usfcsgba2r
     &  + 4.d0*mi2*mj2*t*tsfctsf
     &  + 4.d0*mi2*mj2*t*usfcusf
     &  + 4.d0*mi2*mj2*t*tsfcusf2r
     &  - 4.d0*mi2*mj2*t*tsfcsgba2r
     &  - 4.d0*mi2*mj2*t*usfcsgba2r
     &  - 1.d0/2.d0*mi2*mj2**2*m12*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*mj2**2*m12*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*mj2**2*m22*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*mj2**2*m22*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*mj2**2*s*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*mj2**2*s*usfcsgbb2r
     &  - mi2*m12*m22*s*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*m22*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*m12*m22*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*m22*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*m22*t*sgbacsgbb2r
     &  - 14.d0*mi2*m12*m22*usfcusf
     &  + mi2*m12*m22*tsfcusf2r
     &  + 1.d0/2.d0*mi2*m12*m22*tsfcsgba2r
     &  - 5.d0/2.d0*mi2*m12*m22*usfcsgba2r
     &  - 1.d0/2.d0*mi2*m12*m22*shbcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*m22*sgbbcpoint2r
     &  - 1.d0/2.d0*mi2*m12*m22**2*usfcsgbb2r
     &  + 2.d0*mi2*m12*s*t*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*s*t*sgbacsgbb2r
     &  + 16.d0*mi2*m12*s*usfcusf
     &  - mi2*m12*s*sgbacsgba
     &  - 1.d0/2.d0*mi2*m12*s*tsfcsgba2r
     &  + 5.d0/2.d0*mi2*m12*s*usfcsgba2r
     &  - 1.d0/2.d0*mi2*m12*s*shbcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*s*sgbbcpoint2r
     &  + 3.d0/2.d0*mi2*m12*s**2*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*s**2*sgbacsgbb2r
     &  + 16.d0*mi2*m12*t*usfcusf
     &  - mi2*m12*t*sgbacsgba
     &  - mi2*m12*t*tsfcsgba2r
     &  - mi2*m12*t*usfcsgba2r
     &  + 1.d0/2.d0*mi2*m12*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*tsfcshb2r
     &  + 1.d0/2.d0*mi2*m12*tsfcpoint2r
     &  - 5.d0/2.d0*mi2*m12*usfcshb2r
     &  - 5.d0/2.d0*mi2*m12*usfcpoint2r
     &  - 1.d0/2.d0*mi2*m12*shbcsgba2r
     &  - 1.d0/2.d0*mi2*m12*sgbacpoint2r
     &  + 1.d0/4.d0*mi2*m12**2*m22*sgbacsgbb2r
     &  - 3.d0/2.d0*mi2*m12**2*s*usfcsgbb2r
     &  - 1.d0/4.d0*mi2*m12**2*s*sgbacsgbb2r
     &  - mi2*m12**2*t*usfcsgbb2r
     &  - 8.d0*mi2*m12**2*usfcusf
     &  + 1.d0/2.d0*mi2*m12**2*sgbacsgba
     &  + mi2*m12**2*tsfcusf2r
     &  + 1.d0/4.d0*mi2*m12**2*shbcsgbb2r
     &  + 1.d0/4.d0*mi2*m12**2*sgbbcpoint2r
     &  + 1.d0/2.d0*mi2*m12**3*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*s*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m22*s*t*usfcsgbb2r
     &  + 14.d0*mi2*m22*s*usfcusf
     &  + 4.d0*mi2*m22*s*sgbacsgba
     &  + mi2*m22*s*tsfcusf2r
     &  + 6.d0*mi2*m22*s*usfcsgba2r
     &  - 1.d0/2.d0*mi2*m22*s*shbcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*s*sgbbcpoint2r
     &  + mi2*m22*s**2*usfcsgbb2r
     &  + 1.d0/4.d0*mi2*m22*s**2*sgbacsgbb2r
     &  - 2.d0*mi2*m22*t*tsfctsf
     &  + 14.d0*mi2*m22*t*usfcusf
     &  + mi2*m22*t*sgbacsgba
     &  + 2.d0*mi2*m22*t*tsfcusf2r
     &  + 2.d0*mi2*m22*t*usfcsgba2r
     &  - 1.d0/2.d0*mi2*m22*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*t**2*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*tsfcshb2r
     &  - 1.d0/2.d0*mi2*m22*tsfcpoint2r
     &  - 3.d0/2.d0*mi2*m22*usfcshb2r
     &  - 3.d0/2.d0*mi2*m22*usfcpoint2r
     &  - 1.d0/2.d0*mi2*m22*shbcsgba2r
     &  - 1.d0/2.d0*mi2*m22*sgbacpoint2r
     &  - 1.d0/2.d0*mi2*m22**2*s*usfcsgbb2r
     &  + 1.d0/4.d0*mi2*m22**2*s*sgbacsgbb2r
     &  + 1.d0/2.d0*mi2*m22**2*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m22**2*t*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m22**2*t*sgbacsgbb2r
     &  - 6.d0*mi2*m22**2*usfcusf
     &  - 5.d0/2.d0*mi2*m22**2*sgbacsgba
     &  + 1.d0/2.d0*mi2*m22**2*tsfcsgba2r
     &  - 5.d0/2.d0*mi2*m22**2*usfcsgba2r
     &  + 1.d0/4.d0*mi2*m22**2*shbcsgbb2r
     &  + 1.d0/4.d0*mi2*m22**2*sgbbcpoint2r
     &  - 1.d0/4.d0*mi2*m22**3*sgbacsgbb2r
     &  - 16.d0*mi2*s*t*usfcusf
     &  + mi2*s*t*sgbacsgba
     &  - 4.d0*mi2*s*t*tsfcusf2r
     &  - 2.d0*mi2*s*t*usfcsgba2r
     &  - 1.d0/2.d0*mi2*s*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*s*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*s*tsfcshb2r
     &  + 1.d0/2.d0*mi2*s*tsfcpoint2r
     &  + 3.d0/2.d0*mi2*s*usfcshb2r
     &  + 3.d0/2.d0*mi2*s*usfcpoint2r
     &  + 1.d0/2.d0*mi2*s*shbcsgba2r
     &  + 1.d0/2.d0*mi2*s*sgbacpoint2r
     &  - mi2*s**2*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*s**2*t*sgbacsgbb2r
     &  - 8.d0*mi2*s**2*usfcusf
     &  - 3.d0/2.d0*mi2*s**2*sgbacsgba
     &  - mi2*s**2*tsfcusf2r
     &  - 1.d0/2.d0*mi2*s**2*tsfcsgba2r
     &  - 7.d0/2.d0*mi2*s**2*usfcsgba2r
     &  + 1.d0/4.d0*mi2*s**2*shbcsgbb2r
     &  + 1.d0/4.d0*mi2*s**2*sgbbcpoint2r
     &  - 1.d0/2.d0*mi2*s**3*usfcsgbb2r
     &  - 1.d0/4.d0*mi2*s**3*sgbacsgbb2r
     &  + 2.d0*mi2*t*tsfcshb2r
     &  + 2.d0*mi2*t*tsfcpoint2r
     &  + 2.d0*mi2*t*usfcshb2r
     &  + 2.d0*mi2*t*usfcpoint2r
     &  + 4.d0*mi2*t**2*tsfctsf
     &  - 8.d0*mi2*t**2*usfcusf
     &  - 2.d0*mi2*t**2*tsfcusf2r
     &  + mi2*t**2*tsfcsgba2r
     &  + mi2*t**2*usfcsgba2r
     &  + 1.d0/2.d0*mi2**2*mj2*m12*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*mj2*m12*usfcsgbb2r
     &  + 3.d0/2.d0*mi2**2*mj2*m22*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2**2*mj2*m22*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*mj2*s*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*mj2*s*usfcsgbb2r
     &  - 2.d0*mi2**2*mj2*tsfctsf
     &  - 2.d0*mi2**2*mj2*usfcusf
     &  - 2.d0*mi2**2*mj2*tsfcusf2r
     &  + 2.d0*mi2**2*mj2*tsfcsgba2r
     &  + 2.d0*mi2**2*mj2*usfcsgba2r
     &  - 1.d0/2.d0*mi2**2*m12*m22*s*sgbbcsgbb
     &  + 1.d0/2.d0*mi2**2*m12*m22*tsfcsgbb2r
     &  - 3.d0/2.d0*mi2**2*m12*m22*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*m22*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*m22**2*sgbbcsgbb
     &  - mi2**2*m12*s*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*t*usfcsgbb2r
     &  - 10.d0*mi2**2*m12*usfcusf
     &  - mi2**2*m12*tsfcusf2r
     &  + 1.d0/2.d0*mi2**2*m12*tsfcsgba2r
     &  - 3.d0/2.d0*mi2**2*m12*usfcsgba2r
     &  + 1.d0/4.d0*mi2**2*m12**2*m22*sgbbcsgbb
     &  + 1.d0/2.d0*mi2**2*m12**2*usfcsgbb2r
     &  + 1.d0/4.d0*mi2**2*m12**2*sgbacsgbb2r
     &  + 1.d0/2.d0*mi2**2*m22*s*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m22*s*usfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m22*s*sgbacsgbb2r
     &  + 1.d0/4.d0*mi2**2*m22*s**2*sgbbcsgbb
     &  + 3.d0/2.d0*mi2**2*m22*t*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2**2*m22*t*usfcsgbb2r
     &  + mi2**2*m22*tsfctsf
     &  - 9.d0*mi2**2*m22*usfcusf
     &  - 2.d0*mi2**2*m22*sgbacsgba
     &  - 2.d0*mi2**2*m22*tsfcusf2r
     &  - 3.d0/2.d0*mi2**2*m22*tsfcsgba2r
     &  - 7.d0/2.d0*mi2**2*m22*usfcsgba2r
     &  - 1.d0/2.d0*mi2**2*m22**2*s*sgbbcsgbb
     &  - 1.d0/2.d0*mi2**2*m22**2*tsfcsgbb2r
     &  - mi2**2*m22**2*usfcsgbb2r
     &  - 3.d0/4.d0*mi2**2*m22**2*sgbacsgbb2r
     &  + 1.d0/4.d0*mi2**2*m22**3*sgbbcsgbb
     &  + 1.d0/2.d0*mi2**2*s*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*s*t*usfcsgbb2r
     &  + 10.d0*mi2**2*s*usfcusf
     &  + mi2**2*s*sgbacsgba
     &  + 3.d0*mi2**2*s*tsfcusf2r
     &  + 3.d0/2.d0*mi2**2*s*tsfcsgba2r
     &  + 7.d0/2.d0*mi2**2*s*usfcsgba2r
     &  + 1.d0/2.d0*mi2**2*s**2*usfcsgbb2r
     &  + 1.d0/4.d0*mi2**2*s**2*sgbacsgbb2r
     &  - 2.d0*mi2**2*t*tsfctsf
     &  + 10.d0*mi2**2*t*usfcusf
     &  + 4.d0*mi2**2*t*tsfcusf2r
     &  + mi2**2*t*tsfcsgba2r
     &  + mi2**2*t*usfcsgba2r
     &  - mi2**2*tsfcshb2r
     &  - mi2**2*tsfcpoint2r
     &  - mi2**2*usfcshb2r
     &  - mi2**2*usfcpoint2r
     &  - mi2**3*m22*tsfcsgbb2r
     &  - mi2**3*m22*usfcsgbb2r
     &  - 4.d0*mi2**3*usfcusf
     &  - 2.d0*mi2**3*tsfcusf2r
     &  - mi2**3*tsfcsgba2r
     &  - mi2**3*usfcsgba2r
     &  + mj2*m12*m22*s*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*m22*s*sgbacsgbb2r
     &  + 1.d0/2.d0*mj2*m12*m22*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*m22*t*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*m22*t*sgbacsgbb2r
     &  - 4.d0*mj2*m12*m22*usfcusf
     &  + mj2*m12*m22*sgbacsgba
     &  + mj2*m12*m22*usfcsgba2r
     &  + 1.d0/2.d0*mj2*m12*m22*shbcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*m22*sgbbcpoint2r
     &  + 1.d0/2.d0*mj2*m12*m22**2*usfcsgbb2r
     &  - 2.d0*mj2*m12*s*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*s*t*sgbacsgbb2r
     &  + 4.d0*mj2*m12*s*usfcusf
     &  - mj2*m12*s*sgbacsgba
     &  - mj2*m12*s*usfcsgba2r
     &  + 1.d0/2.d0*mj2*m12*s*shbcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*s*sgbbcpoint2r
     &  - 3.d0/2.d0*mj2*m12*s**2*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*s**2*sgbacsgbb2r
     &  + 4.d0*mj2*m12*t*usfcusf
     &  + mj2*m12*t*sgbacsgba
     &  + 2.d0*mj2*m12*t*tsfcusf2r
     &  - 4.d0*mj2*m12*t*usfcsgba2r
     &  - 1.d0/2.d0*mj2*m12*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*t**2*usfcsgbb2r
     &  + mj2*m12*shbcsgba2r
     &  + mj2*m12*sgbacpoint2r
     &  - 1.d0/4.d0*mj2*m12**2*m22*sgbacsgbb2r
     &  + 3.d0/2.d0*mj2*m12**2*s*usfcsgbb2r
     &  + 1.d0/4.d0*mj2*m12**2*s*sgbacsgbb2r
     &  + mj2*m12**2*t*usfcsgbb2r
     &  - 2.d0*mj2*m12**2*usfcusf
     &  - 1.d0/2.d0*mj2*m12**2*sgbacsgba
     &  + 2.d0*mj2*m12**2*usfcsgba2r
     &  - 1.d0/4.d0*mj2*m12**2*shbcsgbb2r
     &  - 1.d0/4.d0*mj2*m12**2*sgbbcpoint2r
     &  - 1.d0/2.d0*mj2*m12**3*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*s*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m22*s*t*usfcsgbb2r
     &  + 4.d0*mj2*m22*s*usfcusf
     &  + mj2*m22*s*sgbacsgba
     &  + 2.d0*mj2*m22*s*usfcsgba2r
     &  + 1.d0/2.d0*mj2*m22*s*shbcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*s*sgbbcpoint2r
     &  - mj2*m22*s**2*usfcsgbb2r
     &  - 1.d0/4.d0*mj2*m22*s**2*sgbacsgbb2r
     &  + 4.d0*mj2*m22*t*usfcusf
     &  - 3.d0*mj2*m22*t*sgbacsgba
     &  + 2.d0*mj2*m22*t*tsfcusf2r
     &  + mj2*m22*t*tsfcsgba2r
     &  - mj2*m22*t*usfcsgba2r
     &  + 1.d0/2.d0*mj2*m22*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22**2*s*usfcsgbb2r
     &  - 1.d0/4.d0*mj2*m22**2*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2*m22**2*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m22**2*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m22**2*t*sgbacsgbb2r
     &  - 2.d0*mj2*m22**2*usfcusf
     &  - 1.d0/2.d0*mj2*m22**2*sgbacsgba
     &  - mj2*m22**2*usfcsgba2r
     &  - 1.d0/4.d0*mj2*m22**2*shbcsgbb2r
     &  - 1.d0/4.d0*mj2*m22**2*sgbbcpoint2r
     &  + 1.d0/4.d0*mj2*m22**3*sgbacsgbb2r
     &  - 4.d0*mj2*s*t*usfcusf
     &  + 3.d0*mj2*s*t*sgbacsgba
     &  - 2.d0*mj2*s*t*tsfcusf2r
     &  - mj2*s*t*tsfcsgba2r
     &  + mj2*s*t*usfcsgba2r
     &  + 1.d0/2.d0*mj2*s*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*s*t**2*usfcsgbb2r
     &  + mj2*s**2*t*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*s**2*t*sgbacsgbb2r
     &  - 2.d0*mj2*s**2*usfcusf
     &  - 1.d0/2.d0*mj2*s**2*sgbacsgba
     &  - mj2*s**2*usfcsgba2r
     &  - 1.d0/4.d0*mj2*s**2*shbcsgbb2r
     &  - 1.d0/4.d0*mj2*s**2*sgbbcpoint2r
     &  + 1.d0/2.d0*mj2*s**3*usfcsgbb2r
     &  + 1.d0/4.d0*mj2*s**3*sgbacsgbb2r
     &  - 2.d0*mj2*t**2*tsfctsf
     &  - 2.d0*mj2*t**2*usfcusf
     &  - 2.d0*mj2*t**2*tsfcusf2r
     &  + 2.d0*mj2*t**2*tsfcsgba2r
     &  + 2.d0*mj2*t**2*usfcsgba2r
     &  - 1.d0/2.d0*mj2**2*m12*m22*s*sgbbcsgbb
     &  - mj2**2*m12*m22*usfcsgbb2r
     &  - mj2**2*m12*m22*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2**2*m12*m22**2*sgbbcsgbb
     &  + mj2**2*m12*s*usfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m12*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m12*t*usfcsgbb2r
     &  - mj2**2*m12*sgbacsgba
     &  + 1.d0/4.d0*mj2**2*m12**2*m22*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*m12**2*usfcsgbb2r
     &  + 1.d0/4.d0*mj2**2*m12**2*sgbacsgbb2r
     &  + mj2**2*m22*s*usfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m22*s*sgbacsgbb2r
     &  + 1.d0/4.d0*mj2**2*m22*s**2*sgbbcsgbb
     &  + 1.d0/2.d0*mj2**2*m22*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m22*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*m22**2*s*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*m22**2*usfcsgbb2r
     &  - 1.d0/4.d0*mj2**2*m22**2*sgbacsgbb2r
     &  + 1.d0/4.d0*mj2**2*m22**3*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*s*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*s*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*s**2*usfcsgbb2r
     &  - 1.d0/4.d0*mj2**2*s**2*sgbacsgbb2r
     &  + 10.d0*m12*m22*s*usfcusf
     &  - 3.d0/2.d0*m12*m22*s*sgbacsgba
     &  + 10.d0*m12*m22*t*usfcusf
     &  - m12*m22*t*sgbacsgba
     &  - m12*m22*t*tsfcusf2r
     &  - 1.d0/2.d0*m12*m22*t*tsfcsgba2r
     &  - 5.d0/2.d0*m12*m22*t*usfcsgba2r
     &  - 2.d0*m12*m22*usfcshb2r
     &  - 2.d0*m12*m22*usfcpoint2r
     &  - 4.d0*m12*m22**2*usfcusf
     &  + 1.d0/2.d0*m12*m22**2*sgbacsgba
     &  + 1.d0/2.d0*m12*m22**2*usfcsgba2r
     &  - 12.d0*m12*s*t*usfcusf
     &  + 2.d0*m12*s*t*sgbacsgba
     &  + 1.d0/2.d0*m12*s*t*tsfcsgba2r
     &  + 5.d0/2.d0*m12*s*t*usfcsgba2r
     &  + 2.d0*m12*s*usfcshb2r
     &  + 2.d0*m12*s*usfcpoint2r
     &  - 6.d0*m12*s**2*usfcusf
     &  + m12*s**2*sgbacsgba
     &  - 1.d0/2.d0*m12*s**2*usfcsgba2r
     &  - 1.d0/2.d0*m12*t*tsfcshb2r
     &  - 1.d0/2.d0*m12*t*tsfcpoint2r
     &  + 5.d0/2.d0*m12*t*usfcshb2r
     &  + 5.d0/2.d0*m12*t*usfcpoint2r
     &  - 1.d0/2.d0*m12*t*shbcsgba2r
     &  - 1.d0/2.d0*m12*t*sgbacpoint2r
     &  - 6.d0*m12*t**2*usfcusf
     &  + m12*t**2*tsfcusf2r
     &  + 1.d0/2.d0*m12*t**2*tsfcsgba2r
     &  + 5.d0/2.d0*m12*t**2*usfcsgba2r
     &  - m12*shbcshb
     &  - m12*pointcpoint
     &  - m12*shbcpoint2r
     &  - 5.d0*m12**2*m22*usfcusf
     &  + 1.d0/4.d0*m12**2*m22*sgbacsgba
     &  + m12**2*m22*usfcsgba2r
     &  + 6.d0*m12**2*s*usfcusf
     &  - 1.d0/2.d0*m12**2*s*sgbacsgba
     &  - m12**2*s*usfcsgba2r
     &  + 6.d0*m12**2*t*usfcusf
     &  - m12**2*t*tsfcusf2r
     &  - 2.d0*m12**2*t*usfcsgba2r
     &  - 3.d0/2.d0*m12**2*usfcshb2r
     &  - 3.d0/2.d0*m12**2*usfcpoint2r
     &  + 1.d0/4.d0*m12**2*shbcsgba2r
     &  + 1.d0/4.d0*m12**2*sgbacpoint2r
     &  - 2.d0*m12**3*usfcusf
     &  + 1.d0/2.d0*m12**3*usfcsgba2r
     &  - 10.d0*m22*s*t*usfcusf
     &  + 3.d0*m22*s*t*sgbacsgba
     &  - m22*s*t*tsfcusf2r
     &  + m22*s*usfcshb2r
     &  + m22*s*usfcpoint2r
     &  + 1.d0/2.d0*m22*s*shbcsgba2r
     &  + 1.d0/2.d0*m22*s*sgbacpoint2r
     &  - 5.d0*m22*s**2*usfcusf
     &  - 3.d0/4.d0*m22*s**2*sgbacsgba
     &  - 2.d0*m22*s**2*usfcsgba2r
     &  + 1.d0/2.d0*m22*t*tsfcshb2r
     &  + 1.d0/2.d0*m22*t*tsfcpoint2r
     &  + 3.d0/2.d0*m22*t*usfcshb2r
     &  + 3.d0/2.d0*m22*t*usfcpoint2r
     &  + 1.d0/2.d0*m22*t*shbcsgba2r
     &  + 1.d0/2.d0*m22*t*sgbacpoint2r
     &  + m22*t**2*tsfctsf
     &  - 5.d0*m22*t**2*usfcusf
     &  + m22*t**2*sgbacsgba
     &  + 3.d0/2.d0*m22*t**2*tsfcsgba2r
     &  + 3.d0/2.d0*m22*t**2*usfcsgba2r
     &  + 4.d0*m22**2*s*usfcusf
     &  + m22**2*s*usfcsgba2r
     &  + 4.d0*m22**2*t*usfcusf
     &  - m22**2*t*sgbacsgba
     &  - 1.d0/2.d0*m22**2*t*tsfcsgba2r
     &  - 1.d0/2.d0*m22**2*t*usfcsgba2r
     &  - 1.d0/2.d0*m22**2*usfcshb2r
     &  - 1.d0/2.d0*m22**2*usfcpoint2r
     &  - 1.d0/4.d0*m22**2*shbcsgba2r
     &  - 1.d0/4.d0*m22**2*sgbacpoint2r
     &  - m22**3*usfcusf
     &  + 1.d0/4.d0*m22**3*sgbacsgba
     &  - 1.d0/2.d0*s*t*tsfcshb2r
     &  - 1.d0/2.d0*s*t*tsfcpoint2r
     &  - 3.d0/2.d0*s*t*usfcshb2r
     &  - 3.d0/2.d0*s*t*usfcpoint2r
     &  - 1.d0/2.d0*s*t*shbcsgba2r
     &  - 1.d0/2.d0*s*t*sgbacpoint2r
     &  + 6.d0*s*t**2*usfcusf
     &  - 2.d0*s*t**2*sgbacsgba
     &  + s*t**2*tsfcusf2r
     &  - 3.d0/2.d0*s*t**2*tsfcsgba2r
     &  - 3.d0/2.d0*s*t**2*usfcsgba2r
     &  + 6.d0*s**2*t*usfcusf
     &  - 2.d0*s**2*t*sgbacsgba
     &  + s**2*t*tsfcusf2r
     &  + 1.d0/2.d0*s**2*t*tsfcsgba2r
     &  + 1.d0/2.d0*s**2*t*usfcsgba2r
     &  - 1.d0/2.d0*s**2*usfcshb2r
     &  - 1.d0/2.d0*s**2*usfcpoint2r
     &  - 1.d0/4.d0*s**2*shbcsgba2r
     &  - 1.d0/4.d0*s**2*sgbacpoint2r
     &  + 2.d0*s**3*usfcusf
     &  + 1.d0/2.d0*s**3*sgbacsgba
     &  + s**3*usfcsgba2r
     &  - t**2*tsfcshb2r
     &  - t**2*tsfcpoint2r
     &  - t**2*usfcshb2r
     &  - t**2*usfcpoint2r
     &  - 2.d0*t**3*tsfctsf
     &  + 2.d0*t**3*usfcusf
     &  - t**3*tsfcsgba2r
     &  - t**3*usfcsgba2r
     &  )
      cmog2 = mog2*(
     &  + mi2*mj2*m12*m22*s*sgbbcsgbb
     &  - 1.d0/2.d0*mi2*mj2*m12*m22*tsfcsgbb2r
     &  + 5.d0/2.d0*mi2*mj2*m12*m22*usfcsgbb2r
     &  + 3.d0/2.d0*mi2*mj2*m12*m22*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m12*m22**2*sgbbcsgbb
     &  - 1.d0/2.d0*mi2*mj2*m12*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mi2*mj2*m12*s*usfcsgbb2r
     &  - mi2*mj2*m12*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m12*s**2*sgbbcsgbb
     &  - 2.d0*mi2*mj2*m12*t*tsfcsgbb2r
     &  - 2.d0*mi2*mj2*m12*t*usfcsgbb2r
     &  - 4.d0*mi2*mj2*m12*usfcusf
     &  + 3.d0*mi2*mj2*m12*sgbacsgba
     &  - 2.d0*mi2*mj2*m12*tsfcusf2r
     &  - mi2*mj2*m12*tsfcsgba2r
     &  + mi2*mj2*m12*usfcsgba2r
     &  + mi2*mj2*m12**2*m22*sgbbcsgbb
     &  + mi2*mj2*m12**2*s*sgbbcsgbb
     &  + 1.d0/2.d0*mi2*mj2*m12**2*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2*mj2*m12**2*usfcsgbb2r
     &  + mi2*mj2*m12**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*mj2*m12**3*sgbbcsgbb
     &  + 1.d0/2.d0*mi2*mj2*m22*s*sgbacsgbb2r
     &  - 4.d0*mi2*mj2*m22*usfcusf
     &  + mi2*mj2*m22*sgbacsgba
     &  - 2.d0*mi2*mj2*m22*tsfcusf2r
     &  + 4.d0*mi2*mj2*m22*usfcsgba2r
     &  - 1.d0/2.d0*mi2*mj2*m22**2*sgbacsgbb2r
     &  + 4.d0*mi2*mj2*s*usfcusf
     &  - 3.d0*mi2*mj2*s*sgbacsgba
     &  + 2.d0*mi2*mj2*s*tsfcusf2r
     &  + mi2*mj2*s*tsfcsgba2r
     &  - mi2*mj2*s*usfcsgba2r
     &  + 4.d0*mi2*mj2*t*tsfctsf
     &  + 4.d0*mi2*mj2*t*usfcusf
     &  + 4.d0*mi2*mj2*t*tsfcusf2r
     &  - 4.d0*mi2*mj2*t*tsfcsgba2r
     &  - 4.d0*mi2*mj2*t*usfcsgba2r
     &  + 3.d0/2.d0*mi2*mj2**2*m12*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2*mj2**2*m12*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*mj2**2*m22*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*mj2**2*m22*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*mj2**2*s*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*mj2**2*s*usfcsgbb2r
     &  - 2.d0*mi2*mj2**2*tsfctsf
     &  - 2.d0*mi2*mj2**2*usfcusf
     &  - 2.d0*mi2*mj2**2*tsfcusf2r
     &  + 2.d0*mi2*mj2**2*tsfcsgba2r
     &  + 2.d0*mi2*mj2**2*usfcsgba2r
     &  + mi2*m12*m22*s*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*m22*s*sgbacsgbb2r
     &  + 1.d0/2.d0*mi2*m12*m22*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*m22*t*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*m22*t*sgbacsgbb2r
     &  - 4.d0*mi2*m12*m22*usfcusf
     &  + mi2*m12*m22*sgbacsgba
     &  + mi2*m12*m22*usfcsgba2r
     &  + 1.d0/2.d0*mi2*m12*m22*shbcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*m22*sgbbcpoint2r
     &  - 1.d0/4.d0*mi2*m12*m22**2*sgbacsgbb2r
     &  + 1.d0/2.d0*mi2*m12*s*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12*s*t*usfcsgbb2r
     &  + 4.d0*mi2*m12*s*usfcusf
     &  + mi2*m12*s*sgbacsgba
     &  + 2.d0*mi2*m12*s*usfcsgba2r
     &  + 1.d0/2.d0*mi2*m12*s*shbcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*s*sgbbcpoint2r
     &  - mi2*m12*s**2*usfcsgbb2r
     &  - 1.d0/4.d0*mi2*m12*s**2*sgbacsgbb2r
     &  + 4.d0*mi2*m12*t*usfcusf
     &  - 3.d0*mi2*m12*t*sgbacsgba
     &  + 2.d0*mi2*m12*t*tsfcusf2r
     &  + mi2*m12*t*tsfcsgba2r
     &  - mi2*m12*t*usfcsgba2r
     &  + 1.d0/2.d0*mi2*m12*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12**2*m22*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12**2*s*usfcsgbb2r
     &  - 1.d0/4.d0*mi2*m12**2*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*m12**2*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12**2*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12**2*t*sgbacsgbb2r
     &  - 2.d0*mi2*m12**2*usfcusf
     &  - 1.d0/2.d0*mi2*m12**2*sgbacsgba
     &  - mi2*m12**2*usfcsgba2r
     &  - 1.d0/4.d0*mi2*m12**2*shbcsgbb2r
     &  - 1.d0/4.d0*mi2*m12**2*sgbbcpoint2r
     &  + 1.d0/4.d0*mi2*m12**3*sgbacsgbb2r
     &  - 2.d0*mi2*m22*s*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*s*t*sgbacsgbb2r
     &  + 4.d0*mi2*m22*s*usfcusf
     &  - mi2*m22*s*sgbacsgba
     &  - mi2*m22*s*usfcsgba2r
     &  + 1.d0/2.d0*mi2*m22*s*shbcsgbb2r
     &  + 1.d0/2.d0*mi2*m22*s*sgbbcpoint2r
     &  - 3.d0/2.d0*mi2*m22*s**2*usfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*s**2*sgbacsgbb2r
     &  + 4.d0*mi2*m22*t*usfcusf
     &  + mi2*m22*t*sgbacsgba
     &  + 2.d0*mi2*m22*t*tsfcusf2r
     &  - 4.d0*mi2*m22*t*usfcsgba2r
     &  - 1.d0/2.d0*mi2*m22*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m22*t**2*usfcsgbb2r
     &  + mi2*m22*shbcsgba2r
     &  + mi2*m22*sgbacpoint2r
     &  + 3.d0/2.d0*mi2*m22**2*s*usfcsgbb2r
     &  + 1.d0/4.d0*mi2*m22**2*s*sgbacsgbb2r
     &  + mi2*m22**2*t*usfcsgbb2r
     &  - 2.d0*mi2*m22**2*usfcusf
     &  - 1.d0/2.d0*mi2*m22**2*sgbacsgba
     &  + 2.d0*mi2*m22**2*usfcsgba2r
     &  - 1.d0/4.d0*mi2*m22**2*shbcsgbb2r
     &  - 1.d0/4.d0*mi2*m22**2*sgbbcpoint2r
     &  - 1.d0/2.d0*mi2*m22**3*usfcsgbb2r
     &  - 4.d0*mi2*s*t*usfcusf
     &  + 3.d0*mi2*s*t*sgbacsgba
     &  - 2.d0*mi2*s*t*tsfcusf2r
     &  - mi2*s*t*tsfcsgba2r
     &  + mi2*s*t*usfcsgba2r
     &  + 1.d0/2.d0*mi2*s*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*s*t**2*usfcsgbb2r
     &  + mi2*s**2*t*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*s**2*t*sgbacsgbb2r
     &  - 2.d0*mi2*s**2*usfcusf
     &  - 1.d0/2.d0*mi2*s**2*sgbacsgba
     &  - mi2*s**2*usfcsgba2r
     &  - 1.d0/4.d0*mi2*s**2*shbcsgbb2r
     &  - 1.d0/4.d0*mi2*s**2*sgbbcpoint2r
     &  + 1.d0/2.d0*mi2*s**3*usfcsgbb2r
     &  + 1.d0/4.d0*mi2*s**3*sgbacsgbb2r
     &  - 2.d0*mi2*t**2*tsfctsf
     &  - 2.d0*mi2*t**2*usfcusf
     &  - 2.d0*mi2*t**2*tsfcusf2r
     &  + 2.d0*mi2*t**2*tsfcsgba2r
     &  + 2.d0*mi2*t**2*usfcsgba2r
     &  - 1.d0/2.d0*mi2**2*mj2*m12*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*mj2*m12*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*mj2*m22*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*mj2*m22*usfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*mj2*s*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*mj2*s*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*m22*s*sgbbcsgbb
     &  - mi2**2*m12*m22*usfcsgbb2r
     &  - mi2**2*m12*m22*sgbacsgbb2r
     &  + 1.d0/4.d0*mi2**2*m12*m22**2*sgbbcsgbb
     &  + mi2**2*m12*s*usfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m12*s*sgbacsgbb2r
     &  + 1.d0/4.d0*mi2**2*m12*s**2*sgbbcsgbb
     &  + 1.d0/2.d0*mi2**2*m12*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m12*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12**2*m22*sgbbcsgbb
     &  - 1.d0/2.d0*mi2**2*m12**2*s*sgbbcsgbb
     &  - 1.d0/2.d0*mi2**2*m12**2*usfcsgbb2r
     &  - 1.d0/4.d0*mi2**2*m12**2*sgbacsgbb2r
     &  + 1.d0/4.d0*mi2**2*m12**3*sgbbcsgbb
     &  + mi2**2*m22*s*usfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m22*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2**2*m22*t*usfcsgbb2r
     &  - mi2**2*m22*sgbacsgba
     &  - 1.d0/2.d0*mi2**2*m22**2*usfcsgbb2r
     &  + 1.d0/4.d0*mi2**2*m22**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2**2*s*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*s*t*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*s**2*usfcsgbb2r
     &  - 1.d0/4.d0*mi2**2*s**2*sgbacsgbb2r
     &  - mj2*m12*m22*s*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*m22*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2*m12*m22*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*m22*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*m22*t*sgbacsgbb2r
     &  - 14.d0*mj2*m12*m22*usfcusf
     &  + mj2*m12*m22*tsfcusf2r
     &  + 1.d0/2.d0*mj2*m12*m22*tsfcsgba2r
     &  - 5.d0/2.d0*mj2*m12*m22*usfcsgba2r
     &  - 1.d0/2.d0*mj2*m12*m22*shbcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*m22*sgbbcpoint2r
     &  + 1.d0/4.d0*mj2*m12*m22**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2*m12*s*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12*s*t*usfcsgbb2r
     &  + 14.d0*mj2*m12*s*usfcusf
     &  + 4.d0*mj2*m12*s*sgbacsgba
     &  + mj2*m12*s*tsfcusf2r
     &  + 6.d0*mj2*m12*s*usfcsgba2r
     &  - 1.d0/2.d0*mj2*m12*s*shbcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*s*sgbbcpoint2r
     &  + mj2*m12*s**2*usfcsgbb2r
     &  + 1.d0/4.d0*mj2*m12*s**2*sgbacsgbb2r
     &  - 2.d0*mj2*m12*t*tsfctsf
     &  + 14.d0*mj2*m12*t*usfcusf
     &  + mj2*m12*t*sgbacsgba
     &  + 2.d0*mj2*m12*t*tsfcusf2r
     &  + 2.d0*mj2*m12*t*usfcsgba2r
     &  - 1.d0/2.d0*mj2*m12*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*t**2*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*tsfcshb2r
     &  - 1.d0/2.d0*mj2*m12*tsfcpoint2r
     &  - 3.d0/2.d0*mj2*m12*usfcshb2r
     &  - 3.d0/2.d0*mj2*m12*usfcpoint2r
     &  - 1.d0/2.d0*mj2*m12*shbcsgba2r
     &  - 1.d0/2.d0*mj2*m12*sgbacpoint2r
     &  - 1.d0/2.d0*mj2*m12**2*m22*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12**2*s*usfcsgbb2r
     &  + 1.d0/4.d0*mj2*m12**2*s*sgbacsgbb2r
     &  + 1.d0/2.d0*mj2*m12**2*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12**2*t*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12**2*t*sgbacsgbb2r
     &  - 6.d0*mj2*m12**2*usfcusf
     &  - 5.d0/2.d0*mj2*m12**2*sgbacsgba
     &  + 1.d0/2.d0*mj2*m12**2*tsfcsgba2r
     &  - 5.d0/2.d0*mj2*m12**2*usfcsgba2r
     &  + 1.d0/4.d0*mj2*m12**2*shbcsgbb2r
     &  + 1.d0/4.d0*mj2*m12**2*sgbbcpoint2r
     &  - 1.d0/4.d0*mj2*m12**3*sgbacsgbb2r
     &  + 2.d0*mj2*m22*s*t*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*s*t*sgbacsgbb2r
     &  + 16.d0*mj2*m22*s*usfcusf
     &  - mj2*m22*s*sgbacsgba
     &  - 1.d0/2.d0*mj2*m22*s*tsfcsgba2r
     &  + 5.d0/2.d0*mj2*m22*s*usfcsgba2r
     &  - 1.d0/2.d0*mj2*m22*s*shbcsgbb2r
     &  - 1.d0/2.d0*mj2*m22*s*sgbbcpoint2r
     &  + 3.d0/2.d0*mj2*m22*s**2*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*s**2*sgbacsgbb2r
     &  + 16.d0*mj2*m22*t*usfcusf
     &  - mj2*m22*t*sgbacsgba
     &  - mj2*m22*t*tsfcsgba2r
     &  - mj2*m22*t*usfcsgba2r
     &  + 1.d0/2.d0*mj2*m22*t**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*m22*tsfcshb2r
     &  + 1.d0/2.d0*mj2*m22*tsfcpoint2r
     &  - 5.d0/2.d0*mj2*m22*usfcshb2r
     &  - 5.d0/2.d0*mj2*m22*usfcpoint2r
     &  - 1.d0/2.d0*mj2*m22*shbcsgba2r
     &  - 1.d0/2.d0*mj2*m22*sgbacpoint2r
     &  - 3.d0/2.d0*mj2*m22**2*s*usfcsgbb2r
     &  - 1.d0/4.d0*mj2*m22**2*s*sgbacsgbb2r
     &  - mj2*m22**2*t*usfcsgbb2r
     &  - 8.d0*mj2*m22**2*usfcusf
     &  + 1.d0/2.d0*mj2*m22**2*sgbacsgba
     &  + mj2*m22**2*tsfcusf2r
     &  + 1.d0/4.d0*mj2*m22**2*shbcsgbb2r
     &  + 1.d0/4.d0*mj2*m22**2*sgbbcpoint2r
     &  + 1.d0/2.d0*mj2*m22**3*usfcsgbb2r
     &  - 16.d0*mj2*s*t*usfcusf
     &  + mj2*s*t*sgbacsgba
     &  - 4.d0*mj2*s*t*tsfcusf2r
     &  - 2.d0*mj2*s*t*usfcsgba2r
     &  - 1.d0/2.d0*mj2*s*t**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*s*t**2*usfcsgbb2r
     &  + 1.d0/2.d0*mj2*s*tsfcshb2r
     &  + 1.d0/2.d0*mj2*s*tsfcpoint2r
     &  + 3.d0/2.d0*mj2*s*usfcshb2r
     &  + 3.d0/2.d0*mj2*s*usfcpoint2r
     &  + 1.d0/2.d0*mj2*s*shbcsgba2r
     &  + 1.d0/2.d0*mj2*s*sgbacpoint2r
     &  - mj2*s**2*t*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*s**2*t*sgbacsgbb2r
     &  - 8.d0*mj2*s**2*usfcusf
     &  - 3.d0/2.d0*mj2*s**2*sgbacsgba
     &  - mj2*s**2*tsfcusf2r
     &  - 1.d0/2.d0*mj2*s**2*tsfcsgba2r
     &  - 7.d0/2.d0*mj2*s**2*usfcsgba2r
     &  + 1.d0/4.d0*mj2*s**2*shbcsgbb2r
     &  + 1.d0/4.d0*mj2*s**2*sgbbcpoint2r
     &  - 1.d0/2.d0*mj2*s**3*usfcsgbb2r
     &  - 1.d0/4.d0*mj2*s**3*sgbacsgbb2r
     &  + 2.d0*mj2*t*tsfcshb2r
     &  + 2.d0*mj2*t*tsfcpoint2r
     &  + 2.d0*mj2*t*usfcshb2r
     &  + 2.d0*mj2*t*usfcpoint2r
     &  + 4.d0*mj2*t**2*tsfctsf
     &  - 8.d0*mj2*t**2*usfcusf
     &  - 2.d0*mj2*t**2*tsfcusf2r
     &  + mj2*t**2*tsfcsgba2r
     &  + mj2*t**2*usfcsgba2r
     &  - 1.d0/2.d0*mj2**2*m12*m22*s*sgbbcsgbb
     &  + 1.d0/2.d0*mj2**2*m12*m22*tsfcsgbb2r
     &  - 3.d0/2.d0*mj2**2*m12*m22*usfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*m12*m22*sgbacsgbb2r
     &  + 1.d0/4.d0*mj2**2*m12*m22**2*sgbbcsgbb
     &  + 1.d0/2.d0*mj2**2*m12*s*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m12*s*usfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*m12*s*sgbacsgbb2r
     &  + 1.d0/4.d0*mj2**2*m12*s**2*sgbbcsgbb
     &  + 3.d0/2.d0*mj2**2*m12*t*tsfcsgbb2r
     &  + 3.d0/2.d0*mj2**2*m12*t*usfcsgbb2r
     &  + mj2**2*m12*tsfctsf
     &  - 9.d0*mj2**2*m12*usfcusf
     &  - 2.d0*mj2**2*m12*sgbacsgba
     &  - 2.d0*mj2**2*m12*tsfcusf2r
     &  - 3.d0/2.d0*mj2**2*m12*tsfcsgba2r
     &  - 7.d0/2.d0*mj2**2*m12*usfcsgba2r
     &  - 1.d0/2.d0*mj2**2*m12**2*m22*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*m12**2*s*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*m12**2*tsfcsgbb2r
     &  - mj2**2*m12**2*usfcsgbb2r
     &  - 3.d0/4.d0*mj2**2*m12**2*sgbacsgbb2r
     &  + 1.d0/4.d0*mj2**2*m12**3*sgbbcsgbb
     &  - mj2**2*m22*s*usfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*m22*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2**2*m22*t*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*m22*t*usfcsgbb2r
     &  - 10.d0*mj2**2*m22*usfcusf
     &  - mj2**2*m22*tsfcusf2r
     &  + 1.d0/2.d0*mj2**2*m22*tsfcsgba2r
     &  - 3.d0/2.d0*mj2**2*m22*usfcsgba2r
     &  + 1.d0/2.d0*mj2**2*m22**2*usfcsgbb2r
     &  + 1.d0/4.d0*mj2**2*m22**2*sgbacsgbb2r
     &  + 1.d0/2.d0*mj2**2*s*t*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2**2*s*t*usfcsgbb2r
     &  + 10.d0*mj2**2*s*usfcusf
     &  + mj2**2*s*sgbacsgba
     &  + 3.d0*mj2**2*s*tsfcusf2r
     &  + 3.d0/2.d0*mj2**2*s*tsfcsgba2r
     &  + 7.d0/2.d0*mj2**2*s*usfcsgba2r
     &  + 1.d0/2.d0*mj2**2*s**2*usfcsgbb2r
     &  + 1.d0/4.d0*mj2**2*s**2*sgbacsgbb2r
     &  - 2.d0*mj2**2*t*tsfctsf
     &  + 10.d0*mj2**2*t*usfcusf
     &  + 4.d0*mj2**2*t*tsfcusf2r
     &  + mj2**2*t*tsfcsgba2r
     &  + mj2**2*t*usfcsgba2r
     &  - mj2**2*tsfcshb2r
     &  - mj2**2*tsfcpoint2r
     &  - mj2**2*usfcshb2r
     &  - mj2**2*usfcpoint2r
     &  - mj2**3*m12*tsfcsgbb2r
     &  - mj2**3*m12*usfcsgbb2r
     &  - 4.d0*mj2**3*usfcusf
     &  - 2.d0*mj2**3*tsfcusf2r
     &  - mj2**3*tsfcsgba2r
     &  - mj2**3*usfcsgba2r
     &  + 10.d0*m12*m22*s*usfcusf
     &  - 3.d0/2.d0*m12*m22*s*sgbacsgba
     &  + 10.d0*m12*m22*t*usfcusf
     &  - m12*m22*t*sgbacsgba
     &  - m12*m22*t*tsfcusf2r
     &  - 1.d0/2.d0*m12*m22*t*tsfcsgba2r
     &  - 5.d0/2.d0*m12*m22*t*usfcsgba2r
     &  - 2.d0*m12*m22*usfcshb2r
     &  - 2.d0*m12*m22*usfcpoint2r
     &  - 5.d0*m12*m22**2*usfcusf
     &  + 1.d0/4.d0*m12*m22**2*sgbacsgba
     &  + m12*m22**2*usfcsgba2r
     &  - 10.d0*m12*s*t*usfcusf
     &  + 3.d0*m12*s*t*sgbacsgba
     &  - m12*s*t*tsfcusf2r
     &  + m12*s*usfcshb2r
     &  + m12*s*usfcpoint2r
     &  + 1.d0/2.d0*m12*s*shbcsgba2r
     &  + 1.d0/2.d0*m12*s*sgbacpoint2r
     &  - 5.d0*m12*s**2*usfcusf
     &  - 3.d0/4.d0*m12*s**2*sgbacsgba
     &  - 2.d0*m12*s**2*usfcsgba2r
     &  + 1.d0/2.d0*m12*t*tsfcshb2r
     &  + 1.d0/2.d0*m12*t*tsfcpoint2r
     &  + 3.d0/2.d0*m12*t*usfcshb2r
     &  + 3.d0/2.d0*m12*t*usfcpoint2r
     &  + 1.d0/2.d0*m12*t*shbcsgba2r
     &  + 1.d0/2.d0*m12*t*sgbacpoint2r
     &  + m12*t**2*tsfctsf
     &  - 5.d0*m12*t**2*usfcusf
     &  + m12*t**2*sgbacsgba
     &  + 3.d0/2.d0*m12*t**2*tsfcsgba2r
     &  + 3.d0/2.d0*m12*t**2*usfcsgba2r
     &  - 4.d0*m12**2*m22*usfcusf
     &  + 1.d0/2.d0*m12**2*m22*sgbacsgba
     &  + 1.d0/2.d0*m12**2*m22*usfcsgba2r
     &  + 4.d0*m12**2*s*usfcusf
     &  + m12**2*s*usfcsgba2r
     &  + 4.d0*m12**2*t*usfcusf
     &  - m12**2*t*sgbacsgba
     &  - 1.d0/2.d0*m12**2*t*tsfcsgba2r
     &  - 1.d0/2.d0*m12**2*t*usfcsgba2r
     &  - 1.d0/2.d0*m12**2*usfcshb2r
     &  - 1.d0/2.d0*m12**2*usfcpoint2r
     &  - 1.d0/4.d0*m12**2*shbcsgba2r
     &  - 1.d0/4.d0*m12**2*sgbacpoint2r
     &  - m12**3*usfcusf
     &  + 1.d0/4.d0*m12**3*sgbacsgba
     &  - 12.d0*m22*s*t*usfcusf
     &  + 2.d0*m22*s*t*sgbacsgba
     &  + 1.d0/2.d0*m22*s*t*tsfcsgba2r
     &  + 5.d0/2.d0*m22*s*t*usfcsgba2r
     &  + 2.d0*m22*s*usfcshb2r
     &  + 2.d0*m22*s*usfcpoint2r
     &  - 6.d0*m22*s**2*usfcusf
     &  + m22*s**2*sgbacsgba
     &  - 1.d0/2.d0*m22*s**2*usfcsgba2r
     &  - 1.d0/2.d0*m22*t*tsfcshb2r
     &  - 1.d0/2.d0*m22*t*tsfcpoint2r
     &  + 5.d0/2.d0*m22*t*usfcshb2r
     &  + 5.d0/2.d0*m22*t*usfcpoint2r
     &  - 1.d0/2.d0*m22*t*shbcsgba2r
     &  - 1.d0/2.d0*m22*t*sgbacpoint2r
     &  - 6.d0*m22*t**2*usfcusf
     &  + m22*t**2*tsfcusf2r
     &  + 1.d0/2.d0*m22*t**2*tsfcsgba2r
     &  + 5.d0/2.d0*m22*t**2*usfcsgba2r
     &  - m22*shbcshb
     &  - m22*pointcpoint
     &  - m22*shbcpoint2r
     &  + 6.d0*m22**2*s*usfcusf
     &  - 1.d0/2.d0*m22**2*s*sgbacsgba
     &  - m22**2*s*usfcsgba2r
     &  + 6.d0*m22**2*t*usfcusf
     &  - m22**2*t*tsfcusf2r
     &  - 2.d0*m22**2*t*usfcsgba2r
     &  - 3.d0/2.d0*m22**2*usfcshb2r
     &  - 3.d0/2.d0*m22**2*usfcpoint2r
     &  + 1.d0/4.d0*m22**2*shbcsgba2r
     &  + 1.d0/4.d0*m22**2*sgbacpoint2r
     &  - 2.d0*m22**3*usfcusf
     &  + 1.d0/2.d0*m22**3*usfcsgba2r
     &  - 1.d0/2.d0*s*t*tsfcshb2r
     &  - 1.d0/2.d0*s*t*tsfcpoint2r
     &  - 3.d0/2.d0*s*t*usfcshb2r
     &  - 3.d0/2.d0*s*t*usfcpoint2r
     &  - 1.d0/2.d0*s*t*shbcsgba2r
     &  - 1.d0/2.d0*s*t*sgbacpoint2r
     &  + 6.d0*s*t**2*usfcusf
     &  - 2.d0*s*t**2*sgbacsgba
     &  + s*t**2*tsfcusf2r
     &  - 3.d0/2.d0*s*t**2*tsfcsgba2r
     &  - 3.d0/2.d0*s*t**2*usfcsgba2r
     &  + 6.d0*s**2*t*usfcusf
     &  - 2.d0*s**2*t*sgbacsgba
     &  + s**2*t*tsfcusf2r
     &  + 1.d0/2.d0*s**2*t*tsfcsgba2r
     &  + 1.d0/2.d0*s**2*t*usfcsgba2r
     &  - 1.d0/2.d0*s**2*usfcshb2r
     &  - 1.d0/2.d0*s**2*usfcpoint2r
     &  - 1.d0/4.d0*s**2*shbcsgba2r
     &  - 1.d0/4.d0*s**2*sgbacpoint2r
     &  + 2.d0*s**3*usfcusf
     &  + 1.d0/2.d0*s**3*sgbacsgba
     &  + s**3*usfcsgba2r
     &  - t**2*tsfcshb2r
     &  - t**2*tsfcpoint2r
     &  - t**2*usfcshb2r
     &  - t**2*usfcpoint2r
     &  - 2.d0*t**3*tsfctsf
     &  + 2.d0*t**3*usfcusf
     &  - t**3*tsfcsgba2r
     &  - t**3*usfcsgba2r
     &  )
      cmog0 =
     &  + 10.d0*mi2*mj2*m12*m22*sgbbcsgbb
     &  - 2.d0*mi2*mj2*m12*s*sgbbcsgbb
     &  - 4.d0*mi2*mj2*m12*usfcsgbb2r
     &  - 2.d0*mi2*mj2*m12*sgbacsgbb2r
     &  - 5.d0*mi2*mj2*m12**2*sgbbcsgbb
     &  - 2.d0*mi2*mj2*m22*s*sgbbcsgbb
     &  - 4.d0*mi2*mj2*m22*usfcsgbb2r
     &  - 2.d0*mi2*mj2*m22*sgbacsgbb2r
     &  - 5.d0*mi2*mj2*m22**2*sgbbcsgbb
     &  + mi2*mj2*s*tsfcsgbb2r
     &  + 3.d0*mi2*mj2*s*usfcsgbb2r
     &  + mi2*mj2*s*sgbacsgbb2r
     &  + mi2*mj2*s**2*sgbbcsgbb
     &  + 4.d0*mi2*mj2*t*tsfcsgbb2r
     &  + 4.d0*mi2*mj2*t*usfcsgbb2r
     &  + 4.d0*mi2*mj2*tsfctsf
     &  + 20.d0*mi2*mj2*usfcusf
     &  + 14.d0*mi2*mj2*sgbacsgba
     &  + 8.d0*mi2*mj2*tsfcusf2r
     &  - 2.d0*mi2*mj2*tsfcsgba2r
     &  + 2.d0*mi2*mj2*usfcsgba2r
     &  - mi2*mj2**2*tsfcsgbb2r
     &  - mi2*mj2**2*usfcsgbb2r
     &  + 1.d0/2.d0*mi2*m12*s*tsfcsgbb2r
     &  + 3.d0/2.d0*mi2*m12*s*usfcsgbb2r
     &  - 5.d0/2.d0*mi2*m12*s*sgbacsgbb2r
     &  - mi2*m12*t*tsfcsgbb2r
     &  + mi2*m12*t*usfcsgbb2r
     &  - 5.d0*mi2*m12*t*sgbacsgbb2r
     &  + 8.d0*mi2*m12*usfcusf
     &  + 7.d0*mi2*m12*sgbacsgba
     &  - 2.d0*mi2*m12*tsfcusf2r
     &  - 2.d0*mi2*m12*tsfcsgba2r
     &  - 3.d0*mi2*m12*shbcsgbb2r
     &  - 3.d0*mi2*m12*sgbbcpoint2r
     &  + 1.d0/2.d0*mi2*m12**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2*m12**2*usfcsgbb2r
     &  + 5.d0/2.d0*mi2*m12**2*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2*m22*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mi2*m22*s*usfcsgbb2r
     &  + 5.d0/2.d0*mi2*m22*s*sgbacsgbb2r
     &  + mi2*m22*t*tsfcsgbb2r
     &  - mi2*m22*t*usfcsgbb2r
     &  + 5.d0*mi2*m22*t*sgbacsgbb2r
     &  - 2.d0*mi2*m22*tsfctsf
     &  + 10.d0*mi2*m22*usfcusf
     &  + 7.d0*mi2*m22*sgbacsgba
     &  - 2.d0*mi2*m22*tsfcusf2r
     &  + 2.d0*mi2*m22*usfcsgba2r
     &  + 3.d0*mi2*m22*shbcsgbb2r
     &  + 3.d0*mi2*m22*sgbbcpoint2r
     &  - 1.d0/2.d0*mi2*m22**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mi2*m22**2*usfcsgbb2r
     &  - 5.d0/2.d0*mi2*m22**2*sgbacsgbb2r
     &  - 12.d0*mi2*s*usfcusf
     &  + 3.d0*mi2*s*sgbacsgba
     &  - 4.d0*mi2*s*tsfcusf2r
     &  - 2.d0*mi2*s*tsfcsgba2r
     &  + 4.d0*mi2*s*usfcsgba2r
     &  + 4.d0*mi2*t*tsfctsf
     &  - 12.d0*mi2*t*usfcusf
     &  - 10.d0*mi2*t*sgbacsgba
     &  - 2.d0*mi2*t*tsfcsgba2r
     &  + 2.d0*mi2*t*usfcsgba2r
     &  + mi2*tsfcshb2r
     &  + mi2*tsfcpoint2r
     &  + 3.d0*mi2*usfcshb2r
     &  + 3.d0*mi2*usfcpoint2r
     &  - 3.d0*mi2*shbcsgba2r
     &  - 3.d0*mi2*sgbacpoint2r
     &  - mi2**2*mj2*tsfcsgbb2r
     &  - mi2**2*mj2*usfcsgbb2r
     &  - 5.d0*mi2**2*m12*m22*sgbbcsgbb
     &  + mi2**2*m12*s*sgbbcsgbb
     &  - 3.d0/2.d0*mi2**2*m12*tsfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*m12*usfcsgbb2r
     &  + 7.d0/2.d0*mi2**2*m12*sgbacsgbb2r
     &  + 5.d0/2.d0*mi2**2*m12**2*sgbbcsgbb
     &  + mi2**2*m22*s*sgbbcsgbb
     &  + 3.d0/2.d0*mi2**2*m22*tsfcsgbb2r
     &  + 9.d0/2.d0*mi2**2*m22*usfcsgbb2r
     &  - 3.d0/2.d0*mi2**2*m22*sgbacsgbb2r
     &  + 5.d0/2.d0*mi2**2*m22**2*sgbbcsgbb
     &  - 1.d0/2.d0*mi2**2*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mi2**2*s*usfcsgbb2r
     &  - 1.d0/2.d0*mi2**2*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mi2**2*s**2*sgbbcsgbb
     &  - 2.d0*mi2**2*t*tsfcsgbb2r
     &  - 2.d0*mi2**2*t*usfcsgbb2r
     &  + 8.d0*mi2**2*usfcusf
     &  - 2.d0*mi2**2*sgbacsgba
     &  + 4.d0*mi2**2*tsfcusf2r
     &  + 2.d0*mi2**2*tsfcsgba2r
     &  - 2.d0*mi2**2*usfcsgba2r
     &  + mi2**3*tsfcsgbb2r
     &  + mi2**3*usfcsgbb2r
     &  - 1.d0/2.d0*mj2*m12*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mj2*m12*s*usfcsgbb2r
     &  + 5.d0/2.d0*mj2*m12*s*sgbacsgbb2r
     &  + mj2*m12*t*tsfcsgbb2r
     &  - mj2*m12*t*usfcsgbb2r
     &  + 5.d0*mj2*m12*t*sgbacsgbb2r
     &  - 2.d0*mj2*m12*tsfctsf
     &  + 10.d0*mj2*m12*usfcusf
     &  + 7.d0*mj2*m12*sgbacsgba
     &  - 2.d0*mj2*m12*tsfcusf2r
     &  + 2.d0*mj2*m12*usfcsgba2r
     &  + 3.d0*mj2*m12*shbcsgbb2r
     &  + 3.d0*mj2*m12*sgbbcpoint2r
     &  - 1.d0/2.d0*mj2*m12**2*tsfcsgbb2r
     &  + 1.d0/2.d0*mj2*m12**2*usfcsgbb2r
     &  - 5.d0/2.d0*mj2*m12**2*sgbacsgbb2r
     &  + 1.d0/2.d0*mj2*m22*s*tsfcsgbb2r
     &  + 3.d0/2.d0*mj2*m22*s*usfcsgbb2r
     &  - 5.d0/2.d0*mj2*m22*s*sgbacsgbb2r
     &  - mj2*m22*t*tsfcsgbb2r
     &  + mj2*m22*t*usfcsgbb2r
     &  - 5.d0*mj2*m22*t*sgbacsgbb2r
     &  + 8.d0*mj2*m22*usfcusf
     &  + 7.d0*mj2*m22*sgbacsgba
     &  - 2.d0*mj2*m22*tsfcusf2r
     &  - 2.d0*mj2*m22*tsfcsgba2r
     &  - 3.d0*mj2*m22*shbcsgbb2r
     &  - 3.d0*mj2*m22*sgbbcpoint2r
     &  + 1.d0/2.d0*mj2*m22**2*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2*m22**2*usfcsgbb2r
     &  + 5.d0/2.d0*mj2*m22**2*sgbacsgbb2r
     &  - 12.d0*mj2*s*usfcusf
     &  + 3.d0*mj2*s*sgbacsgba
     &  - 4.d0*mj2*s*tsfcusf2r
     &  - 2.d0*mj2*s*tsfcsgba2r
     &  + 4.d0*mj2*s*usfcsgba2r
     &  + 4.d0*mj2*t*tsfctsf
     &  - 12.d0*mj2*t*usfcusf
     &  - 10.d0*mj2*t*sgbacsgba
     &  - 2.d0*mj2*t*tsfcsgba2r
     &  + 2.d0*mj2*t*usfcsgba2r
     &  + mj2*tsfcshb2r
     &  + mj2*tsfcpoint2r
     &  + 3.d0*mj2*usfcshb2r
     &  + 3.d0*mj2*usfcpoint2r
     &  - 3.d0*mj2*shbcsgba2r
     &  - 3.d0*mj2*sgbacpoint2r
     &  - 5.d0*mj2**2*m12*m22*sgbbcsgbb
     &  + mj2**2*m12*s*sgbbcsgbb
     &  + 3.d0/2.d0*mj2**2*m12*tsfcsgbb2r
     &  + 9.d0/2.d0*mj2**2*m12*usfcsgbb2r
     &  - 3.d0/2.d0*mj2**2*m12*sgbacsgbb2r
     &  + 5.d0/2.d0*mj2**2*m12**2*sgbbcsgbb
     &  + mj2**2*m22*s*sgbbcsgbb
     &  - 3.d0/2.d0*mj2**2*m22*tsfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*m22*usfcsgbb2r
     &  + 7.d0/2.d0*mj2**2*m22*sgbacsgbb2r
     &  + 5.d0/2.d0*mj2**2*m22**2*sgbbcsgbb
     &  - 1.d0/2.d0*mj2**2*s*tsfcsgbb2r
     &  - 3.d0/2.d0*mj2**2*s*usfcsgbb2r
     &  - 1.d0/2.d0*mj2**2*s*sgbacsgbb2r
     &  - 1.d0/2.d0*mj2**2*s**2*sgbbcsgbb
     &  - 2.d0*mj2**2*t*tsfcsgbb2r
     &  - 2.d0*mj2**2*t*usfcsgbb2r
     &  + 8.d0*mj2**2*usfcusf
     &  - 2.d0*mj2**2*sgbacsgba
     &  + 4.d0*mj2**2*tsfcusf2r
     &  + 2.d0*mj2**2*tsfcsgba2r
     &  - 2.d0*mj2**2*usfcsgba2r
     &  + mj2**3*tsfcsgbb2r
     &  + mj2**3*usfcsgbb2r
     &  + m12*m22*tsfctsf
     &  + 5.d0*m12*m22*usfcusf
     &  + 5.d0*m12*m22*sgbacsgba
     &  + m12*m22*tsfcusf2r
     &  + m12*m22*tsfcsgba2r
     &  - m12*m22*usfcsgba2r
     &  - 6.d0*m12*s*usfcusf
     &  - 6.d0*m12*s*sgbacsgba
     &  + m12*s*tsfcusf2r
     &  + m12*s*tsfcsgba2r
     &  + m12*s*usfcsgba2r
     &  - 2.d0*m12*t*tsfctsf
     &  - 6.d0*m12*t*usfcusf
     &  - 10.d0*m12*t*sgbacsgba
     &  - 2.d0*m12*t*tsfcsgba2r
     &  + 2.d0*m12*t*usfcsgba2r
     &  - 1.d0/2.d0*m12*tsfcshb2r
     &  - 1.d0/2.d0*m12*tsfcpoint2r
     &  + 3.d0/2.d0*m12*usfcshb2r
     &  + 3.d0/2.d0*m12*usfcpoint2r
     &  - 3.d0*m12*shbcsgba2r
     &  - 3.d0*m12*sgbacpoint2r
     &  + 2.d0*m12**2*usfcusf
     &  + 5.d0/2.d0*m12**2*sgbacsgba
     &  + 1.d0/2.d0*m12**2*tsfcsgba2r
     &  - 1.d0/2.d0*m12**2*usfcsgba2r
     &  - 6.d0*m22*s*usfcusf
     &  - 6.d0*m22*s*sgbacsgba
     &  + m22*s*tsfcusf2r
     &  + m22*s*tsfcsgba2r
     &  + m22*s*usfcsgba2r
     &  - 2.d0*m22*t*tsfctsf
     &  - 6.d0*m22*t*usfcusf
     &  - 10.d0*m22*t*sgbacsgba
     &  - 2.d0*m22*t*tsfcsgba2r
     &  + 2.d0*m22*t*usfcsgba2r
     &  - 1.d0/2.d0*m22*tsfcshb2r
     &  - 1.d0/2.d0*m22*tsfcpoint2r
     &  + 3.d0/2.d0*m22*usfcshb2r
     &  + 3.d0/2.d0*m22*usfcpoint2r
     &  - 3.d0*m22*shbcsgba2r
     &  - 3.d0*m22*sgbacpoint2r
     &  + 2.d0*m22**2*usfcusf
     &  + 5.d0/2.d0*m22**2*sgbacsgba
     &  + 1.d0/2.d0*m22**2*tsfcsgba2r
     &  - 1.d0/2.d0*m22**2*usfcsgba2r
     &  + 8.d0*s*t*usfcusf
     &  + 10.d0*s*t*sgbacsgba
     &  - 4.d0*s*t*usfcsgba2r
     &  - 1.d0/2.d0*s*tsfcshb2r
     &  - 1.d0/2.d0*s*tsfcpoint2r
     &  - 5.d0/2.d0*s*usfcshb2r
     &  - 5.d0/2.d0*s*usfcpoint2r
     &  + 3.d0*s*shbcsgba2r
     &  + 3.d0*s*sgbacpoint2r
     &  + 4.d0*s**2*usfcusf
     &  - 3.d0/2.d0*s**2*sgbacsgba
     &  + s**2*tsfcusf2r
     &  + 1.d0/2.d0*s**2*tsfcsgba2r
     &  - 5.d0/2.d0*s**2*usfcsgba2r
     &  + 2.d0*t*tsfcshb2r
     &  + 2.d0*t*tsfcpoint2r
     &  - 2.d0*t*usfcshb2r
     &  - 2.d0*t*usfcpoint2r
     &  + 6.d0*t*shbcsgba2r
     &  + 6.d0*t*sgbacpoint2r
     &  + 4.d0*t**2*tsfctsf
     &  + 4.d0*t**2*usfcusf
     &  + 10.d0*t**2*sgbacsgba
     &  + 2.d0*t**2*tsfcsgba2r
     &  - 2.d0*t**2*usfcsgba2r
     &  + 4.d0*shbcshb
     &  + 4.d0*pointcpoint
     &  + 4.d0*shbcpoint2r

      dsasgbgbas = cmog12+cmog1+cmog2+cmog0


      if (dsasgbgbas.lt.0.0d0) then
        if(aszeroprint) then
          write(*,*) ' '
          write(*,*) 'DS: negative cross section in dsasgbgbas:'
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
          write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
          write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
          write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
          write(*,*) 'DS: dsasgbgbas = ',dsasgbgbas
          write(*,*) 'DS: mass1=',mass1
          write(*,*) 'DS: mass2=',mass2
          write(*,*) 'DS: mass3=',mass3
          write(*,*) 'DS: mass4=',mass4
          write(*,*) 'DS: s=',s
          write(*,*) 'DS: t=',t
          write(*,*) 'DS: u=',u
          write(*,*) 'DS: p12=',p12
          write(*,*) 'DS: costheta=',costheta
          write(*,*) 'DS: c factors : ',csum,
     &                    c1,c2,c3,c4,c6,c8,c11,c12,c13,c15
          write(*,*) 'DS: c ratios : ',c1/csum,c2/csum,c3/csum,c4/csum,
     &                    c6/csum,c8/csum,c11/csum,c12/csum,c13/csum,
     &                    c15/csum
          write(*,*) 'DS: partial res : ',cmog12,cmog1,cmog2,cmog0
          write(*,*) 'DS: partial res ratios : ',cmog12/dsasgbgbas,
     &                    cmog1/dsasgbgbas,cmog2/dsasgbgbas,
     &                    cmog0/dsasgbgbas
        endif
        dsasgbgbas=0.0d0
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsasgbgbas)*k34/(8.0d0*pi*gg1*gg2*s34*dsqrt(s))

      return
      end

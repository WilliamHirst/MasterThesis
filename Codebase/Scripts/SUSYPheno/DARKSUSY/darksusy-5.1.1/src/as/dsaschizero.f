**************************************************************
*** SUBROUTINE dsaschizero                                 ***
*** computes dW_{ij}/dcostheta                             ***
*** sfermion(i) + neutralino(j) -> gamma/gluon + fermion   ***
*** ampl2 obtained with sum over physical polarizations    ***
***                                                        *** 
*** input askin variables: p12,costheta                    ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 02-06-13                                         ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      SUBROUTINE dsaschizero(kp1,kp2,kp3,kp4,icase,result)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase 
      integer kgb,kfer,kfers,ksfert,ksfer1,ksfer2,kfercheck,kchiu
      real*8 result,colfact
      double complex dsasdepro,ampl2,den
      double complex A0L,A0R,A1L,A1R,A2L1,A2R1,A2L2,A2R2,
     & A2L3,A2R3,A0Ls,A0Rs,A1Ls,A1Rs,A2L1s,A2R1s,A2L2s,A2R2s,
     & A2L3s,A2R3s
      real*8 e1p1,e1p2,k1p1,k1p2,p1p2,mass2b,mass4b
      integer i,j
      real*8 pi
      parameter (pi=3.141592653589793238d0)
*****
c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        result=0.d0
        return
      endif      

      kgb=kp3
      kfer=kp4

*****
      colfact=1.d0
      if(iifam(1).ge.ivfam(ksu1).and.iifam(1).le.ivfam(ksb2))
     &  colfact=3.d0 
*****
***** non equal particle final state: 
      s34=1.d0
***** masses in final state:
      mass3=mass(kgb)
      mass4=mass(kfer)
***** define the kinematic variables
      call dsaskinset2
*****
      p1p2=(Svar-(mass1**2+mass2**2))/2.d0
      k1p1=ep1*ek3-p12*k34*costheta
      k1p2=ep2*ek3+p12*k34*costheta
      mass2b=mass2
      mass4b=mass4
      ampl2=(0.d0,0.d0)
*****
*****
***** sum over two polarizations
      do j=1,2
      if(j.eq.1) then   
***** first polarization:
        e1p1=p12*dsqrt(1.d0-costheta**2)
        e1p2=-p12*dsqrt(1.d0-costheta**2)
***** second polarization:
      elseif(j.eq.2) then  
        e1p1=0.d0
        e1p2=0.d0
      else
        write(*,*) 'wrong polarization index'
        stop
      endif  
*****
***** initialization  
      A0L=(0.d0,0.d0)
      A0R=(0.d0,0.d0)
      A1L=(0.d0,0.d0)
      A1R=(0.d0,0.d0)
      A2L1=(0.d0,0.d0)
      A2R1=(0.d0,0.d0)
      A2L2=(0.d0,0.d0)
      A2R2=(0.d0,0.d0)
      A2L3=(0.d0,0.d0)
      A2R3=(0.d0,0.d0)
*****
*****
***** 
***** the seventh case   
***** sfermion(i) + neutralino(j) -> photon + fermion
***** (not allowed for neutrinos)
*****
      if(icase.eq.7) then
*****
***** f in s-channel ******************************************
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R1=A2R1-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A2L2=A2L2-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R2=A2R2-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A1L=A1L-den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
        A1R=A1R-den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L-e1p1*den
     &      *2.d0*gl(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
          A0R=A0R-e1p1*den
     &      *2.d0*gr(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
        enddo
        goto 100
      endif  
*****
*****
*****
***** the nineth case
***** squark(i) + neutralino(j) 
***** -> gluon + quark
      if(icase.eq.9) then
        colfact=4.d0
*****
***** f in s-channel ******************************************
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R1=A2R1-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A2L2=A2L2-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R2=A2R2-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A1L=A1L-den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
        A1R=A1R-den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L-e1p1*den
     &      *2.d0*gl(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
          A0R=A0R-e1p1*den
     &      *2.d0*gr(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
        enddo
        goto 100
      endif  
*****
*****
***** 
***** the eighth case    
***** down-type-sfermion(i) + chargino^+(j) 
***** -> photon + up-type-fermion
      if(icase.eq.8) then
*****
***** f in s-channel ******************************************
        if(dabs(echarg(kfer)).gt.1.d-16) then
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R1=A2R1-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A2L2=A2L2-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R2=A2R2-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A1L=A1L-den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
        A1R=A1R-den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
        endif 
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L-e1p1*den
     &      *2.d0*gl(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
          A0R=A0R-e1p1*den
     &      *2.d0*gr(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
        enddo
*****
***** kcha(1) in u-channel ************************************
        kchiu=kcha1
        den=dsasdepro(Uvar,kchiu)
        A0L=A0L-2.d0*e1p2*den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A0R=A0R-2.d0*e1p2*den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
        A1L=A1L-den*gr(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A1R=A1R-den*gl(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A2L2=A2L2+den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A2R2=A2R2+den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
        A2L3=A2L3-den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A2R3=A2R3-den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
***** kcha(2) in u-channel ************************************
        kchiu=kcha2
        den=dsasdepro(Uvar,kchiu)
        A0L=A0L-2.d0*e1p2*den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A0R=A0R-2.d0*e1p2*den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
        A1L=A1L-den*gr(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A1R=A1R-den*gl(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A2L2=A2L2+den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A2R2=A2R2+den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
        A2L3=A2L3-den*gl(kp1,kfer,kchiu)*gl(kgb,kchiu,kp2)
        A2R3=A2R3-den*gr(kp1,kfer,kchiu)*gr(kgb,kchiu,kp2)
*****
        goto 100
      endif
*****
*****
*****
***** the tenth case    
***** down-type-squark(i) + chargino^+(j) 
***** -> gluon + up-type-quark
      if(icase.eq.10) then
***** set the color factors
        colfact=4.d0
*****
***** f in s-channel ******************************************
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R1=A2R1-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A2L2=A2L2-den
     &         *gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2)
        A2R2=A2R2-den
     &         *gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2)
        A1L=A1L-den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
        A1R=A1R-den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L-e1p1*den
     &      *2.d0*gl(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1) 
          A0R=A0R-e1p1*den
     &      *2.d0*gr(ksfert,kfer,kp2)*gl(kgb,ksfert,kp1)
        enddo
        goto 100
      endif    
*****
*****
*****
***** the fifth case 
***** anti-up-type-sfermion(i) + \chi^+(j) 
***** -> photon + anti-down-type-fermion
      if(icase.eq.5) then
*****
***** f in s-channel ******************************************
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1+den
     &         *gl(kp2,kfers,kp1)*gl(kgb,kfers,kfer)
        A2R1=A2R1+den
     &         *gr(kp2,kfers,kp1)*gr(kgb,kfers,kfer)
        A2L2=A2L2+den
     &         *gl(kp2,kfers,kp1)*gl(kgb,kfers,kfer)
        A2R2=A2R2+den
     &         *gr(kp2,kfers,kp1)*gr(kgb,kfers,kfer)
        A1L=A1L+den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
        A1R=A1R+den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
*****
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L+e1p1*den
     &      *2.d0*gl(kp2,kfer,ksfert)*gl(kgb,kp1,ksfert) 
          A0R=A0R+e1p1*den
     &      *2.d0*gr(kp2,kfer,ksfert)*gl(kgb,kp1,ksfert) 
        enddo
*****
***** kcha(1) in u-channel ************************************
        kchiu=kcha1
        den=dsasdepro(Uvar,kchiu)
        A0L=A0L-2.d0*e1p2*den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A0R=A0R-2.d0*e1p2*den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
        A1L=A1L-den*gr(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A1R=A1R-den*gl(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A2L2=A2L2+den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A2R2=A2R2+den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
        A2L3=A2L3-den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A2R3=A2R3-den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
***** kcha(2) in u-channel ************************************
        kchiu=kcha2
        den=dsasdepro(Uvar,kchiu)
        A0L=A0L-2.d0*e1p2*den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A0R=A0R-2.d0*e1p2*den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
        A1L=A1L-den*gr(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A1R=A1R-den*gl(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
     &       *mass(kchiu)
        A2L2=A2L2+den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A2R2=A2R2+den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
        A2L3=A2L3-den*gl(kchiu,kfer,kp1)*gl(kgb,kchiu,kp2)
        A2R3=A2R3-den*gr(kchiu,kfer,kp1)*gr(kgb,kchiu,kp2)
        goto 100
      endif
*****
*****
*****
***** the sixth case 
***** anti-up-type-squark(i) + \chi^+(j) 
***** -> gluon + anti-down-type-quark
      if(icase.eq.6) then
***** set the color factors
        colfact=4.d0
*****
***** f in s-channel ******************************************
        kfers=kfer
        den=dsasdepro(Svar,kfers)
        A2L1=A2L1+den
     &         *gl(kp2,kfers,kp1)*gl(kgb,kfers,kfer)
        A2R1=A2R1+den
     &         *gr(kp2,kfers,kp1)*gr(kgb,kfers,kfer)
        A2L2=A2L2+den
     &         *gl(kp2,kfers,kp1)*gl(kgb,kfers,kfer)
        A2R2=A2R2+den
     &         *gr(kp2,kfers,kp1)*gr(kgb,kfers,kfer)
        A1L=A1L+den
     &      *(mass(kfers)*gl(kgb,kfer,kfers)*gr(kp1,kfers,kp2))
        A1R=A1R+den
     &      *(mass(kfers)*gr(kgb,kfer,kfers)*gl(kp1,kfers,kp2))
*****
***** ksfer1 and ksfer2 in t-channel **************************
        do i=1,ncsfert
          ksfert=kcsfertn(i)
          den=dsasdepro(Tvar,ksfert)
          A0L=A0L+e1p1*den
     &      *2.d0*gl(kp2,kfer,ksfert)*gl(kgb,kp1,ksfert) 
          A0R=A0R+e1p1*den
     &      *2.d0*gr(kp2,kfer,ksfert)*gl(kgb,kp1,ksfert) 
        enddo
        goto 100
      endif
*****
*****
*****


 100  continue

      A0Ls=dconjg(A0L)
      A0Rs=dconjg(A0R)
      A1Ls=dconjg(A1L)
      A1Rs=dconjg(A1R)
      A2L1s=dconjg(A2L1)
      A2L2s=dconjg(A2L2)
      A2L3s=dconjg(A2L3)
      A2R1s=dconjg(A2R1)
      A2R2s=dconjg(A2R2)
      A2R3s=dconjg(A2R3)



*****
***** amplitude squared *************************************** 
*****

      ampl2=ampl2
     &  - 2.d0*mass2b*mass4b*mass1**2*A2L1*A2R1s
     &  - 2.d0*mass2b*mass4b*mass1**2*A2L1s*A2R1
     &  - 2.d0*mass2b*mass4b*A1L*A1Rs
     &  - 2.d0*mass2b*mass4b*A1R*A1Ls
     &  + 2.d0*mass2b*mass4b*A0L*A0Rs
     &  + 2.d0*mass2b*mass4b*A0L*e1p1*A2R1s
     &  - 2.d0*mass2b*mass4b*A0L*e1p1*A2R2s
     &  + 2.d0*mass2b*mass4b*A0R*A0Ls
     &  + 2.d0*mass2b*mass4b*A0R*e1p1*A2L1s
     &  - 2.d0*mass2b*mass4b*A0R*e1p1*A2L2s
     &  + 2.d0*mass2b*mass4b*A0Ls*e1p1*A2R1
     &  - 2.d0*mass2b*mass4b*A0Ls*e1p1*A2R2
     &  + 2.d0*mass2b*mass4b*A0Rs*e1p1*A2L1
     &  - 2.d0*mass2b*mass4b*A0Rs*e1p1*A2L2
     &  - 2.d0*mass2b*mass4b*p1p2*A2L1*A2R2s
     &  - 2.d0*mass2b*mass4b*p1p2*A2L1s*A2R2
     &  - 2.d0*mass2b*mass4b*p1p2*A2L2*A2R1s
     &  - 2.d0*mass2b*mass4b*p1p2*A2L2s*A2R1
     &  - 2.d0*mass2b*mass4b*k1p1*A2L1*A2R3s
     &  - 2.d0*mass2b*mass4b*k1p1*A2L1s*A2R3
     &  - 2.d0*mass2b*mass4b*k1p1*A2L3*A2R1s
     &  - 2.d0*mass2b*mass4b*k1p1*A2L3s*A2R1
     &  - 2.d0*mass2b*mass4b*k1p2*A2L2*A2R3s
     &  - 2.d0*mass2b*mass4b*k1p2*A2L2s*A2R3
     &  - 2.d0*mass2b*mass4b*k1p2*A2L3*A2R2s
     &  - 2.d0*mass2b*mass4b*k1p2*A2L3s*A2R2
     &  + 2.d0*mass2b*mass1**2*A1L*A2R1s
     &  + 2.d0*mass2b*mass1**2*A1R*A2L1s
     &  + 2.d0*mass2b*mass1**2*A1Ls*A2R1
     &  + 2.d0*mass2b*mass1**2*A1Rs*A2L1
     &  + 2.d0*mass2b*A1L*p1p2*A2R1s
     &  + 2.d0*mass2b*A1L*p1p2*A2R2s
     &  - 2.d0*mass2b*A1L*k1p1*A2R1s
     &  + 2.d0*mass2b*A1L*k1p1*A2R3s
     &  - 2.d0*mass2b*A1L*k1p2*A2R2s
     &  + 2.d0*mass2b*A1L*k1p2*A2R3s
     &  + 2.d0*mass2b*A1R*p1p2*A2L1s
     &  + 2.d0*mass2b*A1R*p1p2*A2L2s
     &  - 2.d0*mass2b*A1R*k1p1*A2L1s
     &  + 2.d0*mass2b*A1R*k1p1*A2L3s
     &  - 2.d0*mass2b*A1R*k1p2*A2L2s
     &  + 2.d0*mass2b*A1R*k1p2*A2L3s
     &  + 2.d0*mass2b*A1Ls*p1p2*A2R1
     &  + 2.d0*mass2b*A1Ls*p1p2*A2R2
     &  - 2.d0*mass2b*A1Ls*k1p1*A2R1
     &  + 2.d0*mass2b*A1Ls*k1p1*A2R3
     &  - 2.d0*mass2b*A1Ls*k1p2*A2R2
     &  + 2.d0*mass2b*A1Ls*k1p2*A2R3
     &  + 2.d0*mass2b*A1Rs*p1p2*A2L1
     &  + 2.d0*mass2b*A1Rs*p1p2*A2L2
     &  - 2.d0*mass2b*A1Rs*k1p1*A2L1
     &  + 2.d0*mass2b*A1Rs*k1p1*A2L3
     &  - 2.d0*mass2b*A1Rs*k1p2*A2L2
     &  + 2.d0*mass2b*A1Rs*k1p2*A2L3
     &  - 2.d0*mass2b**2*mass4b*A1L*A2L2s
     &  - 2.d0*mass2b**2*mass4b*A1R*A2R2s
     &  - 2.d0*mass2b**2*mass4b*A1Ls*A2L2
     &  - 2.d0*mass2b**2*mass4b*A1Rs*A2R2
     &  - 2.d0*mass2b**2*mass1**2*A2L1*A2L1s
     &  + 2.d0*mass2b**2*mass1**2*A2L1*A2L2s
     &  + 2.d0*mass2b**2*mass1**2*A2L1s*A2L2
     &  - 2.d0*mass2b**2*mass1**2*A2R1*A2R1s
     &  + 2.d0*mass2b**2*mass1**2*A2R1*A2R2s
     &  + 2.d0*mass2b**2*mass1**2*A2R1s*A2R2
     &  + 2.d0*mass2b**2*A1L*A1Ls
     &  + 2.d0*mass2b**2*A1R*A1Rs
     &  + 2.d0*mass2b**2*A0L*A0Ls
     &  + 2.d0*mass2b**2*A0L*e1p1*A2L1s
     &  + 2.d0*mass2b**2*A0R*A0Rs
     &  + 2.d0*mass2b**2*A0R*e1p1*A2R1s
     &  + 2.d0*mass2b**2*A0Ls*e1p1*A2L1
     &  + 2.d0*mass2b**2*A0Rs*e1p1*A2R1
     &  + 2.d0*mass2b**2*p1p2*A2L1*A2L2s
     &  + 2.d0*mass2b**2*p1p2*A2L1s*A2L2
     &  + 2.d0*mass2b**2*p1p2*A2L2*A2L2s
     &  + 2.d0*mass2b**2*p1p2*A2R1*A2R2s
     &  + 2.d0*mass2b**2*p1p2*A2R1s*A2R2
     &  + 2.d0*mass2b**2*p1p2*A2R2*A2R2s
     &  - 2.d0*mass2b**2*k1p1*A2L1*A2L2s
     &  - 2.d0*mass2b**2*k1p1*A2L1*A2L3s
     &  - 2.d0*mass2b**2*k1p1*A2L1s*A2L2
     &  - 2.d0*mass2b**2*k1p1*A2L1s*A2L3
     &  + 2.d0*mass2b**2*k1p1*A2L2*A2L3s
     &  + 2.d0*mass2b**2*k1p1*A2L2s*A2L3
     &  - 2.d0*mass2b**2*k1p1*A2R1*A2R2s
     &  - 2.d0*mass2b**2*k1p1*A2R1*A2R3s
     &  - 2.d0*mass2b**2*k1p1*A2R1s*A2R2
     &  - 2.d0*mass2b**2*k1p1*A2R1s*A2R3
     &  + 2.d0*mass2b**2*k1p1*A2R2*A2R3s
     &  + 2.d0*mass2b**2*k1p1*A2R2s*A2R3
     &  - 2.d0*mass2b**2*k1p2*A2L2*A2L2s
     &  + 2.d0*mass2b**2*k1p2*A2L2*A2L3s
     &  + 2.d0*mass2b**2*k1p2*A2L2s*A2L3
     &  - 2.d0*mass2b**2*k1p2*A2R2*A2R2s
     &  + 2.d0*mass2b**2*k1p2*A2R2*A2R3s
     &  + 2.d0*mass2b**2*k1p2*A2R2s*A2R3
     &  - 2.d0*mass2b**3*mass4b*A2L2*A2R2s
     &  - 2.d0*mass2b**3*mass4b*A2L2s*A2R2
     &  + 2.d0*mass2b**3*A1L*A2R2s
     &  + 2.d0*mass2b**3*A1R*A2L2s
     &  + 2.d0*mass2b**3*A1Ls*A2R2
     &  + 2.d0*mass2b**3*A1Rs*A2L2
     &  + 2.d0*mass2b**4*A2L2*A2L2s
     &  + 2.d0*mass2b**4*A2R2*A2R2s
     &  - 2.d0*mass4b*A1L*A0Ls*e1p1
     &  - 2.d0*mass4b*A1L*p1p2*A2L1s
     &  - 2.d0*mass4b*A1L*k1p2*A2L3s
     &  - 2.d0*mass4b*A1R*A0Rs*e1p1
     &  - 2.d0*mass4b*A1R*p1p2*A2R1s
     &  - 2.d0*mass4b*A1R*k1p2*A2R3s
     &  - 2.d0*mass4b*A1Ls*A0L*e1p1
     &  - 2.d0*mass4b*A1Ls*p1p2*A2L1
     &  - 2.d0*mass4b*A1Ls*k1p2*A2L3
     &  - 2.d0*mass4b*A1Rs*A0R*e1p1
     &  - 2.d0*mass4b*A1Rs*p1p2*A2R1
     &  - 2.d0*mass4b*A1Rs*k1p2*A2R3
     &  + 2.d0*mass1**2*A0L*e1p1*A2L1s
     &  + 2.d0*mass1**2*A0R*e1p1*A2R1s
     &  + 2.d0*mass1**2*A0Ls*e1p1*A2L1
     &  + 2.d0*mass1**2*A0Rs*e1p1*A2R1
     &  + 2.d0*mass1**2*p1p2*A2L1*A2L1s
     &  + 2.d0*mass1**2*p1p2*A2R1*A2R1s
     &  + 2.d0*mass1**2*k1p2*A2L1*A2L1s
     &  + 2.d0*mass1**2*k1p2*A2L1*A2L3s
     &  + 2.d0*mass1**2*k1p2*A2L1s*A2L3
     &  + 2.d0*mass1**2*k1p2*A2R1*A2R1s
     &  + 2.d0*mass1**2*k1p2*A2R1*A2R3s
     &  + 2.d0*mass1**2*k1p2*A2R1s*A2R3
     &  + 2.d0*A1L*A1Ls*p1p2
     &  - 2.d0*A1L*A1Ls*k1p2
     &  + 2.d0*A1R*A1Rs*p1p2
     &  - 2.d0*A1R*A1Rs*k1p2
     &  + 2.d0*A0L*A0Ls*p1p2
     &  - 2.d0*A0L*A0Ls*k1p2
     &  + 4.d0*A0L*p1p2*e1p1*A2L1s
     &  - 2.d0*A0L*e1p1*k1p1*A2L1s
     &  + 2.d0*A0L*e1p1*k1p1*A2L3s
     &  - 2.d0*A0L*e1p1*k1p2*A2L1s
     &  + 2.d0*A0L*e1p1*k1p2*A2L3s
     &  + 2.d0*A0R*A0Rs*p1p2
     &  - 2.d0*A0R*A0Rs*k1p2
     &  + 4.d0*A0R*p1p2*e1p1*A2R1s
     &  - 2.d0*A0R*e1p1*k1p1*A2R1s
     &  + 2.d0*A0R*e1p1*k1p1*A2R3s
     &  - 2.d0*A0R*e1p1*k1p2*A2R1s
     &  + 2.d0*A0R*e1p1*k1p2*A2R3s
     &  + 4.d0*A0Ls*p1p2*e1p1*A2L1
     &  - 2.d0*A0Ls*e1p1*k1p1*A2L1
     &  + 2.d0*A0Ls*e1p1*k1p1*A2L3
     &  - 2.d0*A0Ls*e1p1*k1p2*A2L1
     &  + 2.d0*A0Ls*e1p1*k1p2*A2L3
     &  + 4.d0*A0Rs*p1p2*e1p1*A2R1
     &  - 2.d0*A0Rs*e1p1*k1p1*A2R1
     &  + 2.d0*A0Rs*e1p1*k1p1*A2R3
     &  - 2.d0*A0Rs*e1p1*k1p2*A2R1
     &  + 2.d0*A0Rs*e1p1*k1p2*A2R3
     &  - 4.d0*p1p2*k1p1*A2L1*A2L1s
     &  - 4.d0*p1p2*k1p1*A2R1*A2R1s
     &  + 4.d0*p1p2*k1p2*A2L1*A2L3s
     &  + 4.d0*p1p2*k1p2*A2L1s*A2L3
     &  + 4.d0*p1p2*k1p2*A2R1*A2R3s
     &  + 4.d0*p1p2*k1p2*A2R1s*A2R3
     &  + 4.d0*p1p2**2*A2L1*A2L1s
     &  + 4.d0*p1p2**2*A2R1*A2R1s
     &  + 4.d0*k1p1*k1p2*A2L3*A2L3s
     &  + 4.d0*k1p1*k1p2*A2R3*A2R3s
     &  + 4.d0*k1p2**2*A2L3*A2L3s
     &  + 4.d0*k1p2**2*A2R3*A2R3s
*****
      enddo
*****
      if (dble(ampl2).lt.0.0d0) then     
        if(aszeroprint) then
        write(*,*) ' '
        write(*,*) 
     &  'ERROR IN dsaschizero with negative cross section:'
        write(*,*) 'Model: ',idtag
        write(*,*) 'kp1=',kp1,'  ',pname(kp1)
        write(*,*) 'kp2=',kp2,'  ',pname(kp2)
        write(*,*) 'kp3=',kp3,'  ',pname(kp3)
        write(*,*) 'kp4=',kp4,'  ',pname(kp4)
        write(*,*) 'amplitude sqaured =',ampl2
        write(*,*) 'mass1=',mass1
        write(*,*) 'mass2=',mass2
        write(*,*) 'mass3=',mass3
        write(*,*) 'mass4=',mass4
        write(*,*) 's=',Svar
        write(*,*) 't=',Tvar
        write(*,*) 'u=',Uvar
        write(*,*) 'p12=',p12
        write(*,*) 'costheta=',costheta
        endif
        ampl2=dcmplx(0.0d0,0.0d0)
      endif
      result=dble(ampl2)*k34/(8.d0*pi*gg1c*gg2c*s34*dsqrt(Svar))
     &       *colfact
      return
      end

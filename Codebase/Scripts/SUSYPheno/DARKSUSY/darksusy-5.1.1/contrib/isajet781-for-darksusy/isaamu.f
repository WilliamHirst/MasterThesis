CDECK  ID>, ISAAMU.
C---------------------------------------------------------------------
       SUBROUTINE ISAAMU (
     ,        SRV2V1 ,SALEMI ,SGAMMAL,SGAMMAR,STWOM1 ,SAAL   ,SAMMLSS,
     ,        SAMMRSS,SAMZ1SS,SAMZ2SS,SAMZ3SS,SAMZ4SS,SAMW1SS,SAMW2SS,
     ,        SAMN2SS,SZMIXSS,iprt,DAMU
     ,                        )
C---------------------------------------------------------------------
      IMPLICIT NONE
      
      REAL*8 DAMU
      
      Real    SRV2V1 ,SALEMI ,SGAMMAL,SGAMMAR,STWOM1,SAAL,      SAMMLSS,
     ,        SAMMRSS,SAMZ1SS,SAMZ2SS,SAMZ3SS,SAMZ4SS,SAMW1SS  ,SAMW2SS,
     ,	      SAMN2SS,SZMIXSS(4,4)
C      SUSY parameters
C      SRV2V1		    = ratio v2/v1 of vev's
c      SALEMI		    = 1/alfa_em
C      SGAMMAL,SGAMMAR      = Wi left, right mixing angles
C      STWOM1		    = Higgsino mass = - mu
C      SAAL		    = stau trilinear term
C      SAMZiSS  	    = signed mass of Zi
C      SAMWiSS  	    = signed Wi mass
C      SAMNiSS  	    = sneutrino mass for generation i
C      SZMIXSS  	    = Zi mixing matrix
c      
      Real*8  RV2V1,ALEMI,GAMMAL,GAMMAR,TWOM1,AAL,AMMLSS,AMMRSS,
     ,        AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,AMW1SS,AMW2SS,AMN2SS,
     ,        ZMIXSS(4,4),ANWI,AMUZI
     
      COMPLEX*16 N1L(4),N1R(4),N2L(4),N2R(4)
 
      REAL*8 AMZISS(4),TH(4),XM1(4),XM2(4),AM1Z(4),AM2Z(4)
 
      COMPLEX*16 ZI
      REAL*8 AMW,AMZ,SN2THW,AMMU
      DATA ZI/(0.d0,1.d0)/,AMW/80.41d0/,AMZ/91.19d0/,SN2THW/.232d0/,
     ,     AMMU/.10566d0/
           
      REAL*8 ALFAEM,PI,SR2,BE,G,GP,FMU,SINB,COSB,TANB,COS2B,XM,YM,THX,
     ,       THY,AMM1SS,AMM2SS,TANTHM,THETAM,COSM,SINM,CW1L,CW2L,CW1R,
     ,       CW2R,AM1ZI,AM2ZI, XN1,XN2,TW1,TW2,ANW1,ANW2
      
      INTEGER I,J,iprt
      real*8 SM11,SM12,SM21,SM22,PNOM1,PNOM2
      
      real*8 amml,ammr,nom1,nom2,den,
     &       MSTAUL2,MSTAUR2,SQRTSTAU,MSTAU12,MSTAU22
      
CsB >
      Real*8 eps 
      eps=0.005
CsB <
     
      RV2V1=SRV2V1
      ALEMI=SALEMI
      GAMMAL=SGAMMAL
      GAMMAR=SGAMMAR
      TWOM1=STWOM1
      AAL=SAAL
      AMMLSS=SAMMLSS
      AMMRSS=SAMMRSS
      
      AMZ1SS=SAMZ1SS
      AMZ2SS=SAMZ2SS
      AMZ3SS=SAMZ3SS
      AMZ4SS=SAMZ4SS
      AMW1SS=SAMW1SS
      AMW2SS=SAMW2SS
      AMN2SS=SAMN2SS
      
      AMW=AMZ*DSQRT(1.D0-SN2THW)
      
      DO I=1,4
        DO J=1,4
          ZMIXSS(I,J)=SZMIXSS(I,J)
	ENDDO
      ENDDO
      
c----------------------------------------
      
      ALFAEM=1.d0/ALEMI
      PI=4.d0*DATAN(1.d0)
      SR2=DSQRT(2.d0)
      BE=DATAN(1.d0/RV2V1)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*DSQRT(SN2THW/(1.d0-SN2THW))
      FMU=G*AMMU/SR2/AMW/DCOS(BE)
      SINB=DSIN(BE)
      COSB=DCOS(BE)
      TANB=SINB/COSB
      COS2B=DCOS(2.d0*BE)
      XM=1.d0/TAN(GAMMAL)
      YM=1.d0/TAN(GAMMAR)
      THX=DSIGN(1.d0,XM)
      THY=DSIGN(1.d0,YM)
      AMZISS(1)=AMZ1SS
      AMZISS(2)=AMZ2SS
      AMZISS(3)=AMZ3SS
      AMZISS(4)=AMZ4SS
     
    
c----- MSTAUL2,MSTAUR2 are 1loop corrected values of soft-masses squared
c----- we just resonstruct them back to use in Moroi formulae
 
      MSTAUL2=AMMLSS**2-AMMU**2-AMZ**2*COS2B*(-.5+SN2THW)
      MSTAUR2=AMMRSS**2-AMMU**2-AMZ**2*COS2B*(-SN2THW)
c------------      

      SQRTSTAU=sqrt( 
     &((MSTAUL2-MSTAUR2)/2d0-COS2B*(4d0*AMW**2-3d0*AMZ**2)/4d0)**2
     &+AMMU**2*(-AAL-TWOM1*SINB/COSB)**2
     &)
      
      MSTAU12=(MSTAUL2+MSTAUR2)/2d0-AMZ**2*COS2B/4d0+AMMU**2-SQRTSTAU
      MSTAU22=(MSTAUL2+MSTAUR2)/2d0-AMZ**2*COS2B/4d0+AMMU**2+SQRTSTAU
      
      AMM1SS=sqrt(min(MSTAU12,MSTAU22))
      AMM2SS=sqrt(max(MSTAU12,MSTAU22))
      
      nom1=(AMM1SS**2-AMMU**2+AMZ**2*COS2B*(.5-SN2THW)-MSTAUL2)
     
      den=AMMU*(TWOM1*SINB/COSB+AAL)
      
      TANTHM =NOM1/den
       
      THETAM=DATAN(TANTHM)

      COSM=   DCOS(THETAM)
      SINM=   DSIN(THETAM)
          
      if(IPRT.ge.4) then
        print *,'AMML=',AMMLSS,' AMM1=',AMM1SS
        print *,'AMMR=',AMMRSS,' AMM2=',AMM2SS
        print *,'tan=', TANTHM,THETAM
        print *,'COSM=',COSM, ',' , 'SINM=',SINM
      ENDIF
      
      CW1L=-FMU*COS(GAMMAL)
      CW2L=FMU*THX*SIN(GAMMAL)
      CW1R=-SIGN(1d0,AMW1SS)*G*SIN(GAMMAR)
      CW2R=-SIGN(1d0,AMW2SS)*G*THY*COS(GAMMAR)
      
      DO I=1,4
      
        IF (AMZISS(I).GT.0.) THEN
         TH(I)=0.d0
        ELSE
         TH(I)=1.d0
        END IF

      N1L(I)=(-ZI)**TH(I)*
     ,(-FMU*ZMIXSS(2,I)*COSM+SR2*GP*ZMIXSS(4,I)*SINM)
    
      N1R(I)=ZI**TH(I)*((G/SR2*ZMIXSS(3,I)+GP/SR2*ZMIXSS(4,I))*COSM
     ,+FMU*ZMIXSS(2,I)*SINM)
    
      N2L(I)=-(-ZI)**TH(I)*
     ,(FMU*ZMIXSS(2,I)*SINM+SR2*GP*ZMIXSS(4,I)*COSM)
    
      N2R(I)=ZI**TH(I)*((G/SR2*ZMIXSS(3,I)+GP/SR2*ZMIXSS(4,I))*SINM
     ,-FMU*ZMIXSS(2,I)*COSM)

      
      XM1(I)=(AMZISS(I)/AMM1SS)**2
      XM2(I)=(AMZISS(I)/AMM2SS)**2

      IF (Abs(XM1(I)-1.d0).LT.eps)
     &  XM1(I)=1.d0+Sign(1.d0,XM2(I)-1.d0)*eps
      IF (Abs(XM2(I)-1.d0).LT.eps)
     &	XM2(I)=1.d0+Sign(1.d0,XM2(I)-1.d0)*eps


      if(iprt.ge.4) THEN
        print *,'==========================='
        print *,'I=',I
        print '(3A15)'	  ,'ZMIXSS(2,I)','ZMIXSS(3,I)','ZMIXSS(4,I)'
        print '(3D15.4)'    , ZMIXSS(2,I),  ZMIXSS(3,I),  ZMIXSS(4,I)
        print *,'--'
        print '(2A20)',     'N1L(I)','N1R(I)'
        print * ,    N1L(I) , N1R(I) 
        print * ,'--'
        print '(2A20)',     'N2L(I)','N2R(I)'
        print * ,    N2L(I) , N2R(I)
        print * ,'==========================='
      endif 
      END DO     
      
C======NEUTRALINO LOOPS================================================
      AM1ZI=0.d0
      AM2ZI=0.d0
      
  
      DO I=1,4
        AM1Z(I)=
     &-AMMU/(6.d0*AMM1SS**2*(1.d0-XM1(I))**4)
     & *(N1L(I)*DCONJG(N1L(I))+N1R(I)*DCONJG(N1R(I)))
     & *(1.d0-6*XM1(I)+
     &   3.d0*XM1(I)**2+2.d0*XM1(I)**3-6.d0*XM1(I)**2*LOG(XM1(I)))
     &-ABS(AMZISS(I))/(AMM1SS**2*(1.d0-XM1(I))**3)
     & *N1L(I)*DCONJG(N1R(I))*(1.d0-XM1(I)**2+2.d0*XM1(I)*LOG(XM1(I)))

        AM1ZI=AM1ZI+AM1Z(I)

        AM2Z(I)=
     &-AMMU/(6.d0*AMM2SS**2*(1.d0-XM2(I))**4)
     & *(N2L(I)*DCONJG(N2L(I))+N2R(I)*DCONJG(N2R(I)))
     & *(1.d0-6.d0*XM2(I)+
     &   3.d0*XM2(I)**2+2.d0*XM2(I)**3-6.d0*XM2(I)**2*LOG(XM2(I)))
     &-ABS(AMZISS(I))/(AMM2SS**2*(1.d0-XM2(I))**3)
     & *N2L(I)*DCONJG(N2R(I))*(1.d0-XM2(I)**2+2.d0*XM2(I)*LOG(XM2(I)))

        AM2ZI=AM2ZI+AM2Z(I)
      
      if(iprt.ge.4) THEN
        print *,'I, AM1Z(I), AM2Z(I):', I, AM1Z(I), AM2Z(I)
      endif 
      
c      print *,I, AMM1SS, AMM2SS, AMZISS(I), XM1(I),XM2(I)
c      print *,I, N1L(I)*DCONJG(N1L(I)), N1R(I)*DCONJG(N1R(I))
c      print *,I, N1L(I)*DCONJG(N1R(I))
c      print *,I, N2L(I)*DCONJG(N2L(I)), N2R(I)*DCONJG(N2R(I))
c      print *,I, N2L(I)*DCONJG(N2R(I)) 

      END DO

      AMUZI=(AM1ZI+AM2ZI)/16.d0/PI**2*AMMU
      
      if(iprt.ge.4) THEN
        print *,'AM1ZI, AM2ZI,AMUZI :',
     &          (AM1ZI)/16.d0/PI**2*AMMU,
     &          (AM2ZI)/16.d0/PI**2*AMMU,
     &          (AM1ZI+AM2ZI)/16.d0/PI**2*AMMU
      endif 
 
c=======CHARGINO LOOPS=================================================
     
      XN1=(AMW1SS/AMN2SS)**2
      XN2=(AMW2SS/AMN2SS)**2

      IF (Abs(XN1-1.d0).LT.eps) XN1=1.d0+Sign(1.d0,XN1-1.d0)*eps
      IF (Abs(XN2-1.d0).LT.eps) XN2=1.d0+Sign(1.d0,XN2-1.d0)*eps

    
      TW1=AMMU/(3.d0*AMN2SS**2*(1.d0-XN1)**4)*(CW1L*CW1L+CW1R*CW1R)*
     ,(1.d0+1.5d0*XN1-3.d0*XN1**2+.5d0*XN1**3+3.d0*XN1*DLOG(XN1))
     
      TW2=-(3.d0*ABS(AMW1SS))/(AMN2SS**2*(1.d0-XN1)**3)*CW1L*CW1R*
     ,(1.d0-4.d0*XN1/3.d0+XN1**2/3.d0+2.d0*DLOG(XN1)/3.d0)
     
      ANW1=TW1+TW2
      
      TW1=AMMU/(3.d0*AMN2SS**2*(1.d0-XN2)**4)*(CW2L*CW2L+CW2R*CW2R)*
     ,(1.d0+1.5d0*XN2-3.d0*XN2**2+.5d0*XN2**3+3.d0*XN2*DLOG(XN2))

      TW2=-(3.d0*ABS(AMW2SS))/(AMN2SS**2*(1.d0-XN2)**3)*CW2L*CW2R*
     ,(1.d0-4.d0*XN2/3.d0+XN2**2/3.d0+2.d0*DLOG(XN2)/3.d0)
      ANW2=TW1+TW2

      if(iprt.ge.4) THEN
        print *,'ANW1  ,ANW2  ,ANW  :', 
     &           ANW1/16.d0/PI**2*AMMU,ANW2/16.d0/PI**2*AMMU,
     &           ANW1/16.d0/PI**2*AMMU+ANW2/16.d0/PI**2*AMMU
      endif
      
      ANWI=(ANW1+ANW2)/16.d0/PI**2*AMMU
      DAMU=ANWI+AMUZI

      if(iprt.ge.4) THEN
      print*, 'AMU=',ANWI+AMUZI
      print *,'==========================================='
      print *,'==========================================='
      endif
      
      Return
      End

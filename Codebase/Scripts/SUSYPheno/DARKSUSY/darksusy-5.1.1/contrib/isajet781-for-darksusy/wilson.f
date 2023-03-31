C---------------------------------------------------------------------- 
      SUBROUTINE WILSON(GK,HK,
     &                  CSM,CH1,CH2,CC111,CC112,CC113,CC121,CC122,CC123,
     &                  CC211,CC212,CC213,CC221,CC222,CC223,
     &                  DSSM,DSH1,DSH2,DSC111,DSC112,DSC113,DSC121,
     &                  DSC122,DSC123,DSC211,DSC212,DSC213,DSC221,
     &                  DSC222,DSC223)
C----------------------------------------------------------------------
C
C    Calculates Wilson coefficients and contributions from decoupling
C    of heavy particles at the scale of their mass.
C
      IMPLICIT NONE
      REAL*8 GK(2,3),HK(2,3)
      REAL*8 CSM(9),CH1(9),CH2(19),
     $       CC111(12),CC112(12),CC113(12),
     $       CC121(12),CC122(12),CC123(12),
     $       CC211(17),CC212(17),CC213(17),
     $       CC221(17),CC222(17),CC223(17),      
     $       DSSM(9),DSH1(9),DSH2(9),
     $       DSC111(9),DSC112(9),DSC113(9),
     $       DSC121(9),DSC122(9),DSC123(9),
     $       DSC211(9),DSC212(9),DSC213(9),
     $       DSC221(9),DSC222(9),DSC223(9)
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX
      SAVE /CHRGN/        
      REAL*8 X,FN1,FN2,FN3,FN4,SCALE,XWA,XWB,XCB,YPS
	
      CSM(1)=-1./2.d0
      CSM(2)=-1./2.d0
      CSM(3)= 1.d0
      CSM(4)= 11./18.d0
      CSM(5)=-8./9.d0
      CSM(6)= 11./18.d0
      CSM(7)= 1./2.d0
      CSM(8)=-9./4.d0
      CSM(9)= 3./2.d0
c
      X=(MT/MW)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSSM(1)=-X*FN4+1./2.d0
      DSSM(2)=-X/2.*(FN3+FN4)+1./2.d0	     
      DSSM(3)=-2.*DSSM(2)
      DSSM(4)= (X+2.d0)/3.*(2.*FN2+FN3+2.*FN4)-11./18.d0	
      DSSM(5)= 2.*(X+2.d0)/3.*(FN2-FN3-2.*FN4)+8./9.d0   
      DSSM(6)= (X+2.d0)/3.*(2.*FN2+FN3+2.*FN4)-11./18.d0	
      DSSM(7)= (X-2.d0)*FN4-1./2.d0   
      DSSM(8)=-3.*((X+2.d0)*(1./2.*FN3+FN4-LOG(X)/6./(X-1.d0))
     $  	   -3./4.d0)	    
      DSSM(9)=-3.*(7./2.d0-2.*X*FN3-5.*X*FN4)
C
      CH1(1)= 1./2.d0
      CH1(2)= 1./2.d0
      CH1(3)=-1.d0
      CH1(4)= 11./18.d0/TANB**2
      CH1(5)=-8./9.d0/TANB**2
      CH1(6)= 11./18.d0/TANB**2
      CH1(7)= 1./2.d0/TANB**2
      CH1(8)=-9./4.d0/TANB**2
      CH1(9)= 3./2.d0/TANB**2	  
c
      X=(MT/MHPLUS)**2
      CH2(1)= 1./2.d0*X
      CH2(2)=-1./2.d0*X
      CH2(3)= X
      CH2(4)=-1./9.d0*X/TANB**2
      CH2(5)= 7./18.d0*X/TANB**2
      CH2(6)=-1./9.d0*X/TANB**2
      CH2(7)= 1./2.d0*X/TANB**2
      CH2(8)= 3./4.d0*X/TANB**2
      CH2(9)= 3./2.d0*X/TANB**2     
      CH2(10)=0.d0
      CH2(12)=0.d0
      CH2(13)=0.d0
      CH2(14)=0.d0
      CH2(15)=0.d0    
      CH2(16)=0.d0
      CH2(18)=0.d0
      CH2(19)=0.d0
c
      X=(MT/MHPLUS)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSH1(1)= X*FN4-1./2.d0
      DSH1(2)= X/2.*(FN3+FN4)-1./2.d0	     
      DSH1(3)=-2.*DSH1(2)
      DSH1(4)= (X/3.d0*(2.d0*FN2+FN3+2*FN4)-11./18.d0)/TANB**2        
      DSH1(5)= (2.*X/3.d0*(FN2-FN3-2*FN4)+8./9.d0)/TANB**2   
      DSH1(6)= (X/3.d0*(2.*FN2+FN3+2*FN4)-11./18.d0)/TANB**2	    
      DSH1(7)= (X*FN4-1./2.d0)/TANB**2   
      DSH1(8)=-3.*(X*(1./2.*FN3+FN4-LOG(X)/6./(X-1.d0))
     $  	   -3./4.d0)/TANB**2	    
      DSH1(9)=3.*DSH1(7)
c
      X=(MT/MHPLUS)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSH2(1)=X*(FN4-1./2.d0)
      DSH2(2)=X/2.d0*(FN3+FN4+1.d0)	   
     $        +(MT/MHPLUS)**2*LOG(MT/MHPLUS)
      DSH2(3)=-2.*DSH2(2)
c              -2.*(MT/MHPLUS)**2*LOG(MT/MHPLUS)
      DSH2(4)=(X/3.d0*(2*FN2+FN3+2*FN4+1./3.d0))/TANB**2	
      DSH2(5)=(2*X/3.d0*(FN2-FN3-2*FN4)-7./18.d0*X)/TANB**2   
      DSH2(6)=(X/3.d0*(2*FN2+FN3+2*FN4+1./3.d0))/TANB**2	
      DSH2(7)=X*(FN4-1./2.d0)/TANB**2	
      DSH2(8)=-3.*(X*(1./2.d0*FN3+FN4)-LOG(X)/6.d0/(X-1.d0)
     $  	   +1./4.d0*X)/TANB**2        
      DSH2(9)=3.*DSH2(7)		

      IF(MCHA(1).GT.MW) THEN
       SCALE=MCHA(1)
      ELSE
       SCALE=MW
      ENDIF
C
      XWA=(MW/MSQU(1))**2
      XCB=MCHA(1)/MBQ
      CC111(1)=-GK(1,1)*HK(1,1)*XCB*XWA
      CC111(2)= 0.d0
      CC111(3)= 3.*GK(1,1)*HK(1,1)*XCB*XWA
      CC111(4)= 5./18.d0*GK(1,1)*GK(1,1)*XWA
      CC111(5)=-2./9.d0*GK(1,1)*GK(1,1)*XWA
      CC111(6)= 5./18.d0*GK(1,1)*GK(1,1)*XWA
      CC111(7)= 0.d0
      CC111(8)=-3./2.d0*GK(1,1)*GK(1,1)*XWA
      CC111(9)=-3.*GK(1,1)*GK(1,1)*XWA
c
      X=(MCHA(1)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
        DSC111(1)=-GK(1,1)*HK(1,1)*XCB*XWA*(2.*FN4-1.d0)
        DSC111(2)= 0.d0
        DSC111(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(1))) 
        DSC111(4)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
        DSC111(5)=2./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
        DSC111(6)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
        DSC111(7)=0.d0
        DSC111(8)=6.*GK(1,1)*GK(1,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(1)))
        DSC111(9)=-6.*GK(1,1)*GK(1,1)*XWA*(FN4-1./2.d0)
      ELSE
        DSC111(1)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC111(1)
        DSC111(2)=0.d0
        DSC111(3)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(1)))-CC111(3) 
        DSC111(4)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC111(4) 
        DSC111(5)=GK(1,1)*GK(1,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC111(5)
        DSC111(6)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC111(6) 
        DSC111(7)=0.d0
        DSC111(8)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(1)))-CC111(8)
        DSC111(9)=GK(1,1)*GK(1,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC111(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC211(1)=-GK(1,1)*HK(1,1)*XCB*XWB
      CC211(2)= 0.d0
      CC211(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWB
      CC211(4)=-5./18.d0*GK(1,1)*GK(1,1)*XWB
      CC211(5)= 11./9.d0*GK(1,1)*GK(1,1)*XWB
      CC211(6)=-5./18.d0*GK(1,1)*GK(1,1)*XWB
      CC211(7)= 0.d0
      CC211(8)= 9./2.d0*GK(1,1)*GK(1,1)*XWB
      CC211(9)=-3.*GK(1,1)*GK(1,1)*XWB      
      CC211(10)=-2.*GK(1,1)*HK(1,1)*XCB*XWB	 
      CC211(11)=0.d0
      CC211(12)=2.*GK(1,1)*GK(1,1)*XWB      
      CC211(13)=2.*GK(1,1)*GK(1,1)*XWB      
      CC211(14)=0.d0
      CC211(15)=0.d0   
      CC211(16)=0.d0
      CC211(17)=0.d0
C
      X=(MCHA(1)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC211(1)=-GK(1,1)*HK(1,1)*XCB*XWA*2.*FN4-CC211(1)
       DSC211(2)= 0.d0
       DSC211(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWA*
     $  	    (FN3+FN4)-CC211(3)
       DSC211(4)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC211(4)
       DSC211(5)=2./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC211(5)
       DSC211(6)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC211(6)
       DSC211(7)=0.d0
       DSC211(8)= 6.*GK(1,1)*GK(1,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC211(8)
       DSC211(9)=-6.*GK(1,1)*GK(1,1)*XWA*FN4-CC211(9)		       
      ELSE
       DSC211(1)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC211(1)
       DSC211(2)=0.d0
       DSC211(3)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC211(3) 
       DSC211(4)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC211(4) 
       DSC211(5)=GK(1,1)*GK(1,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC211(5)
       DSC211(6)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC211(6) 
       DSC211(7)=0.d0
       DSC211(8)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC211(8)
       DSC211(9)=GK(1,1)*GK(1,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC211(9)
      ENDIF
C
      XWA=(MW/MSQU(2))**2
      XCB=MCHA(1)/MBQ
      CC112(1)=-GK(1,2)*HK(1,2)*XCB*XWA
      CC112(2)= 0.d0
      CC112(3)= 3.*GK(1,2)*HK(1,2)*XCB*XWA
      CC112(4)= 5./18.d0*GK(1,2)*GK(1,2)*XWA
      CC112(5)=-2./9.d0*GK(1,2)*GK(1,2)*XWA
      CC112(6)= 5./18.d0*GK(1,2)*GK(1,2)*XWA
      CC112(7)= 0.d0
      CC112(8)=-3./2.d0*GK(1,2)*GK(1,2)*XWA
      CC112(9)=-3.*GK(1,2)*GK(1,2)*XWA         
c
      X=(MCHA(1)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC112(1)=-GK(1,2)*HK(1,2)*XCB*XWA*(2.*FN4-1.d0)
       DSC112(2)= 0.d0
       DSC112(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(2)))
       DSC112(4)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC112(5)=2./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC112(6)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC112(7)=0.d0
       DSC112(8)=6.*GK(1,2)*GK(1,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(2)))
       DSC112(9)=-6.*GK(1,2)*GK(1,2)*XWA*(FN4-1./2.d0)
      ELSE
       DSC112(1)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC112(1)
       DSC112(2)=0.d0
       DSC112(3)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(2)))-CC112(3) 
       DSC112(4)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC112(4) 
       DSC112(5)=GK(1,2)*GK(1,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC112(5)
       DSC112(6)=GK(1,1)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC112(6) 
       DSC112(7)=0.d0
       DSC112(8)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(2)))-CC112(8)
       DSC112(9)=GK(1,2)*GK(1,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC112(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC212(1)=-GK(1,2)*HK(1,2)*XCB*XWB
      CC212(2)= 0.d0
      CC212(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWB
      CC212(4)=-5./18.d0*GK(1,2)*GK(1,2)*XWB
      CC212(5)= 11./9.d0*GK(1,2)*GK(1,2)*XWB
      CC212(6)=-5./18.d0*GK(1,2)*GK(1,2)*XWB
      CC212(7)= 0.d0
      CC212(8)= 9./2.d0*GK(1,2)*GK(1,2)*XWB
      CC212(9)=-3.*GK(1,2)*GK(1,2)*XWB      
      CC212(10)=-2.*GK(1,2)*HK(1,2)*XCB*XWB	 
      CC212(11)=0.d0
      CC212(12)=2.*GK(1,2)*GK(1,2)*XWB      
      CC212(13)=2.*GK(1,2)*GK(1,2)*XWB      
      CC212(14)=0.d0
      CC212(15)=0.d0	
      CC212(16)=0.d0
      CC212(17)=0.d0
C
      X=(MCHA(1)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC212(1)=-GK(1,2)*HK(1,2)*XCB*XWA*2.*FN4-CC212(1)
       DSC212(2)=0.d0
       DSC212(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWA*
     $  	    (FN3+FN4)-CC212(3)
       DSC212(4)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC212(4)
       DSC212(5)=2./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC212(5)
       DSC212(6)=1./3.*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC212(6)
       DSC212(7)=0.d0
       DSC212(8)=6.*GK(1,2)*GK(1,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC212(8)
       DSC212(9)=-6.*GK(1,2)*GK(1,2)*XWA*FN4-CC212(9)	       
      ELSE
       DSC212(1)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC212(1)
       DSC212(2)=0.d0
       DSC212(3)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC212(3) 
       DSC212(4)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC212(4) 
       DSC212(5)=GK(1,2)*GK(1,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC212(5)
       DSC212(6)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC212(6) 
       DSC212(7)=0.d0
       DSC212(8)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC212(8)
       DSC212(9)=GK(1,2)*GK(1,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC212(9)
      ENDIF
c
      XWA=(MW/MSQU(3))**2
      XCB=MCHA(1)/MBQ
      CC113(1)=-HK(1,3)*XCB*XWA
      CC113(2)= 0.d0
      CC113(3)= 3.*HK(1,3)*XCB*XWA
      CC113(4)= 5./18.d0*GK(1,3)*XWA
      CC113(5)=-2./9.d0*GK(1,3)*XWA
      CC113(6)= 5./18.d0*GK(1,3)*XWA
      CC113(7)= 0.d0
      CC113(8)=-3./2.d0*GK(1,3)*XWA
      CC113(9)=-3.*GK(1,3)*XWA         
c
      X=(MCHA(1)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC113(1)=-HK(1,3)*XCB*XWA*(2.*FN4-1.d0)
       DSC113(2)=0.d0
       DSC113(3)=-3.*HK(1,3)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(3)))
       DSC113(4)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC113(5)=2./3.d0*GK(1,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC113(6)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC113(7)=0.d0
       DSC113(8)=6.*GK(1,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(3)))
       DSC113(9)=-6.*GK(1,3)*XWA*(FN4-1./2.d0)
      ELSE
       DSC113(1)=HK(1,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC113(1)
       DSC113(2)=0.d0
       DSC113(3)=HK(1,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(3)))-CC113(3) 
       DSC113(4)=GK(1,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC113(4) 
       DSC113(5)=GK(1,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC113(5)
       DSC113(6)=GK(1,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC113(6) 
       DSC113(7)=0.d0
       DSC113(8)=GK(1,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(3)))-CC113(8)
       DSC113(9)=GK(1,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC113(9)
      ENDIF
C       								     
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC213(1)=-HK(1,3)*XCB*XWB
      CC213(2)= 0.d0
      CC213(3)=-3.*HK(1,3)*XCB*XWB
      CC213(4)=-5./18.d0*GK(1,3)*XWB
      CC213(5)= 11./9.d0*GK(1,3)*XWB
      CC213(6)=-5./18.d0*GK(1,3)*XWB
      CC213(7)= 0.d0
      CC213(8)= 9./2.d0*GK(1,3)*XWB
      CC213(9)=-3.*GK(1,3)*XWB      
      CC213(10)=-2.*HK(1,3)*XCB*XWB	 
      CC213(11)=0.d0
      CC213(12)=2.*GK(1,3)*XWB      
      CC213(13)=2.*GK(1,3)*XWB      
      CC213(14)=0.d0
      CC213(15)=0.d0	
      CC213(16)=0.d0
      CC213(17)=0.d0
C
      X=(MCHA(1)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC213(1)=-HK(1,3)*XCB*XWA*2.*FN4-CC213(1)
       DSC213(2)=0.d0
       DSC213(3)=-3.*HK(1,3)*XCB*XWA*
     $  	    (FN3+FN4)-CC213(3)
       DSC213(4)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC213(4)
       DSC213(5)=2./3.d0*GK(1,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC213(5)
       DSC213(6)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC213(6)
       DSC213(7)=0.d0
       DSC213(8)=6.*GK(1,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC213(8)
       DSC213(9)=-6.*GK(1,3)*XWA*FN4-CC213(9)	       
      ELSE
       DSC213(1)=HK(1,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC213(1)
       DSC213(2)=0.d0
       DSC213(3)=HK(1,3)*XCB*XWA*
     $  	 (-3./2.d0+YPS)-CC213(3) 
       DSC213(4)=GK(1,3)*XWA*
     $  	 (1./6.d0-1./20.d0*YPS)-CC213(4) 
       DSC213(5)=GK(1,3)*XWA*
     $  	 (-1./6.d0+1./30.d0*YPS)-CC213(5)
       DSC213(6)=GK(1,3)*XWA*
     $  	 (1./6.d0-1./20.d0*YPS)-CC213(6) 
       DSC213(7)=0.
       DSC213(8)=GK(1,3)*XWA*
     $  	 (1.d0-3./4.d0*YPS)-CC213(8)
       DSC213(9)=GK(1,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC213(9)
      ENDIF
c
      IF(MCHA(2).GT.MW) THEN
       SCALE=MCHA(2)
      ELSE
       SCALE=MW
      ENDIF
C
      XWA=(MW/MSQU(1))**2
      XCB=MCHA(2)/MBQ
      CC121(1)=-GK(2,1)*HK(2,1)*XCB*XWA
      CC121(2)= 0.d0
      CC121(3)= 3.*GK(2,1)*HK(2,1)*XCB*XWA
      CC121(4)= 5./18.d0*GK(2,1)*GK(2,1)*XWA
      CC121(5)=-2./9.d0*GK(2,1)*GK(2,1)*XWA
      CC121(6)= 5./18.d0*GK(2,1)*GK(2,1)*XWA
      CC121(7)= 0.d0
      CC121(8)=-3./2.d0*GK(2,1)*GK(2,1)*XWA
      CC121(9)=-3.*GK(2,1)*GK(2,1)*XWA         
c
      X=(MCHA(2)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC121(1)=-GK(2,1)*HK(2,1)*XCB*XWA*(2.*FN4-1.d0)
       DSC121(2)=0.
       DSC121(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(1)))
       DSC121(4)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC121(5)=2./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC121(6)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC121(7)=0.d0
       DSC121(8)=6.*GK(2,1)*GK(2,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(1)))
       DSC121(9)=-6.*GK(2,1)*GK(2,1)*XWA*(FN4-1./2.d0)
      ELSE
       DSC121(1)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC121(1)
       DSC121(2)=0.d0
       DSC121(3)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(1)))-CC121(3) 
       DSC121(4)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC121(4) 
       DSC121(5)=GK(2,1)*GK(2,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC121(5)
       DSC121(6)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC121(6) 
       DSC121(7)=0.d0
       DSC121(8)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(1)))-CC121(8)
       DSC121(9)=GK(2,1)*GK(2,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC121(9)
      ENDIF
C       							       
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC221(1)=-GK(2,1)*HK(2,1)*XCB*XWB
      CC221(2)= 0.d0
      CC221(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWB
      CC221(4)=-5./18.d0*GK(2,1)*GK(2,1)*XWB
      CC221(5)= 11./9.d0*GK(2,1)*GK(2,1)*XWB
      CC221(6)=-5./18.d0*GK(2,1)*GK(2,1)*XWB
      CC221(7)= 0.d0
      CC221(8)= 9./2.d0*GK(2,1)*GK(2,1)*XWB
      CC221(9)=-3.*GK(2,1)*GK(2,1)*XWB      
      CC221(10)=-2.*GK(2,1)*HK(2,1)*XCB*XWB	 
      CC221(11)=0.d0
      CC221(12)=2.*GK(2,1)*GK(2,1)*XWB      
      CC221(13)=2.*GK(2,1)*GK(2,1)*XWB      
      CC221(14)=0.d0
      CC221(15)=0.d0	
      CC221(16)=0.d0
      CC221(17)=0.d0
C
      X=(MCHA(2)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC221(1)=-GK(2,1)*HK(2,1)*XCB*XWA*2.*FN4-CC221(1)
       DSC221(2)=0.
       DSC221(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWA*
     $  	    (FN3+FN4)-CC221(3)
       DSC221(4)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC221(4)
       DSC221(5)=2./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC221(5)
       DSC221(6)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC221(6)
       DSC221(7)=0.d0
       DSC221(8)=6.*GK(2,1)*GK(2,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC221(8)
       DSC221(9)=-6.*GK(2,1)*GK(2,1)*XWA*FN4-CC221(9)		       
      ELSE
       DSC221(1)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC221(1)
       DSC221(2)=0.d0
       DSC221(3)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC221(3) 
       DSC221(4)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC221(4) 
       DSC221(5)=GK(2,1)*GK(2,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC221(5)
       DSC221(6)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC221(6) 
       DSC221(7)=0.d0
       DSC221(8)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC221(8)
       DSC221(9)=GK(2,1)*GK(2,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC221(9)
      ENDIF
C
      XWA=(MW/MSQU(2))**2
      XCB=MCHA(2)/MBQ
      CC122(1)=-GK(2,2)*HK(2,2)*XCB*XWA
      CC122(2)= 0.d0
      CC122(3)= 3.*GK(2,2)*HK(2,2)*XCB*XWA
      CC122(4)= 5./18.d0*GK(2,2)*GK(2,2)*XWA
      CC122(5)=-2./9.d0*GK(2,2)*GK(2,2)*XWA
      CC122(6)= 5./18.d0*GK(2,2)*GK(2,2)*XWA
      CC122(7)= 0.d0
      CC122(8)=-3./2.d0*GK(2,2)*GK(2,2)*XWA
      CC122(9)=-3.*GK(2,2)*GK(2,2)*XWA         
c
      X=(MCHA(2)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC122(1)=-GK(2,2)*HK(2,2)*XCB*XWA*(2.*FN4-1.d0)
       DSC122(2)=0.d0
       DSC122(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (FN3+FN4+1.+2.*LOG(SCALE/MSQU(2)))
       DSC122(4)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC122(5)=2./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC122(6)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC122(7)=0.d0
       DSC122(8)=6.*GK(2,2)*GK(2,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(2)))
       DSC122(9)=-6.*GK(2,2)*GK(2,2)*XWA*(FN4-1./2.d0)
      ELSE
       DSC122(1)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC122(1)
       DSC122(2)=0.d0
       DSC122(3)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	 (-3./2.d0+YPS
     $  	  -6.*LOG(SCALE/MSQU(2)))-CC122(3) 
       DSC122(4)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC122(4) 
       DSC122(5)=GK(2,2)*GK(2,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC122(5)
       DSC122(6)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC122(6) 
       DSC122(7)=0.d0
       DSC122(8)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(2)))-CC122(8)
       DSC122(9)=GK(2,2)*GK(2,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC122(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC222(1)=-GK(2,2)*HK(2,2)*XCB*XWB
      CC222(2)= 0.d0
      CC222(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWB
      CC222(4)=-5./18.d0*GK(2,2)*GK(2,2)*XWB
      CC222(5)= 11./9.d0*GK(2,2)*GK(2,2)*XWB
      CC222(6)=-5./18.d0*GK(2,2)*GK(2,2)*XWB
      CC222(7)= 0.d0
      CC222(8)= 9./2.d0*GK(2,2)*GK(2,2)*XWB
      CC222(9)=-3.*GK(2,2)*GK(2,2)*XWB      
      CC222(10)=-2.*GK(2,2)*HK(2,2)*XCB*XWB	 
      CC222(11)=0.d0
      CC222(12)=2.*GK(2,2)*GK(2,2)*XWB      
      CC222(13)=2.*GK(2,2)*GK(2,2)*XWB      
      CC222(14)=0.d0
      CC222(15)=0.d0	
      CC222(16)=0.d0
      CC222(17)=0.d0
C
      X=(MCHA(2)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC222(1)=-GK(2,2)*HK(2,2)*XCB*XWA*2.*FN4-CC222(1)
       DSC222(2)=0.d0
       DSC222(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWA*
     $  	    (FN3+FN4)-CC222(3)
       DSC222(4)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC222(4)
       DSC222(5)=2./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC222(5)
       DSC222(6)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC222(6)
       DSC222(7)=0.d0
       DSC222(8)=6.*GK(2,2)*GK(2,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC222(8)
       DSC222(9)=-6.*GK(2,2)*GK(2,2)*XWA*FN4-CC222(9)	       
      ELSE
       DSC222(1)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC222(1)
       DSC222(2)=0.d0
       DSC222(3)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC222(3) 
       DSC222(4)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.-1./20.*YPS)-CC222(4) 
       DSC222(5)=GK(2,2)*GK(2,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC222(5)
       DSC222(6)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC222(6) 
       DSC222(7)=0.d0
       DSC222(8)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC222(8)
       DSC222(9)=GK(2,2)*GK(2,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC222(9)
      ENDIF
c
      XWA=(MW/MSQU(3))**2
      XCB=MCHA(2)/MBQ
      CC123(1)=-HK(2,3)*XCB*XWA
      CC123(2)= 0.d0
      CC123(3)= 3.*HK(2,3)*XCB*XWA
      CC123(4)= 5./18.d0*GK(2,3)*XWA
      CC123(5)=-2./9.d0*GK(2,3)*XWA
      CC123(6)= 5./18.d0*GK(2,3)*XWA
      CC123(7)= 0.d0
      CC123(8)=-3./2.d0*GK(2,3)*XWA
      CC123(9)=-3.*GK(2,3)*XWA         
c
      X=(MCHA(2)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC123(1)=-HK(2,3)*XCB*XWA*(2.*FN4-1.d0)
       DSC123(2)=0.d0
       DSC123(3)=-3.*HK(2,3)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(3)))
       DSC123(4)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3+2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC123(5)=2./3.d0*GK(2,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC123(6)=1./3.d0*GK(2,3)*XWA*
     $  	  (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC123(7)=0.d0
       DSC123(8)=6.*GK(2,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(3)))
       DSC123(9)=-6.*GK(2,3)*XWA*(FN4-1./2.d0)
      ELSE
       DSC123(1)=HK(2,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC123(1)
       DSC123(2)=0.d0
       DSC123(3)=HK(2,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(3)))-CC123(3) 
       DSC123(4)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC123(4) 
       DSC123(5)=GK(2,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC123(5)
       DSC123(6)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC123(6) 
       DSC123(7)=0.
       DSC123(8)=GK(2,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(3)))-CC123(8)
       DSC123(9)=GK(2,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC123(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC223(1)=-HK(2,3)*XCB*XWB
      CC223(2)=0.d0
      CC223(3)=-3.*HK(2,3)*XCB*XWB
      CC223(4)=-5./18.d0*GK(2,3)*XWB
      CC223(5)=11./9.d0*GK(2,3)*XWB
      CC223(6)=-5./18.d0*GK(2,3)*XWB
      CC223(7)=0.d0
      CC223(8)=9./2.d0*GK(2,3)*XWB
      CC223(9)=-3.*GK(2,3)*XWB      
      CC223(10)=-2.*HK(2,3)*XCB*XWB	 
      CC223(11)=0.d0
      CC223(12)=2.*GK(2,3)*XWB      
      CC223(13)=2.*GK(2,3)*XWB      
      CC223(14)=0.d0
      CC223(15)=0.d0	
      CC223(16)=0.d0
      CC223(17)=0.d0
C
      X=(MCHA(2)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4) 
      IF(ABS(YPS).GT.0.05) THEN
       DSC223(1)=-HK(2,3)*XCB*XWA*2.*FN4-CC223(1)
       DSC223(2)=0.d0
       DSC223(3)=-3.*HK(2,3)*XCB*XWA*
     $  	    (FN3+FN4)-CC223(3)
       DSC223(4)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC223(4)
       DSC223(5)=2./3.d0*GK(2,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC223(5)
       DSC223(6)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC223(6)
       DSC223(7)=0.d0
       DSC223(8)=6.*GK(2,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC223(8)
       DSC223(9)=-6.*GK(2,3)*XWA*FN4-CC223(9)	       
      ELSE
       DSC223(1)=HK(2,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC223(1)
       DSC223(2)=0.d0
       DSC223(3)=HK(2,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC223(3) 
       DSC223(4)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC223(4) 
       DSC223(5)=GK(2,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC223(5)
       DSC223(6)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC223(6) 
       DSC223(7)=0.
       DSC223(8)=GK(2,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC223(8)
       DSC223(9)=GK(2,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC223(9)
      ENDIF

      RETURN
      END

C***********************************************************************
C*====================================                                 *
C* SUBROUTINE SPRGEN( F, MXTRY, NTRY )                                 *
C*====================================                                 *
C*                                                                     *
C*     Generation of events according to the probability density       *
C*     which is stored in a disk file.                                 *
C*                                                                     *
C*    Coded   by S.Kawabata   at July,1980                             *
C*    Update     S.Kawabata   September '84                            *
C*                                                                     *
C***********************************************************************
C                                                                       
       SUBROUTINE SPRGEN(F,MXTRY,NTRY)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      EXTERNAL F                                                        
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)                   
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                    
     .               IG(MXDIM),NCALL                                    
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),      
     .               ND,NG,NPG,MA(MXDIM)                                
                                                                        
      COMMON /SPRNG1/ XND, DXG, XJAC, DXMAX, NSP                        
                                                                        
      DIMENSION Y(MXDIM),KG(MXDIM)                                      
      DATA ONE/1.0D0/                                                   
C                                                                       
C                                                                       
      RX    = DRN(IDUMY)*DXMAX                                           
C                                                                       
C  -------------- Binary Search  --------------------------------       
C                                                                       
      IPMIN = 1                                                         
      IPMAX = NSP                                                       
C                                                                       
 300  IC    = (IPMIN+IPMAX)/2                                           
        IF(RX .LT. DXD(IC)) THEN                                        
          IPMAX = IC                                                    
        ELSE                                                            
          IPMIN = IC                                                    
        ENDIF                                                           
      IF(IPMAX-IPMIN .GT.  2) GO TO 300                                 
C                                                                       
      IC    = IPMIN-1                                                   
 350  IC    = IC+1                                                      
      IF(DXD(IC) .LT. RX) GO TO 350                                     
C                                                                       
C --------------------------------------------------------------------  
C      Identify the hypecube number from sequential number IC           
C --------------------------------------------------------------------  
C                                                                       
       FMAX  = DXP(IC)                                                  
C                                                                       
       IX    = IC-1                                                     
                                                                        
       KG(NWILD) = IX/MA(NWILD) + 1                                     
       IF( NWILD .GT. 1 ) THEN                                          
           DO 400 J = 1,NWILD-1                                         
              NUM   = MOD(IX,MA(J+1))                                   
              KG(J) = NUM/MA(J) + 1                                     
  400      CONTINUE                                                     
       ENDIF                                                            
C                                                                       
C  ------------------------------------------------------------------   
C                     Sample and test a event                           
C  ------------------------------------------------------------------   
C                                                                       
      DO 600 NTRY = 1,MXTRY                                             
        WGT   = XJAC                                                    
        DO 550 J=1,NDIM                                                 
          IF( J .LE. NWILD) THEN                                        
             XN    = (KG(J)-DRN(IDUMY))*DXG+ONE                          
          ELSE                                                          
             XN    = ND*DRN(IDUMY) + ONE                                 
          ENDIF                                                         
          IAJ   = XN                                                    
          IF(IAJ .EQ. 1) THEN                                           
            XO    = XI(IAJ,J)                                           
            RC    = (XN-IAJ)*XO                                         
          ELSE                                                          
            XO    = XI(IAJ,J)-XI(IAJ-1,J)                               
            RC    = XI(IAJ-1,J)+(XN-IAJ)*XO                             
          ENDIF                                                         
          Y(J)  = XL(J) + RC*DX(J)                                      
          WGT   = WGT*XO*XND                                            
  550   CONTINUE                                                        
C                                                                       
*       FX    = F(Y)*WGT                                                
        FF    = F(Y)                                                    
        FX    = FF*WGT                                                  
        FUNCT = FX/FMAX                                                 
C                                                                       
        IF( FX .GT. 0.0D0 ) THEN                                        
*           IF( DRN(IDUMY) .LE. FUNCT ) GO TO 700                        
            XJ = DRN(IDUMY)                                              
            IF( XJ .LE. FUNCT ) GO TO 700                               
*           IF( XJ .LE. FUNCT ) THEN                                    
*               WRITE(6,9999) NTRY,IC,FF,WGT,XJ,FUNCT                   
*9999           FORMAT(1X,'NTRY,IC,FF,WGT,XJ,FUNCT = ',2I5,4E12.4)      
*               GO TO 700                                               
*           ENDIF                                                       
        ELSE                                                            
     .  IF( FX .LT. 0.0D0 ) THEN                                        
            WRITE(6,9100) IC                                            
 9100       FORMAT(                                                     
     .      /5X,'********** FATAL ERROR IN SPRING **********',          
     .      /5X,'* A negative value of function was found  *',          
     .      /5X,'*        in the ',I6,'-th Hypercube.      *',          
     .      /5X,'*******************************************')          
            WRITE(6,9405)                                               
 9405       FORMAT(5X,'------',3('+---------------'),'+')               
            WRITE(6,9410)                                               
 9410       FORMAT(5X,'    i       XL(i)             X       ',         
     .                '     XU(i)')                                     
            WRITE(6,9405)                                               
            DO 450 I = 1,NDIM                                           
                WRITE(6,9420) I,XL(I),Y(I),XU(I)                        
 9420           FORMAT(5X,I5,1P,3('  ',E14.6))                          
  450       CONTINUE                                                    
            WRITE(6,9405)                                               
            STOP                                                        
        ENDIF                                                           
C                                                                       
        CALL SHCLER                                                     
C                                                                       
  600 CONTINUE                                                          
                                                                        
      NTRY  = MXTRY + 1                                                 
                                                                        
  700 RETURN                                                            
      END                                                               

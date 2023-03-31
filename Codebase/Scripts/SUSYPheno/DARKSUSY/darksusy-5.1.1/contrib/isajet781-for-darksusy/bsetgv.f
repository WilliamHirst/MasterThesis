C***********************************************************************
C*                                                                     *
C*========================                                             *
C*    SUBROUTINE BSETGV( IFLAG )                                       *
C*========================                                             *
C*((Function))                                                         *
C*    Refine the grid sizes                                            *
C*                                                                     *
C*    Coded   by S.Kawabata    Aug. 1984 at Nagoya Univ.               *
C*    Last update              Oct. 1985 at KEK                        *
C*                                                                     *
C***********************************************************************
C                                                                       
      SUBROUTINE BSETGV( IFLAG )                                        
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)                   
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                    
     .               IG(MXDIM),NCALL                                    
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),      
     .               ND,NG,NPG,MA(MXDIM)                                
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT                          
      COMMON /BASE6/ D(NDMX,MXDIM),                                     
     .               ALPH,XSAVE(NDMX,MXDIM),XTI,XTSI,XACC,ITSX          
      REAL*4 STIME                                                      
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF                        
*                                                                       
                                                                        
      DIMENSION  XIN(NDMX),R(NDMX),DT(MXDIM),DDX(NDMX)                  
      DATA  ONE/1.0D0/,ZERO/0.0D0/,N0/0/,N1/1/                          
*                                                                       
*========= Save the grid information for the best accuracy ===========  
*                                                                       
      IF( ITSX .GT. 0 ) THEN                                            
          IF( IFLAG .EQ. 0 ) THEN                                       
              IF( IT .GE. 5 ) THEN                                      
                  IF( ( TI .GT. AVGI+SD) .AND. TSI .LT. XTSI ) THEN     
                      DO 400 J = 1, NDIM                                
                      DO 400 I = 1, ND                                  
                         XSAVE(I,J) = XI(I,J)                           
  400                 CONTINUE                                          
                      XACC         = TACC                               
                      ITSX         = IT                                 
                      XTI          = TI                                 
                      XTSI         = TSI                                
                  ENDIF                                                 
              ENDIF                                                     
          ELSE                                                          
              IF( ( XTI .GT. TI) .AND. XTSI .LT. TSI ) THEN             
                  DO 500 J = 1, NDIM                                    
                  DO 500 I = 1, ND                                      
                     XI(I,J) = XSAVE(I,J)                               
  500             CONTINUE                                              
*                ==========                                             
                   RETURN                                               
*                ==========                                             
              ENDIF                                                     
          ENDIF                                                         
      ENDIF                                                             
                                                                        
C======= SMOOTHING THE FUNCTION D(I,J)                                  
C                                                                       
        CLOGE   = 1.0D0/LOG(10.0D0)                                     
                                                                        
        NDM     = ND-1                                                  
        DO 780 J= N1,NDIM                                               
         IF( IG(J) .EQ. 1 ) THEN                                        
          DDX(1)= 0.5D0*(D(1,J) + D(2,J))                               
          DO 710 I=2,NDM                                                
            DDX(I)= (D(I+1,J) + D(I,J) + D(I-1,J))/3.D0                 
  710     CONTINUE                                                      
          DDX(ND) = 0.5D0*(D(NDM,J) + D(ND,J))                          
          DT(J) = 0.D0                                                  
          DO 720 I = 1, ND                                              
             D(I,J) = DDX(I)                                            
             DT(J)  = DT(J)+D(I,J)                                      
  720     CONTINUE                                                      
C                                                                       
C=========== REDEFINE THE GRID                                          
C                                                                       
                                                                        
          DTLOG   = LOG(DT(J))                                          
          DT10    = CLOGE*DTLOG                                         
          RC    = ZERO                                                  
          DO 730 I= N1,ND                                               
            R(I)  = ZERO                                                
            IF(D(I,J) .GT. ZERO) THEN                                   
               DILOG = LOG(D(I,J))                                      
               IF( DT10 - CLOGE*DILOG  .LE. 70.0D0 ) THEN               
                   XO    = DT(J)/D(I,J)                                 
                   R(I)  = ((XO-ONE)/(XO*(DTLOG-DILOG)))**ALPH          
               ELSE                                                     
C                  XO    = DT(J)/D(I,J)                                 
                   R(I)  = (DTLOG-DILOG)**(-ALPH)                       
               ENDIF                                                    
            ENDIF                                                       
            RC    = RC+R(I)                                             
  730     CONTINUE                                                      
          RC    = RC/ND                                                 
          K     = N0                                                    
          XN    = N0                                                    
          DR    = XN                                                    
          I     = K                                                     
  740  K     = K + N1                                                   
          DR    = DR+R(K)                                               
          XO    = XN                                                    
          XN    = XI(K,J)                                               
  750 IF(RC.GT.DR)GO TO 740                                             
          I     = I + N1                                                
          DR    = DR-RC                                                 
          XIN(I)= XN-(XN-XO)*DR/R(K)                                    
      IF(I.LT.NDM)GO TO 750                                             
          DO 760 I= N1,NDM                                              
            XI(I,J)= XIN(I)                                             
  760     CONTINUE                                                      
          XI(ND,J)= ONE                                                 
         ENDIF                                                          
  780   CONTINUE                                                        
C                                                                       
      RETURN                                                            
      END                                                               

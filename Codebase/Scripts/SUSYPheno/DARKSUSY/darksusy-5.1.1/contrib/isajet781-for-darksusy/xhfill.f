************************************************************************
*    ==============================                                    *
      SUBROUTINE XHFILL(ID, DX, FX )                                    
*    ==============================                                    *
* ((Function))                                                         *
*     To fill histograms.                                              *
*   This routine identifies the bin number which is to be updated      *
*   with weight FX*WGT.  Up to five points per histogram are able      *
*   to be stacked before calling BHUPDT or SHUPDT.                     *
* ((Input))                                                            *
*   ID    : Histogram identification number                            *
*   DX    : Input value                                                *
*   FX    : Input value of the function                                *
* ((Author))                                                           *
*   S.Kawabata         June '90 at KEK                                 *
*                                                                      *
************************************************************************
                                                                        
      REAL*8 DX,FX                                                      
      COMMON /BASE0/ IFLAG,IBASES                                       
      REAL*8         SCALLS,WGT,TI,TSI,TACC                             
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT                          
                                                                        
      PARAMETER ( NHS = 50, NSC = 50 )                                  
      INTEGER*4 XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD                
      COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),        
     .              NHIST, MAPL(4,NHS),                                 
     .              NSCAT, MAPD(4,NSC),                                 
     .              NW                                                  
      COMMON /PLOTB/ IBUF( 281*NHS + 2527*NSC )                         
      REAL*4         BUFF( 281*NHS + 2527*NSC )                         
      EQUIVALENCE (IBUF(1),BUFF(1))                                     
                                                                        
C     COMMON/PLOTLU/ LU                                                 
C                                                                       
      IF( NHIST .GT. 0 ) THEN                                           
C                                                                       
          I  = IABS(MOD( ID, 13 )) + 1                                  
          IF( XHASH(1, I) .EQ. 1 ) THEN                                 
            IF( ID .EQ. MAPL( 1, XHASH(2,I))) THEN                      
                IHIST = XHASH(2,I)                                      
                GO TO 200                                               
            ENDIF                                                       
          ELSEIF( XHASH(1, I) .GT. 1 ) THEN                             
            DO 100 K = 2, XHASH(1,I)+1                                  
               IF( ID .EQ. MAPL( 1, XHASH(K,I))) THEN                   
                   IHIST = XHASH(K,I)                                   
                   GO TO 200                                            
               ENDIF                                                    
  100       CONTINUE                                                    
          ENDIF                                                         
      ENDIF                                                             
C     IF( LU .GT. 0 ) THEN                                              
C         WRITE(LU,9000) ID                                             
C     ENDIF                                                             
C9000 FORMAT(1X,'No Histogram corresponds to ID =',I5,                  
C    .      /1X,' This call is neglected.')                             
      RETURN                                                            
C                                                                       
                                                                        
  200 X     = DX*1.0                                                    
                                                                        
          IX    = -1                                                    
          IP1   = MAPL(2,IHIST)                                         
          XMIN  = BUFF(IP1)                                             
          XMAX  = BUFF(IP1+1)                                           
          NXBIN = IBUF(IP1+2)                                           
          DEV   = BUFF(IP1+3)                                           
          IF(     X .LT. XMIN ) THEN                                    
                  IX   = 0                                              
          ELSEIF( X .GT. XMAX ) THEN                                    
                 IX   = NXBIN + 1                                       
          ELSE                                                          
                 IX   = INT((X - XMIN)/DEV + 1.0)                       
                 IF( IX .GT. NXBIN ) IX = NXBIN                         
          ENDIF                                                         
C        PRINT*,'ID, IHIST, IFBASE =',ID,IHIST,(IFBASE(I),I=1,NHIST)    
                                                                        
      IF( IBASES .EQ. 1 ) THEN                                          
                                                                        
          IP2       = MAPL(3,IHIST) + IX                                
          IBUF(IP2) = IBUF(IP2) + 1                                     
          FXWGT     = FX*WGT                                            
          IP2       = IP2 + 52                                          
          BUFF(IP2) = BUFF(IP2) + FXWGT                                 
          IP2       = IP2 + 52                                          
          BUFF(IP2) = BUFF(IP2) + FXWGT*FXWGT                           
*   Add March 1994                                                      
          IFBASE(IHIST) = 1                                             
                                                                        
      ELSE                                                              
C        PRINT*,'ID, IHIST, IFBASE =',ID,IHIST,(IFBASE(I),I=1,NHIST)    
                                                                        
         IP3        =  MAPL(4,IHIST)                                    
         IBUF(IP3)  = IX                                                
                                                                        
      ENDIF                                                             
                                                                        
C                                                                       
      RETURN                                                            
      END                                                               

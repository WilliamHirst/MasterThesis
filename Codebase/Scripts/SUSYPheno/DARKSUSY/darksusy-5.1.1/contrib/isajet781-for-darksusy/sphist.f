************************************************************************
*    =========================                                         *
      SUBROUTINE SPHIST( LU )                                           
*    =========================                                         *
* ((Purpose))                                                          *
*      To print the histogram for event generation                     *
* ((Input))                                                            *
*      LU   : logical unit number for the printer to be printed        *
*                                                                      *
* ((Author))                                                           *
*       S.Kawabata    April 1994                                       *
*                                                                      *
************************************************************************
                                                                        
      PARAMETER ( MXBIN = 51 )                                          
      COMMON/PLOTSP/ NBIN,IBUFSP( MXBIN )                               
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
                                                                        
      REAL  VAL(MXBIN),VLOG(MXBIN)                                      
      CHARACTER*50 CHARR,CHAR1                                          
      CHARACTER*52 SCALE                                                
      CHARACTER*1  BLNK,STAR,OO,AI,CN                                   
      DATA  YMAX / 50/                                                  
      DATA  BLNK /' '/, STAR /'*'/, OO /'O'/, AI /'I'/                  
                                                                        
      CN    = CHAR(12)                                                  
      IF( IPNT .EQ. 0 ) THEN                                            
          WRITE(LU,9000)                                                
 9000     FORMAT(/1H1,/1H )                                             
      ELSE                                                              
          WRITE(LU,9005) CN                                             
 9005     FORMAT(A1)                                                    
      ENDIF                                                             
          WRITE(LU,9102)                                                
 9102     FORMAT(5X,                                                    
     .  '************* Number of trials to get an event *************') 
                                                                        
      XMIN  = 1.0                                                       
      XMAX  = NBIN                                                      
      DEV   = 1.0                                                       
      NBIN1 = NBIN + 1                                                  
                                                                        
          NTOTAL     = IBUFSP(NBIN1)                                    
          VAL(NBIN1) = FLOAT(NTOTAL)                                    
          VMIN       = 0.0                                              
          VMAX       = VMIN                                             
          DO  55 I   = 1,NBIN                                           
              NTR    = IBUFSP(I)                                        
              VLS    = FLOAT( NTR )                                     
              NTOTAL = NTR + NTOTAL                                     
              IF( VMAX .LT. VLS ) VMAX = VLS                            
              VAL(I) = VLS                                              
   55     CONTINUE                                                      
                                                                        
          VLMAX = LOG10(VMAX)                                           
          VLMIN = VLMAX                                                 
                                                                        
           DO  60 I = 1,NBIN1                                           
               IF( VAL(I) .GT. 0.0 ) THEN                               
                   VLS   = LOG10( VAL(I) )                              
                   IF( I .LE. NBIN ) THEN                               
                       IF( VLS .LT. VLMIN ) VLMIN = VLS                 
                   ENDIF                                                
                   VLOG(I)  = VLS                                       
               ENDIF                                                    
   60      CONTINUE                                                     
                                                                        
           IF( VLMIN .LT. 0.0) THEN                                     
               VXMIN = IFIX(VLMIN) - 1.0                                
           ELSE                                                         
               VXMIN = IFIX(VLMIN)                                      
           ENDIF                                                        
           VXMAX = VLMAX                                                
           IFLG  = 1                                                    
           CALL XHRNGE( IFLG, VXMIN, VXMAX, VLMIN, VLMAX, VLSTP )       
           UNITL = (VLMAX-VLMIN)/YMAX                                   
                                                                        
       IFLG   = 0                                                       
           IF( VMIN .GE. 0.0 ) THEN                                     
               VXMAX  = 1.2*VMAX                                        
               VXMIN  = 0.0                                             
               CALL XHRNGE( IFLG, VXMIN, VXMAX, VMIN, VMAX, VSTP )      
           ELSE                                                         
               VXMAX  = 1.1*VMAX                                        
               VXMIN  = 1.1*VMIN                                        
               CALL XHRNGE( IFLG, VXMIN, VXMAX, VMIN, VMAX, VSTP )      
           ENDIF                                                        
                                                                        
       UNIT  = (VMAX-VMIN)/YMAX                                         
                                                                        
       CALL XHSCLE( IFLG, VMIN, VMAX, VSTP, UNIT, SCALE, CHAR1 )        
                                                                        
             WRITE(LU,9210) NTOTAL                                      
 9210        FORMAT(1X,'Total =',I10,' events',                         
     .        3X,'"*" : No. of events in Linear scale.')                
             WRITE(LU,9205) SCALE                                       
 9205        FORMAT(1X,'   x      Lg(dN/dx)  dN/dx',A52)                
             WRITE(LU,9251) CHAR1                                       
 9251        FORMAT(1X,                                                 
     .             '+-------+----------+-------+',                      
     .       A50 )                                                      
                                                                        
       VX    = ABS(XMAX)                                                
       XM    = ABS(XMIN)                                                
       IF( XM .GT. VX ) VX = XM                                         
                                                                        
       CALL XHORDR( VX, F2, ORD, IORD )                                 
                                                                        
       IF( VMIN .LT. 0.0 ) THEN                                         
           V1    = VMIN                                                 
           NUMBL = 1                                                    
           DO 150 I = 1, 80                                             
              V2    = V1 + UNIT                                         
              IF( V1 .LE. 0.0 .AND. V2 .GE. 0.0 ) THEN                  
                  NUMBL  = I                                            
                  GO TO 180                                             
              ENDIF                                                     
              V1    = V2                                                
  150      CONTINUE                                                     
       ENDIF                                                            
                                                                        
  180  DO 300 I = 1,NBIN1                                               
          VX   = VAL(I)                                                 
          IF( VMIN .GE. 0.0 ) THEN                                      
              IF( VX .GT. 0.0 ) THEN                                    
                  NUMBL  = (VLOG(I) - VLMIN)/UNITL + 1.0                
                  NUMB   = VX/UNIT + 1.0                                
              ELSE                                                      
                  NUMBL  = 0                                            
                  NUMB   = 0                                            
              ENDIF                                                     
              IF( NUMB .GT. 50 ) NUMB = 50                              
              IF( NUMBL.GT. 50 ) NUMBL= 50                              
              DO 200 K = 1,50                                           
                 IF(     ( K .GT. NUMBL) .AND. (K .GT. NUMB ) ) THEN    
                           IF( K .EQ. 50 ) THEN                         
                               CHARR(K:K) = AI                          
                           ELSE                                         
                               CHARR(K:K) = BLNK                        
                           ENDIF                                        
                 ELSEIF( ( K .LE. NUMBL) .AND. (K .GT. NUMB )) THEN     
                           CHARR(K:K) = OO                              
                 ELSEIF( ( K .GT. NUMBL) .AND. (K .LE. NUMB )) THEN     
                           CHARR(K:K) = STAR                            
                 ELSEIF( ( K .LE. NUMBL) .AND. (K .LE. NUMB)) THEN      
                           IF( NUMB .GE. NUMBL ) THEN                   
                               CHARR(K:K) = OO                          
                           ELSE                                         
                               CHARR(K:K) = STAR                        
                           ENDIF                                        
                 ENDIF                                                  
  200         CONTINUE                                                  
          ELSE                                                          
                                                                        
              V1          = VMIN                                        
              NHIG        = 1                                           
              DO 220  J = 1, 50                                         
                 V2     = V1 + UNIT                                     
                 IF( VX .GE. V1 .AND. VX .LT. V2 ) THEN                 
                     NHIG   = J                                         
                     GO TO 240                                          
                 ENDIF                                                  
                 V1    = V2                                             
  220         CONTINUE                                                  
  240         NLOW   = NUMBL                                            
              IF( NHIG .LT. NLOW) THEN                                  
                  NX    = NHIG                                          
                  NHIG  = NLOW                                          
                  NLOW  = NX                                            
              ENDIF                                                     
                                                                        
              DO 250 K = 1, 49                                          
                 IF(     K .EQ. NUMBL ) THEN                            
                         CHARR(K:K) = AI                                
                 ELSEIF( K .GT. NHIG ) THEN                             
                         CHARR(K:K) = BLNK                              
                 ELSEIF( K .LT. NLOW ) THEN                             
                         CHARR(K:K) = BLNK                              
                 ELSE                                                   
                     IF( K .EQ. NHIG .AND. K .EQ. NLOW) THEN            
                         CHARR(K:K) = AI                                
                     ELSE                                               
                         CHARR(K:K) = STAR                              
                     ENDIF                                              
                 ENDIF                                                  
  250         CONTINUE                                                  
              CHARR(50:50) = AI                                         
          ENDIF                                                         
                                                                        
             NX  = VAL(I)                                               
             VX     = VAL(I)                                            
             VX1    = VX                                                
             IF( VX .LT. 0.0 ) VX1 = -VX                                
             CALL XHORDR( VX1, F2, ORDER, IORDR )                       
             F2     = VX/ORDER                                          
             IF( I .EQ. NBIN1 ) THEN                                    
                 WRITE(LU,9400) IORD,F2,IORDR,NX,CHARR                  
 9400            FORMAT(1X,'I  E',I3,' I',F6.3,'E',I3,'I',              
     .                                            I7,'I',A50)           
             ELSE                                                       
                   XM  = (XMIN + DEV*(I - 1))/ORD                       
                   WRITE(LU,9440) XM,F2,IORDR,NX,CHARR                  
 9440              FORMAT(1X,'I',F6.3,' I',F6.3,'E',I3,'I',             
     .                                        I7,'I',A50)               
             ENDIF                                                      
                                                                        
  300  CONTINUE                                                         
                                                                        
       IF( VMIN .GE. 0.0 ) THEN                                         
           CALL XHSCLE( 1, VLMIN, VLMAX, VLSTP, UNITL, SCALE, CHAR1)    
           VXMIN  = 10**VLMIN                                           
       ENDIF                                                            
                                                                        
           WRITE(LU,9251) CHAR1                                         
           WRITE(LU,9205) SCALE                                         
           WRITE(LU,9360)                                               
 9360      FORMAT(30X,'"O" : No. of Events in Log. scale.')             
                                                                        
C                                                                       
                                                                        
      RETURN                                                            
      END                                                               

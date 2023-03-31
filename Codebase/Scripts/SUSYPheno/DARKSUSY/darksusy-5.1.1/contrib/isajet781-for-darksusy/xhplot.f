************************************************************************
*    =====================================                             *
      SUBROUTINE XHPLOT( LU, IFG, IHIST )                               
*    =====================================                             *
* ((Purpose))                                                          *
*      To print histograms for BASES and SPRING.                       *
* ((Input))                                                            *
*      IFG  : Flag which indicats whether this is called by BASES      *
*             or SPRING.  IFG = ( 0 / anyother) = ( By BASES/ SPRING)  *
*      IHIST: Serial number of the histogram                           *
*                                                                      *
* ((Author))                                                           *
*       S.Kawabata    June '90 at KEK                                  *
*     Last update     March '94                                        *
*                                                                      *
************************************************************************
                                                                        
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
                                                                        
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
                                                                        
      REAL  VAL(0:51),VLOG(0:51),VERR(0:51)                             
      CHARACTER*50 CHARR,CHAR1                                          
      CHARACTER*52 SCALE                                                
      CHARACTER*1  BLNK,STAR,OO,AI,CN                                   
      DATA  YMAX / 50/                                                  
      DATA  BLNK /' '/, STAR /'*'/, OO /'O'/, AI /'I'/                  
                                                                        
      CN    = CHAR(12)                                                  
      IP3   = MAPL(4,IHIST)                                             
      IF(     IFG .EQ. 0 ) THEN                                         
            IF( IPNT .EQ. 0 ) THEN                                      
                WRITE(LU,9000)                                          
            ELSE                                                        
                WRITE(LU,9005) CN                                       
            ENDIF                                                       
 9000       FORMAT(/1H1)                                                
 9005       FORMAT(A1)                                                  
            WRITE(LU,9100) MAPL(1,IHIST),(BUFF(I),I=IP3+1,IP3+16)       
 9100       FORMAT(1X,'Histogram (ID =',I3,' ) for ',16A4)              
      ELSEIF( IFG .EQ. -10 ) THEN                                       
            IF( IPNT .EQ. 0 ) THEN                                      
                WRITE(LU,9000)                                          
            ELSE                                                        
                WRITE(LU,9005) CN                                       
            ENDIF                                                       
            WRITE(LU,9102) (BUFF(I),I=IP3+1,IP3+16)                     
 9102       FORMAT(5X,16A4)                                             
      ELSE                                                              
            IF( IPNT .EQ. 0 ) THEN                                      
                WRITE(LU,9000)                                          
            ELSE                                                        
                WRITE(LU,9005) CN                                       
            ENDIF                                                       
            WRITE(LU,9105) MAPL(1,IHIST),(BUFF(I),I=IP3+1,IP3+16)       
 9105       FORMAT(                                                     
     .      1X,'Additional Histogram (ID =',I3,' ) for ',16A4)          
      ENDIF                                                             
                                                                        
      IP1   = MAPL(2,IHIST)                                             
      XMIN  = BUFF(IP1)                                                 
      XMAX  = BUFF(IP1+1)                                               
      NXBIN = IBUF(IP1+2) + 1                                           
      DEV   = BUFF(IP1+3)                                               
      IP2   = MAPL(3,IHIST)                                             
                                                                        
      IF( IFG .EQ. 0 ) THEN                                             
C         NTOTAL     = SCALLS                                           
          FACT       = 1./(SCALLS*DEV)                                  
          IPF   = IP2 + 156                                             
          IPF2  = IPF + 52                                              
          VAL(0)     = BUFF(IPF)/SCALLS                                 
          VAL(NXBIN) = BUFF(IPF+NXBIN)/SCALLS                           
          VMAX       = FACT*BUFF(IPF+1)                                 
          VMIN       = VMAX                                             
          DO  50 I   = 1,NXBIN-1                                        
              TX     = BUFF(I+IPF)                                      
              NX     = IBUF(I+IP2)                                      
              VLS    = TX*FACT                                          
              IF( VMAX .LT. VLS ) VMAX = VLS                            
              IF( VMIN .GT. VLS ) VMIN = VLS                            
              VAL(I) = VLS                                              
              IF( NX .GT. 1 ) THEN                                      
                  DEV2   =  NX*BUFF(I+IPF2)-TX*TX                       
                  IF( DEV2 .LE. 0.0 ) THEN                              
                      VERR(I)= 0.0                                      
                  ELSE                                                  
                      VERR(I)= FACT*SQRT( DEV2/( NX-1 ))                
                  ENDIF                                                 
              ELSEIF( NX .EQ. 1 ) THEN                                  
                  VERR(I)= VLS                                          
              ELSE                                                      
                  VERR(I)= 0.0                                          
              ENDIF                                                     
                                                                        
   50     CONTINUE                                                      
      ELSE                                                              
          IPX   = IP2 + 52                                              
          VAL(0)     = BUFF(IPX)                                        
          VAL(NXBIN) = BUFF(IPX+NXBIN)                                  
          NTOTAL     = INT(VAL(0)) + INT(VAL(NXBIN))                    
          VMIN       = 0.0                                              
          VMAX       = VMIN                                             
          DO  55 I   = 1,NXBIN-1                                        
              VLS    = BUFF(I+IPX)                                      
              NTOTAL = INT(VLS) + NTOTAL                                
              IF( VMAX .LT. VLS ) VMAX = VLS                            
              VAL(I) = VLS                                              
              IF( VLS .GT. 0.0 ) THEN                                   
                  VERR(I) = SQRT(VLS)                                   
              ELSE                                                      
                  VERR(I) = 0.0                                         
              ENDIF                                                     
   55     CONTINUE                                                      
       ENDIF                                                            
***                                                                     
       IF( VMAX .EQ. 0.0 .AND. VMIN .EQ. 0.0) THEN                      
           V0 = VAL(0)                                                  
           VM = VAL(NXBIN)                                              
           IF( V0 .GE. 0.0 .AND. VM .GE. 0.0 ) THEN                     
               VMIN  = 0.0                                              
               IF( V0 .GT. VM  ) THEN                                   
                   VMAX = V0                                            
               ELSE                                                     
                   VMAX = VM                                            
               ENDIF                                                    
           ELSEIF( V0 .LT. 0.0 .AND. VM .LT. 0.0 ) THEN                 
               VMAX  = 0.0                                              
               IF( V0 .LT. VM ) THEN                                    
                   VMIN  = V0                                           
               ELSE                                                     
                   VMIN  = VM                                           
               ENDIF                                                    
           ELSEIF( V0 .GT. VM ) THEN                                    
               VMAX  = V0                                               
               VMIN  = VM                                               
           ELSE                                                         
               VMAX  = VM                                               
               VMIN  = V0                                               
           ENDIF                                                        
       ENDIF                                                            
***                                                                     
       IF( VMIN .GE. 0.0 ) THEN                                         
C//VV                                                                   
           IF( VMAX .GT. 0.0 ) THEN                                     
               VLMAX = LOG10(VMAX)                                      
           ELSE                                                         
               VLMAX = 2.0                                              
           ENDIF                                                        
C//                                                                     
           VLMIN = VLMAX                                                
           DO  60 I = 0,NXBIN                                           
               IF( VAL(I) .GT. 0.0 ) THEN                               
                   VLS   = LOG10( VAL(I) )                              
                   IF( I .GT. 0 .AND. I .LT. NXBIN ) THEN               
                       IF( VLS .LT. VLMIN ) VLMIN = VLS                 
                   ENDIF                                                
                   VLOG(I)  = VLS                                       
C//VV                                                                   
C              ELSE                                                     
C                  VLOG(I)  = 0.0                                       
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
                                                                        
       ENDIF                                                            
                                                                        
       IFLG   = 0                                                       
       IF( VMAX .GT. 0.0 ) THEN                                         
           IF( VMIN .GE. 0.0 ) THEN                                     
               VXMAX  = 1.2*VMAX                                        
               VXMIN  = 0.0                                             
               CALL XHRNGE( IFLG, VXMIN, VXMAX, VMIN, VMAX, VSTP )      
           ELSE                                                         
               VXMAX  = 1.1*VMAX                                        
               VXMIN  = 1.1*VMIN                                        
               CALL XHRNGE( IFLG, VXMIN, VXMAX, VMIN, VMAX, VSTP )      
           ENDIF                                                        
       ELSE                                                             
          VXMAX  = 0.0                                                  
          VXMIN  = 1.1*VMIN                                             
          CALL XHRNGE( IFLG, VXMIN, VXMAX, VMIN, VMAX, VSTP )           
       ENDIF                                                            
                                                                        
       UNIT  = (VMAX-VMIN)/YMAX                                         
                                                                        
       CALL XHSCLE( IFLG, VMIN, VMAX, VSTP, UNIT, SCALE, CHAR1 )        
C                                                                       
C                                                                       
       IF( IFG .EQ. 0 ) THEN                                            
           WRITE(LU,9150)                                               
 9150      FORMAT(30X,'Linear Scale indicated by "*"')                  
           WRITE(LU,9200) SCALE                                         
 9200      FORMAT(1X,'    x      d(Sigma)/dx    ',A52)                  
           WRITE(LU,9250) CHAR1                                         
 9250      FORMAT(1X,                                                   
     .                '+-------+------------------+',                   
     .           A50 )                                                  
       ELSE                                                             
             WRITE(LU,9210) NTOTAL                                      
 9210        FORMAT(1X,'Total =',I10,' events',                         
     .        3X,'"*" : No. of events in Linear scale.')                
             WRITE(LU,9205) SCALE                                       
 9205        FORMAT(1X,'   x      Lg(dN/dx)  dN/dx',A52)                
             WRITE(LU,9251) CHAR1                                       
 9251        FORMAT(1X,                                                 
     .             '+-------+----------+-------+',                      
     .       A50 )                                                      
       ENDIF                                                            
                                                                        
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
                                                                        
  180  DO 300 I = 0,NXBIN                                               
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
                                                                        
          IF( IFG .EQ. 0 ) THEN                                         
                                                                        
              NX     = IBUF(I+IP2)                                      
              VX     = VAL(I)                                           
              VX1    = VX                                               
              IF( VX .LT. 0.0 ) VX1 = -VX                               
                                                                        
                                                                        
              IF( I .EQ. 0 .OR. I. EQ. NXBIN ) THEN                     
                  CALL XHORDR( VX1, F2, ORDER, IORDR )                  
                  F2     = VX/ORDER                                     
                  WRITE(LU,9300) IORD,F2,IORDR,CHARR                    
 9300             FORMAT(1X,'I  E',I3,' I',F6.3,8X,'E',I3,              
     .                                     'I',A50)                     
              ELSE                                                      
                  XM    = (XMIN + DEV*(I-1))/ORD                        
                  VE     = VERR(I)                                      
                  IF( VE .GT. VX1 ) THEN                                
                      CALL XHORDR(  VE, F2, ORDER, IORDR )              
                  ELSE                                                  
                      CALL XHORDR( VX1, F2, ORDER, IORDR )              
                  ENDIF                                                 
                  F2   = VX/ORDER                                       
                  VE   = VE/ORDER                                       
                  WRITE(LU,9340) XM,F2,VE,IORDR,CHARR                   
 9340             FORMAT(1X,'I', F6.3,' I',F6.3,'+-',F5.3,' E',I3,      
     .                                    'I',A50)                      
             ENDIF                                                      
          ELSE                                                          
             NX  = VAL(I)                                               
             VX     = VAL(I)                                            
             VX1    = VX                                                
             IF( VX .LT. 0.0 ) VX1 = -VX                                
             CALL XHORDR( VX1, F2, ORDER, IORDR )                       
             F2     = VX/ORDER                                          
             IF( I .EQ. 0 .OR. I .EQ. NXBIN ) THEN                      
                 WRITE(LU,9400) IORD,F2,IORDR,NX,CHARR                  
 9400            FORMAT(1X,'I  E',I3,' I',F6.3,'E',I3,'I',              
     .                                            I7,'I',A50)           
             ELSE                                                       
                   XM  = (XMIN + DEV*(I - 1))/ORD                       
                   WRITE(LU,9440) XM,F2,IORDR,NX,CHARR                  
 9440              FORMAT(1X,'I',F6.3,' I',F6.3,'E',I3,'I',             
     .                                        I7,'I',A50)               
             ENDIF                                                      
          ENDIF                                                         
  300  CONTINUE                                                         
                                                                        
       IF( VMIN .GE. 0.0 ) THEN                                         
           CALL XHSCLE( 1, VLMIN, VLMAX, VLSTP, UNITL, SCALE, CHAR1)    
           VXMIN  = 10**VLMIN                                           
       ENDIF                                                            
                                                                        
       IF( IFG .EQ. 0 ) THEN                                            
           WRITE(LU,9250) CHAR1                                         
           IF( VMIN .GE. 0.0 ) THEN                                     
               WRITE(LU,9200) SCALE                                     
               WRITE(LU,9260)                                           
 9260          FORMAT(30X,'Logarithmic Scale indicated by "O"')         
           ELSE                                                         
               WRITE(LU,9200) SCALE                                     
           ENDIF                                                        
       ELSE                                                             
           WRITE(LU,9251) CHAR1                                         
           WRITE(LU,9205) SCALE                                         
           WRITE(LU,9360)                                               
 9360      FORMAT(30X,'"O" : No. of Events in Log. scale.')             
       ENDIF                                                            
                                                                        
C                                                                       
  500  CONTINUE                                                         
                                                                        
      RETURN                                                            
      END                                                               

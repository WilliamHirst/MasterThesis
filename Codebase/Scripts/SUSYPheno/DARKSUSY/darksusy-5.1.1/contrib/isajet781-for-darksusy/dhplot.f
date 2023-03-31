************************************************************************
*     =========================                                        *
       SUBROUTINE DHPLOT( LU )                                          
*     =========================                                        *
* ((Purpose))                                                          *
*      To print scatter plots for BASES and SPRING                     *
*                                                                      *
* ((Author))                                                           *
*       S.Kawabata    June '90                                         *
*                                                                      *
************************************************************************
                                                                        
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
                                                                        
      CHARACTER*1  PLUS,MINUS,BLNK,STAR,NUM(0:9),NEG(0:9),SHARP,PNT     
      REAL*4       X(50)                                                
      CHARACTER*1 CHARR(50), CN                                         
      CHARACTER*80 FORM1,FORM                                           
      DATA  PLUS /'+'/, MINUS /'-'/, BLNK /' '/, STAR /'*'/             
      DATA  SHARP /'#'/,  PNT /'.'/                                     
      DATA  NUM  / '0','1','2','3','4','5','6','7','8','9'/             
      DATA  NEG  / '-','a','b','c','d','e','f','g','h','i'/             
*                                                                      *
*--------------------------- Entry point ------------------------------*
*                                                                      *
      CN   = CHAR(12)                                                   
                                                                        
      IF( NSCAT .GT. 0 ) THEN                                           
         DO 900 ISCAT = 1, NSCAT                                        
            IP3   = MAPD(4,ISCAT)                                       
            IF( IPNT .EQ. 0 ) THEN                                      
                WRITE(LU,9010)                                          
            ELSE                                                        
                WRITE(LU,9020) CN                                       
            ENDIF                                                       
 9010       FORMAT(/1H1)                                                
 9020       FORMAT(A1)                                                  
            WRITE(LU,9100) MAPD(1,ISCAT),(BUFF(I), I=IP3+2,IP3+17)      
 9100       FORMAT(/5X,'Scat_Plot (ID =',I3,' ) for ',16A4,/)           
                                                                        
            IP1   = MAPD(2,ISCAT)                                       
            XL    = BUFF(IP1)                                           
            XU    = BUFF(IP1+1)                                         
            NX    = IBUF(IP1+2)                                         
            DX    = BUFF(IP1+3)                                         
            XM    = ABS(XU)                                             
            XX    = ABS(XL)                                             
            IF( XX .GT. XM ) XM = XX                                    
            CALL XHORDR( XU, FX, XORD, IXORD)                           
            YL    = BUFF(IP1+4)                                         
            YU    = BUFF(IP1+5)                                         
            NY    = IBUF(IP1+6)                                         
            DY    = BUFF(IP1+7)                                         
            MIDY  = NY/2                                                
            IF( MIDY .EQ. 0 ) MIDY = 1                                  
            YM    = ABS(YU)                                             
            YY    = ABS(YL)                                             
            IF( YY .GT. YM ) YM = YY                                    
            CALL XHORDR( YM, FY, YORD, IYORD)                           
            IP2   = MAPD(3,ISCAT)                                       
            NTOTAL= IBUF(IP2)                                           
            VMAX  = BUFF(IP2+1)                                         
            VMIN  = VMAX                                                
            DO 100 J = 0,NY-1                                           
               IB    = NX*J + IP2                                       
               DO 100 I = 1,NX                                          
                  VLS    = BUFF( I + IB )                               
                  IF( VLS .GT. VMAX ) VMAX = VLS                        
                  IF( VLS .LT. VMIN ) VMIN = VLS                        
  100       CONTINUE                                                    
***                                                                     
            IF( VMAX .EQ. 0.0 .AND. VMIN .EQ. 0.0 ) THEN                
                VMAX  = 10.0                                            
                VMIN  = 0.0                                             
            ENDIF                                                       
***                                                                     
            IF( VMAX .GT. -VMIN ) THEN                                  
                UNIT = ABS(VMAX)/11.0                                   
            ELSE                                                        
                UNIT = ABS(VMIN)/11.0                                   
            ENDIF                                                       
            WRITE(FORM1,9200) NX                                        
*9200       FORMAT('(7X,''E'',I3,3X,''+'',',I2,'(''--''),''-+'')')      
 9200       FORMAT('(7X,''E'',I3,3X,''+'',',I2,'(''-''),''+'')')        
            WRITE(LU,FORM1) IYORD                                       
            DO 300 L = NY-1,0,-1                                        
               IB     = NX*L + IP2                                      
               DO 200 I = 1,NX                                          
                 XNUM   = BUFF( I + IB )/UNIT                           
                 IF( XNUM .LT. 0 0 ) THEN                               
                     NUMB   = XNUM - 1.0                                
                     IF(     NUMB .GE. -1 )THEN                         
                             CHARR(I) = MINUS                           
                     ELSEIF( NUMB .GE. -10 ) THEN                       
                            CHARR(I) = NEG(-NUMB-1)                     
                     ELSE                                               
                            CHARR(I) = SHARP                            
                     ENDIF                                              
                 ELSE                                                   
                     NUMB   = XNUM + 1.0                                
                     IF(     XNUM .EQ. 0.0 ) THEN                       
                             CHARR(I) = BLNK                            
                     ELSEIF( NUMB .LE.  1 ) THEN                        
                             CHARR(I) = PLUS                            
                             IF( VMIN .GE. 0.0 ) CHARR(I) = PNT         
                     ELSEIF( NUMB .LE. 10 ) THEN                        
                             CHARR(I) = NUM(NUMB-1)                     
                     ELSE                                               
                             CHARR(I) = STAR                            
                     ENDIF                                              
                 ENDIF                                                  
  200          CONTINUE                                                 
                                                                        
               Y   = (L*DY + YL)/YORD                                   
               IF( L .EQ. MIDY ) THEN                                   
                   WRITE(FORM,9300) NX                                  
*9300              FORMAT('(5X,F6.3,'' Y I'',',I2,'(1X,A1),'' I'')')    
 9300              FORMAT('(5X,F6.3,'' Y I'',',I2,'A1,''I'')')          
               ELSE                                                     
                   WRITE(FORM,9310) NX                                  
*9310              FORMAT('(5X,F6.3,''   I'',',I2,'(1X,A1),'' I'')')    
 9310              FORMAT('(5X,F6.3,''   I'',',I2,'A1,''I'')')          
               ENDIF                                                    
               WRITE(LU,FORM) Y,(CHARR(M),M=1,NX)                       
                                                                        
  300       CONTINUE                                                    
                                                                        
            WRITE(LU,FORM1) IYORD                                       
                                                                        
            NXH   = NX/2                                                
            IF( NXH .EQ. 0 ) NXH = 1                                    
            WRITE(FORM,9400) NXH                                        
                                                                        
*           WRITE(FORM,9400) NX                                         
 9400       FORMAT('(6X,''Low-'',5X,',I2,'X,''X'')')                    
            WRITE(LU,FORM)                                              
                                                                        
            XORD     = XORD*10.                                         
            DO 400 I = 1, NX                                            
               X(I)  = ((I-1)*DX + XL)/XORD                             
               IF( X(I) .LT. 0.0 ) THEN                                 
                   CHARR(I)  = MINUS                                    
                   X(I)      = -X(I)                                    
               ELSE                                                     
                   CHARR(I)  = BLNK                                     
               ENDIF                                                    
  400       CONTINUE                                                    
            WRITE(FORM1,9500) NX                                        
*9500       FORMAT('(6X,''Edge'',5X,',I2,'(1X,A1))')                    
 9500       FORMAT('(6X,''Edge'',5X,',I2,'A1)')                         
            WRITE(LU,FORM1) (CHARR(M),M=1,NX)                           
                                                                        
            XORD      = 1.0                                             
            DO 600 I  = 1,5                                             
               IF( I .EQ. 2 ) THEN                                      
                   WRITE(FORM,9602) NX                                  
 9602              FORMAT('(7X,''E'',I3,4X',I2,                         
     .                    '(''.''))')                                   
                   WRITE(LU,FORM) IXORD                                 
               ELSE                                                     
                   DO 500 J = 1, NX                                     
                      XX        = X(J)*10.0                             
                      NUMB      = XX                                    
                      CHARR(J)  = NUM(NUMB)                             
                      X(J)      = XX - FLOAT(NUMB)                      
  500              CONTINUE                                             
                   IF(     I .EQ. 4 ) THEN                              
                           WRITE(FORM,9604) NX                          
 9604                      FORMAT('(7X,''Low-'',4X,',I2,                
     .                            'A1)')                                
                   ELSEIF( I .EQ. 5 ) THEN                              
                           WRITE(FORM,9605) NX                          
 9605                      FORMAT('(7X,''Edge'',4X,',I2,                
     .                            'A1)')                                
                   ELSE                                                 
                           WRITE(FORM,9601) NX                          
 9601                      FORMAT('(15X,',I2,                           
     .                            'A1)')                                
                   ENDIF                                                
                   WRITE(LU,FORM) (CHARR(M),M=1,NX)                     
               ENDIF                                                    
  600       CONTINUE                                                    
                                                                        
  900    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               

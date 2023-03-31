************************************************************************
*    ============================================                      *
      SUBROUTINE XHINIT(ID,DXMIN,DXMAX,NBIN,TNAME)                      
*    ============================================                      *
* ((Function))                                                         *
*     To define a histogram.                                           *
* ((Input))                                                            *
*    ID   : Histogram identification number                            *
*    DXMIN: Lower limit of the histogram                               *
*    DXMAX: Upper limit of the histogram                               *
*    NBIN : Number of bins for the histogram (Max. is 50 )             *
*    TNAME: Title of the histogram in the character string (upto 64    *
*            characters)                                               *
* ((Author))                                                           *
*    S.Kawabata    June '90                                            *
*                                                                      *
************************************************************************
                                                                        
      REAL*8 DXMIN, DXMAX                                               
      CHARACTER*(*) TNAME                                               
      CHARACTER*68  NAME                                                
                                                                        
      PARAMETER ( NHS = 50, NSC = 50 )                                  
      INTEGER*4 XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD                
      COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),        
     .              NHIST, MAPL(4,NHS),                                 
     .              NSCAT, MAPD(4,NSC),                                 
     .              NW                                                  
      COMMON /PLOTB/ IBUF( 281*NHS + 2527*NSC )                         
      REAL*4         BUFF( 281*NHS + 2527*NSC )                         
      EQUIVALENCE (IBUF(1),BUFF(1))                                     
                                                                        
*     COMMON/XHCNTL/ LOCK                                               
      COMMON/PLOTLU/ LU                                                 
                                                                        
      IF( NHIST .GE. NHS ) THEN                                         
*         IF( LOCK .NE. 0 ) RETURN                                      
          IF( LU .GT. 0 ) THEN                                          
              WRITE(LU,9000) NHIST,ID                                   
          ENDIF                                                         
 9000     FORMAT(1X,'Number of Histograms exceeds ',I3,' at ID = ',I3,  
     .            /1X,'This call is neglected.')                        
          RETURN                                                        
      ENDIF                                                             
                                                                        
      IEXIST = 0                                                        
      I  = IABS(MOD( ID, 13 )) + 1                                      
      NH = XHASH(1,I)                                                   
                                                                        
      IF( NH .EQ. 1 ) THEN                                              
            IF( ID .EQ. MAPL( 1, XHASH(2,I))) THEN                      
*               IF( LOCK .NE. 0 ) RETURN                                
                IEXIST = XHASH(2,I)                                     
            ENDIF                                                       
      ELSEIF( NH .GT. 1 ) THEN                                          
          DO 100 K = 2, NH+1                                            
            IF( ID .EQ. MAPL( 1, XHASH(K,I))) THEN                      
*               IF( LOCK .NE. 0 ) RETURN                                
                IEXIST = XHASH(K,I)                                     
                GO TO 110                                               
            ENDIF                                                       
  100    CONTINUE                                                       
  110    CONTINUE                                                       
      ENDIF                                                             
      XMIN  = DXMIN*1.0                                                 
      XMAX  = DXMAX*1.0                                                 
                                                                        
      IF( IEXIST .GT. 0 ) THEN                                          
          IF( LU .GT. 0 ) THEN                                          
              WRITE(LU,9100) ID                                         
          ENDIF                                                         
 9100     FORMAT(1X,'Histogram ID (',I3,' ) exists already.')           
          IP1    =  MAPL(2,IEXIST)                                      
          IF(( XMIN .EQ. BUFF(IP1))   .AND.                             
     .       ( XMAX .EQ. BUFF(IP1+1)) .AND.                             
     .       ( NBIN .EQ. IBUF(IP1+2)) )    THEN                         
               IF( LU .GT. 0 ) THEN                                     
                   WRITE(LU,9110)                                       
               ENDIF                                                    
 9110          FORMAT(1X,' This call is neglected.')                    
               RETURN                                                   
          ENDIF                                                         
          IF( LU .GT. 0 ) THEN                                          
              WRITE(LU,9120) ID,XMIN,XMAX,NBIN                          
          ENDIF                                                         
 9120     FORMAT(1X,'Histogram ( ID =',I3,' ) parameters are replaced', 
     .          /1X,'by the following new parameters :',                
     .          /1X,' XMIN(',E12.5,')  XMAX(',E12.5,' )  NBIN(',I4,' )')
      ENDIF                                                             
                                                                        
      IF((NHIST .GE. NHS) .AND. (ID .GT. 0) ) THEN                      
*         IF( LOCK .NE. 0 ) RETURN                                      
          IF( LU .GT. 0 ) THEN                                          
              WRITE(LU,9000) NHS,ID                                     
          ENDIF                                                         
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF(NBIN  .GT. 50 ) THEN                                           
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9200) NBIN,ID                                     
         ENDIF                                                          
 9200    FORMAT(1X,'Bin size (',I3,' )  exceeds 50 at ID =',I5,         
     .         /1X,' This call is neglected.')                          
         RETURN                                                         
      ENDIF                                                             
      IF(XMIN  .GE. XMAX ) THEN                                         
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9300) ID                                          
         ENDIF                                                          
 9300    FORMAT(1X,'Lower limit is larger than upper at ID =',I5,       
     .         /1X,' This call is neglected.')                          
         RETURN                                                         
      ENDIF                                                             
      IF(XHASH(1,I) .GE. NHS) THEN                                      
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9400) I                                           
         ENDIF                                                          
 9400    FORMAT(1X,I5,'-th Hash table overflow',                        
     .         /1X,' This call is neglected.')                          
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF( IEXIST .GT. 0 ) THEN                                          
          NHST     = IEXIST                                             
      ELSE                                                              
          NHIST        = NHIST + 1                                      
          XHASH(1,I)   = XHASH(1,I) + 1                                 
          K            = XHASH(1,I) + 1                                 
          XHASH(K,I)   = NHIST                                          
          NHST         = NHIST                                          
          IP1    = NW + 1                                               
          NW  = NW + 281                                                
          MAPL(1,NHST)  = ID                                            
          MAPL(2,NHST)  = IP1                                           
      ENDIF                                                             
         BUFF(IP1     ) = XMIN                                          
         BUFF(IP1 +  1) = XMAX                                          
         IBUF(IP1 +  2) = NBIN                                          
         DEV            = XMAX - XMIN                                   
         BUFF(IP1 +  3) = DEV/NBIN                                      
      IP2   = IP1 + 4                                                   
         MAPL(3,NHST)  = IP2                                            
      IP3   = IP1 + 264                                                 
         MAPL(4,NHST)  = IP3                                            
         IBUF(IP3)     = -1                                             
                                                                        
         I1   = IP3 + 1                                                 
         I2   = I1 + 15                                                 
         NAME = TNAME                                                   
         READ(NAME,9800) (BUFF(I),I=I1,I2)                              
 9800    FORMAT(16A4)                                                   
                                                                        
C                                                                       
 1000 CONTINUE                                                          
      RETURN                                                            
      END                                                               

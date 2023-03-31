************************************************************************
*  =================================================================== *
      SUBROUTINE DHINIT(ID,DXMIN,DXMAX,NXBIN,DYMIN,DYMAX,NYBIN,TNAME)   
*  =================================================================== *
* ((Function))                                                         *
*     To define a scatter plot                                         *
* ((Input))                                                            *
*    ID   : scatter plot identification number                         *
*    DXMIN: Lower limit of X for the scatter plot                      *
*    DXMAX: Upper limit of X for the scatter plot                      *
*    NXBIN: Number of bins of X for the plot (Max. is 50 )             *
*    DYMIN: Lower limit of Y for the scatter plot                      *
*    DYMAX: Upper limit of Y for the scatter plot                      *
*    NYBIN: Number of bins of Y for the plot (Max. is 50 )             *
*    TNAME: Title of the plot in the character string (upto 64         *
*            characters)                                               *
* ((Author))                                                           *
*    S.Kawabata     June '90 at KEK                                    *
*                                                                      *
************************************************************************
                                                                        
      REAL*8 DXMIN,DXMAX,DYMIN,DYMAX                                    
      CHARACTER*(*) TNAME                                               
      CHARACTER*64 NAME                                                 
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
*                                                                      *
*--------------------------- Entry point ------------------------------*
*                                                                      *
*======================================================================*
*               Find the scatter plot ID in the table                  *
*======================================================================*
*                                                                      *
      IF( NSCAT .GE. NSC ) THEN                                         
*         IF( LOCK .NE. 0 ) RETURN                                      
          IF( LU .GT. 0 ) THEN                                          
            WRITE(LU,9000) NSCAT,ID                                     
 9000       FORMAT(1X,'Numberof Scat_plots exceeds ',I3,' at ID = ',I3, 
     .            /1X,'This call is neglected.')                        
          ENDIF                                                         
          RETURN                                                        
      ENDIF                                                             
                                                                        
      IEXIST = 0                                                        
      I  = IABS(MOD( ID, 13 )) + 1                                      
      NS     = DHASH(1, I)                                              
                                                                        
      IF( NS .EQ. 1 ) THEN                                              
            IF( ID .EQ. MAPD( 1, DHASH(2,I))) THEN                      
*               IF( LOCK .NE. 0 ) RETURN                                
                IEXIST = DHASH(2,I)                                     
            ENDIF                                                       
      ELSEIF( NS .GT. 1 ) THEN                                          
          DO 100 K = 2, DHASH(1,I)+1                                    
            IF( ID .EQ. MAPD( 1, DHASH(K,I))) THEN                      
*               IF( LOCK .NE. 0 ) RETURN                                
                IEXIST = DHASH(K,I)                                     
                GO TO 110                                               
            ENDIF                                                       
  100    CONTINUE                                                       
  110    CONTINUE                                                       
      ENDIF                                                             
      XMIN  = DXMIN*1.0                                                 
      XMAX  = DXMAX*1.0                                                 
      YMIN  = DYMIN*1.0                                                 
      YMAX  = DYMAX*1.0                                                 
                                                                        
      IF( IEXIST .GT. 0 ) THEN                                          
          IF( LU .GT. 0 ) THEN                                          
            WRITE(LU,9100) ID                                           
          ENDIF                                                         
 9100     FORMAT(1X,'Scat_Plot ID (',I3,' ) exists already.')           
          IP1    =  MAPD(2,IEXIST)                                      
          IF(( XMIN .EQ. BUFF(IP1))   .AND.                             
     .       ( XMAX .EQ. BUFF(IP1+1)) .AND.                             
     .       ( NXBIN .EQ. IBUF(IP1+2)) )    THEN                        
             IF(( YMIN .EQ. BUFF(IP1+4))   .AND.                        
     .          ( YMAX .EQ. BUFF(IP1+5)) .AND.                          
     .          ( NYBIN .EQ. IBUF(IP1+6)) )    THEN                     
                  IF( LU .GT. 0 ) THEN                                  
                      WRITE(LU,9110)                                    
                  ENDIF                                                 
 9110             FORMAT(1X,' This call is neglected.')                 
                  RETURN                                                
             ENDIF                                                      
          ENDIF                                                         
          IF( LU .GT. 0 ) THEN                                          
              WRITE(LU,9120) ID,XMIN,XMAX,NXBIN,YMIN,YMAX,NYBIN         
          ENDIF                                                         
 9120     FORMAT(1X,'Scat_Plot ( ID =',I3,' ) parameters are replaced', 
     .          /1X,'by the following new parameters :',                
     .          /1X,' XMIN(',E12.5,')  XMAX(',E12.5,' )  XBIN(',I4,' )',
     .          /1X,' YMIN(',E12.5,')  YMAX(',E12.5,' )  YBIN(',I4,' )')
      ENDIF                                                             
      IF(NXBIN .GT. 50 .OR. NYBIN .GT. 50 ) THEN                        
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9300) NXBIN,NYBIN,ID                              
         ENDIF                                                          
 9300    FORMAT(1X,'Bin size (',2I3,' )  exceeds 50 at ID =',I5,        
     .         /1X,' This call is neglected .')                         
         RETURN                                                         
      ELSEIF((XMIN .GE. XMAX) .OR. (YMIN .GE. YMAX)) THEN               
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9400) ID,XMIN,XMAX,YMIN,YMAX                      
         ENDIF                                                          
 9400    FORMAT(1X,'Lower limit is larger than upper at SC_PL ID =',I5, 
     .         /1X,' This call is neglected .',                         
     .         /1X,' XMIN =',G13.4,' XMAX =',G13.4,                     
     .         /1X,' YMIN =',G13.4,' YMAX =',G13.4)                     
         RETURN                                                         
      ENDIF                                                             
      IF(DHASH(1,I) .GE. NSC ) THEN                                     
         IF( LU .GT. 0 ) THEN                                           
             WRITE(LU,9500) I                                           
         ENDIF                                                          
 9500    FORMAT(1X,I5,'-th Hash table overflow',                        
     .         /1X,' This call is neglected.')                          
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IF( IEXIST .GT. 0 ) THEN                                          
          NSCT     = IEXIST                                             
      ELSE                                                              
          NSCAT        = NSCAT + 1                                      
          DHASH(1,I)   = DHASH(1,I) + 1                                 
          K            = DHASH(1,I) + 1                                 
          DHASH(K,I)   = NSCAT                                          
          NSCT         = NSCAT                                          
          IP1    = NW + 1                                               
          NW  = NW + 2527                                               
          MAPD(1,NSCT)  = ID                                            
          MAPD(2,NSCT)  = IP1                                           
      ENDIF                                                             
                                                                        
         BUFF(IP1     ) = XMIN                                          
         BUFF(IP1 +  1) = XMAX                                          
         IBUF(IP1 +  2) = NXBIN                                         
         DEV            = XMAX - XMIN                                   
         BUFF(IP1 +  3) = DEV/NXBIN                                     
         BUFF(IP1 +  4) = YMIN                                          
         BUFF(IP1 +  5) = YMAX                                          
         IBUF(IP1 +  6) = NYBIN                                         
         DEV            = YMAX - YMIN                                   
         BUFF(IP1 +  7) = DEV/NYBIN                                     
      IP2   = IP1 + 8                                                   
         MAPD(3,NSCT)  = IP2                                            
         IBUF(IP2     ) = 0                                             
      IP3   = IP1 + 2509                                                
         MAPD(4,NSCT)  = IP3                                            
         IBUF(IP3     ) =  0                                            
         IBUF(IP3 +  1) =  0                                            
                                                                        
         I1   = IP3 + 2                                                 
         I2   = I1 + 15                                                 
         NAME = TNAME                                                   
         READ(NAME,9800) (BUFF(I),I=I1,I2)                              
 9800    FORMAT(16A4)                                                   
                                                                        
      RETURN                                                            
      END                                                               

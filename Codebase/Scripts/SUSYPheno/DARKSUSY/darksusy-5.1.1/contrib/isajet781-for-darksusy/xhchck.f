************************************************************************
*    =======================                                           *
      SUBROUTINE XHCHCK(LU)                                             
*    =======================                                           *
* ((Purpose))                                                          *
*      To check the contents of the histogram table                    *
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
      CHARACTER*1 CN                                                    
                                                                        
      CN  = CHAR(12)                                                    
                                                                        
      IF( IPNT .EQ. 0 ) THEN                                            
          WRITE(LU,9000)                                                
      ELSE                                                              
          WRITE(LU,9010) CN                                             
      ENDIF                                                             
 9000 FORMAT(/1H1)                                                      
 9010 FORMAT(A1)                                                        
                                                                        
      WRITE(LU,9050) NW                                                 
 9050 FORMAT(                                                           
     . //5X,'*********  Contents of the histogram Header *********',    
     .     //1X,'(1) Actual Buffer size     = ',I6,' Words')            
      WRITE(LU,9100) NHS,NHIST                                          
 9100 FORMAT(1X,'(2) Contents of Histograms ',                          
     .      /1X,'    Max. No. of Histograms = ',I6,                     
     .      /1X,'    Number   of Histograms = ',I6)                     
                                                                        
      IF( NHIST .GT. 0 ) THEN                                           
          WRITE(LU,9200)                                                
 9200     FORMAT(1X,'   ID     X_min        X_max    X_bin',            
     .              ' Hash Hst#')                                       
          DO 200 I = 1, 13                                              
             NT    = XHASH(1,I)                                         
             IF( NT .GT. 0 ) THEN                                       
                 DO 100 J = 2, NT+1                                     
                    K     = XHASH(J,I)                                  
                    IP1   = MAPL(2,K)                                   
                    IP3   = MAPL(4,K)                                   
                    XMIN  = BUFF(IP1)                                   
                    XMAX  = BUFF(IP1+1)                                 
                    NBIN  = IBUF(IP1+2)                                 
                    WRITE(LU,9300) MAPL(1,K),XMIN,XMAX,NBIN,I,NT,K      
 9300               FORMAT(1X,I5,1X,1PE12.4,1X,E12.4,I5,2I3,I5)         
  100            CONTINUE                                               
             ENDIF                                                      
  200     CONTINUE                                                      
      ENDIF                                                             
                                                                        
      WRITE(LU,9400) NSC,NSCAT                                          
 9400 FORMAT(1X,'(3) Contents of Scatter Plots',                        
     .      /1X,'    Max. No. of Scat_Plots = ',I6,                     
     .      /1X,'    Number   of Scat_Plots = ',I6)                     
                                                                        
      IF( NSCAT .GT. 0 ) THEN                                           
          WRITE(LU,9500)                                                
 9500     FORMAT(1X,'   ID      X_min   ',                              
     .              '     X_max   X-Bin    Y_min   ',                   
     .              '     Y_max   Y_Bin Hash Hst#')                     
          DO 400 I = 1, 13                                              
             NT    = DHASH(1,I)                                         
             IF( NT .GT. 0 ) THEN                                       
                 DO 300 J = 2, NT+1                                     
                    K     = DHASH(J,I)                                  
                    IP1   = MAPD(2,K)                                   
                    IP3   = MAPD(4,K)                                   
                    XMIN  = BUFF(IP1)                                   
                    XMAX  = BUFF(IP1+1)                                 
                    NXBN  = IBUF(IP1+2)                                 
                    YMIN  = BUFF(IP1+4)                                 
                    YMAX  = BUFF(IP1+5)                                 
                    NYBN  = IBUF(IP1+6)                                 
                    WRITE(LU,9600) MAPD(1,K),XMIN,XMAX,NXBN,            
     .                            YMIN,YMAX,NYBN,I,NT,K                 
 9600               FORMAT(1X,I5,1X,1PE12.4,1X,E12.4,I5,                
     .                                 E12.4,1X,E12.4,I5,2I3,I5)        
  300            CONTINUE                                               
             ENDIF                                                      
  400     CONTINUE                                                      
      ENDIF                                                             
      RETURN                                                            
      END                                                               

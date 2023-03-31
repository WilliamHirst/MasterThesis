************************************************************************
*    ====================                                              *
      SUBROUTINE SHRSET                                                 
*    ====================                                              *
* ((Function))                                                         *
*     To reset the content of histograms and scatter plots.            *
* ((Author))                                                           *
*     S.Kawabata   June '90 at KEK                                     *
*                                                                      *
* **********************************************************************
                                                                        
      PARAMETER ( NHS = 50, NSC = 50 )                                  
      INTEGER*4 XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD                
      COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),        
     .              NHIST, MAPL(4,NHS),                                 
     .              NSCAT, MAPD(4,NSC),                                 
     .              NW                                                  
      COMMON /PLOTB/ IBUF( 281*NHS + 2527*NSC )                         
      REAL*4         BUFF( 281*NHS + 2527*NSC )                         
      EQUIVALENCE (IBUF(1),BUFF(1))                                     
                                                                        
      IF( NHIST .GT. 0 ) THEN                                           
         DO 100 IHIST = 1, NHIST                                        
            IP2       = MAPL(3,IHIST) + 52                              
            IP3       = MAPL(4,IHIST)                                   
            IBUF(IP3) = -1                                              
            DO 100 I = 0,51                                             
               BUFF(I+IP2) = 0.0                                        
  100      CONTINUE                                                     
      ENDIF                                                             
C                                                                       
      IF( NSCAT .GT. 0 ) THEN                                           
         DO 400   ISCAT = 1, NSCAT                                      
            IP3         = MAPD(4,ISCAT)                                 
            IBUF(IP3)   = 0                                             
            IBUF(IP3+1) = 0                                             
            IP2         = MAPD(3,ISCAT)                                 
            IBUF(IP2)   = 0                                             
            DO 400   I  = IP2+1,IP2+2500                                
               BUFF(I)  = 0.0                                           
  400      CONTINUE                                                     
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               

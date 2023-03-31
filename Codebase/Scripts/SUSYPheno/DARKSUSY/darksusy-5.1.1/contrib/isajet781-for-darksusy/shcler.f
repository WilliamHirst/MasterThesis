************************************************************************
*    ===================                                               *
      SUBROUTINE SHCLER                                                 
*    ===================                                               *
* ((FUNCTION))                                                         *
*     To cancel the update of histograms and scatter plots in case     *
*   of the trial was rejected.                                         *
* ((Author))                                                           *
*     S.Kawabata June '90 at KEK                                       *
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
                                                                        
      IF( NHIST .GT. 0 ) THEN                                           
         DO 200  J   = 1, NHIST                                         
           IP3       = MAPL(3,J)                                        
           IBUF(IP3) = -1                                               
  200    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      IF( NSCAT .GT. 0 ) THEN                                           
         DO 500   K    = 1, NSCAT                                       
           IP3         = MAPD(4,K)                                      
           IBUF(IP3)   =  0                                             
           IBUF(IP3+1) =  0                                             
  500    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               

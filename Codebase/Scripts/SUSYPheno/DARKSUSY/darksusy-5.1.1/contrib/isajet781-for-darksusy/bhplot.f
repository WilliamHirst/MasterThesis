************************************************************************
*     =========================                                        *
       SUBROUTINE BHPLOT( LU )                                          
*     =========================                                        *
* ((Purpose))                                                          *
*     Interface routine to print histograms and scatter plots.         *
*     Routines XHPLOT and DHPLOT are called to print them.             *
* ((Author))                                                           *
*     S.Kawabata  June '90  at KEK                                     *
************************************************************************
                                                                        
      PARAMETER ( NHS = 50, NSC = 50 )                                  
      INTEGER*4 XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD                
      COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),        
     .              NHIST, MAPL(4,NHS),                                 
     .              NSCAT, MAPD(4,NSC),                                 
     .              NW                                                  
      REAL*8 AVGI,SD,CHI2A                                              
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF                        
*                                                                      *
*--------------------------- Entry point ------------------------------*
*                                                                       
      IF( ITF .LE. 0 ) RETURN                                           
*    ===================                                               *
      CALL XHCHCK( LU )                                                 
*    ===================                                                
                                                                        
      IF( NHIST .LE. 0 ) THEN                                           
         WRITE(LU,9000)                                                 
 9000    FORMAT(1X,'No Histogram')                                      
      ELSE                                                              
         DO 500 J = 1, NHIST                                            
            IFBASE(J) = 1                                               
*          =====================                                        
            CALL XHPLOT(LU, 0, J )                                      
*          =====================                                        
  500    CONTINUE                                                       
      ENDIF                                                             
                                                                        
*    ===================                                                
      CALL DHPLOT( LU )                                                 
*    ===================                                                
                                                                        
      RETURN                                                            
      END                                                               

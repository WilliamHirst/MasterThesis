************************************************************************
*    ============================                                      *
      SUBROUTINE SHINIT( MXTRY )                                        
*    ============================                                      *
* ((Function))                                                         *
*     To clear the histogram buffer for generation efficiency          *
* ((Input))                                                            *
*    MXTRY: Maximum number of trials for one event generation          *
* ((Author))                                                           *
*    S.Kawabata    April 1994                                          *
*                                                                      *
************************************************************************
                                                                        
      INTEGER MXTRY                                                     
      PARAMETER ( MXBIN = 51 )                                          
      COMMON/PLOTSP/ NBIN,IBUFSP( MXBIN )                               
                                                                        
      IF( MXTRY .GT. 50 ) THEN                                          
          NBIN  = 50                                                    
      ELSE                                                              
          NBIN  = MXTRY                                                 
      ENDIF                                                             
                                                                        
      DO 100 I = 1,NBIN+1                                               
         IBUFSP(I) = 0                                                  
  100 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               

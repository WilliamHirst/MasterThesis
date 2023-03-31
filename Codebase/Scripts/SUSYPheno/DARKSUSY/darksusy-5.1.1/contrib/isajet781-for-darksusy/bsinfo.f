*********************************************************************** 
*============================                                         * 
      SUBROUTINE BSINFO( LU )                                           
*============================                                         * 
*((Purpose))                                                          * 
*    Print the information for                                        * 
*        (1) BASES parameters                                         * 
*        (2) Computer time information                                * 
*        (3) Convergency behavior of the Grid optimization step       * 
*        (4) Convergency behavior of the integration step             * 
*(( Input ))                                                          * 
*    LU  :  Logical unit number of printer                            * 
*                                                                     * 
*           by S.Kawabata    March 1994 at KEK                          
*                                                                     * 
*********************************************************************** 
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*4 STIME                                                      
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF                        
                                                                        
*  Print Title                                                          
                                                                        
      CALL BSPRNT( LU, 1, IDUM1, IDUM2 )                                
                                                                        
*  Print Bases parameters                                               
                                                                        
      CALL BSPRNT( LU, 4, IDUM1, IDUM2 )                                
                                                                        
*  Print Computing time information                                     
                                                                        
      CALL BSPRNT( LU, 3, IDUM1, IDUM2 )                                
                                                                        
*  Print Convergency Behaviors                                          
                                                                        
      DO 100 ISTEP = 0, 1                                               
         ITX  = ITG                                                     
         IF( ISTEP .EQ. 1 ) ITX = ITF                                   
                                                                        
      IF( ITX .GT. 0 ) THEN                                             
                                                                        
         CALL BSPRNT( LU, 8, ITX, ISTEP )
	                                                                
                                                                        
      ENDIF                                                             
  100 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               

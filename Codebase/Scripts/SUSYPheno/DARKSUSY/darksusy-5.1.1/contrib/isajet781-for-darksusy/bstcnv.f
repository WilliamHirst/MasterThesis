************************************************************************
*=================================================                      
      SUBROUTINE BSTCNV( TIME, IH, MN, IS1, IS2 )                       
*=================================================                      
* (Purpose)                                                             
*    Resolve TIME in second into IH, MN, IS1, IS2                       
* (Input)                                                               
*    TIME : in the unit of second                                       
* (Output)                                                              
*    IH   : Hours                                                       
*    MN   : Minute                                                      
*    IS1  : Second                                                      
*    IS2  : 0.xx Second                                                 
* (Author)                                                              
*    S.Kawabata 1992 June 15                                            
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*4 TIME                                                       
      INTEGER  HOUR                                                     
      DATA HOUR, MINUT, N100/ 360000, 6000, 100 /                       
                                                                        
      ISEC  = TIME*N100                                                 
      IH    = 0                                                         
      MN    = IH                                                        
      IF( ISEC .GE. MINUT ) THEN                                        
          ITIME = ISEC                                                  
          IF( ISEC .GE. HOUR ) THEN                                     
              IH    = ITIME/HOUR                                        
              IHX   = IH*HOUR                                           
              ITIME = ITIME - IHX                                       
              ISEC  = ISEC - IHX                                        
          ENDIF                                                         
          MN    = ITIME/MINUT                                           
          ISEC  = ISEC - MN*MINUT                                       
      ENDIF                                                             
      IS1  = ISEC/N100                                                  
      IS2  = MOD( ISEC, N100)                                           
                                                                        
      RETURN                                                            
      END                                                               

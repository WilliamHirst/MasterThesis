C********************************************************************** 
C*======================                                              * 
C* FUNCTION DRN( ISEED)                                               * 
C*======================                                              * 
C*  Machine-independent Random number generator                       * 
C*     General purpose Version,  OK as long as >= 32 bits             * 
C*((Arguement))                                                       * 
C*  ISEED: Seed                                                       * 
C*                                                                    * 
C********************************************************************** 
                                                                        
      REAL*8 FUNCTION DRN(ISEED)                                        
                                                                        
      COMMON/RANDM/RDM(31),RM1,RM2,IA1,IC1,M1,IX1,                      
     .                             IA2,IC2,M2,IX2,                      
     .                             IA3,IC3,M3,IX3                       
                                                                        
C Generate Next number in sequence                                      
                                                                        
      IX1    = MOD( IA1*IX1+IC1, M1 )                                   
      IX2    = MOD( IA2*IX2+IC2, M2 )                                   
      IX3    = MOD( IA3*IX3+IC3, M3 )                                   
      J      = 1 + (31*IX3)/M3                                          
      DRN    = RDM(J)                                                   
      RDM(J) = ( FLOAT(IX1)+FLOAT(IX2)*RM2 )*RM1                        
                                                                        
C Omit following statement if function arguement passed by value:       
                                                                        
      ISEED = IX1                                                       
      RETURN                                                            
      END                                                               

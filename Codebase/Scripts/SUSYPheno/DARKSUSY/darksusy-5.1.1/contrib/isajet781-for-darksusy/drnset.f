C********************************************************************** 
C*============================                                        * 
C* Subroutine DRNSET( ISEED )                                         * 
C*============================                                        * 
C*((Purpose))                                                         * 
C*  Initialization routine of                                         * 
C*         Machine-independent Random number generator                * 
C*         General purpose Version,  OK as long as >= 32 bits         * 
C*((Arguement))                                                       * 
C*  ISEED: SEED                                                       * 
C*                                                                    * 
C********************************************************************** 
                                                                        
      SUBROUTINE DRNSET( ISEED )                                        
                                                                        
      COMMON/RANDM/RDM(31),RM1,RM2,IA1,IC1,M1,IX1,                      
     .                             IA2,IC2,M2,IX2,                      
     .                             IA3,IC3,M3,IX3                       
                                                                        
      IA1 =    1279                                                     
      IC1 =  351762                                                     
      M1  = 1664557                                                     
      IA2 =    2011                                                     
      IC2 =  221592                                                     
      M2  = 1048583                                                     
      IA3 =   15091                                                     
      IC3 =    6171                                                     
      M3  =   29201                                                     
                                                                        
C Initialization                                                        
                                                                        
      IX1  = MOD( ISEED, M1 )                                           
      IX1  = MOD( IA1*IX1+IC1, M1 )                                     
      IX2  = MOD( IX1, M2 )                                             
      IX1  = MOD( IA1*IX1+IC1, M1 )                                     
      IX3  = MOD( IX1,M3)                                               
      RM1  = 1./FLOAT(M1)                                               
      RM2  = 1./FLOAT(M2)                                               
      DO 100 J = 1,31                                                   
         IX1   = MOD( IA1*IX1+IC1, M1 )                                 
         IX2   = MOD( IA2*IX2+IC2, M2 )                                 
         RDM(J)= ( FLOAT(IX1)+FLOAT(IX2)*RM2 )*RM1                      
  100 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               

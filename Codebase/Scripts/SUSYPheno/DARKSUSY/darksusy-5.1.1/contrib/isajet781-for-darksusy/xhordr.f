C***********************************************************************
C*                                                                     *
C*=============================================                        *
C*    SUBROUTINE XHORDR( VAL, F2, ORDER, IORDR)                        *
C*=============================================                        *
C*((Function))                                                         *
C*    To resolve the real number VAL into mantester and exponent parts.*
C*  When VAL = 1230.0 is given, output are                             *
C*        F2 = 1.2  and ORDER = 4.0.                                   *
C*((Input))                                                            *
C*  VAL  : Real*4 value                                                *
C*((Output))                                                           *
C*  F2   : The upper two digits is given                               *
C*  ORDER: Order is given                                              *
C*  IORDR: Exponent is given                                           *
C*((Author))                                                           *
C*  S.Kawabata                                                         *
C*                                                                     *
C***********************************************************************
                                                                        
      SUBROUTINE XHORDR(VAL, F2, ORDER, IORDR)                          
                                                                        
      IF( VAL .NE. 0.0 ) THEN                                           
          ORDER    =  LOG10( VAL )                                      
          IORDR    =  INT( ORDER )                                      
          IF( ORDER .LT. 0.0 ) IORDR = IORDR - 1                        
          ORDER  = 10.0**IORDR                                          
          F2     = VAL/ORDER                                            
      ELSE                                                              
          IORDR  = 0                                                    
          ORDER  = 1.0                                                  
          F2    = 0.0                                                   
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

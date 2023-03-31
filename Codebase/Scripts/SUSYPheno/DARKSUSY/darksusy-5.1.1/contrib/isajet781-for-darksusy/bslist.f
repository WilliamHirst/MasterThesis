*********************************************************************** 
*    ===================================                              * 
      SUBROUTINE BSLIST( LU, I, ISTEP )                                 
*    ===================================                              * 
* ((purpose))                                                         * 
*     Print out results of each iteration and cumulative result       * 
* ((Argument))                                                        * 
*  (Input)                                                            * 
*     LU      : Logical unit number for the printer                   * 
*     I       : Address in the arrays of common /BASE5/               * 
*     ISTEP   : The Set-Identifier                                    * 
*               ISTEP = ( 0 / 1 ) = ( Grid opt. / Integration step )  * 
*                                                                     * 
*     S. Kawabata   March '94                                         * 
*********************************************************************** 
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (ITM = 50)                                              
      REAL*4 TIME, EFF, WRONG, TRSLT, TSTD, PCNT                        
      COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),         
     .               WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),      
     .               TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)         
                                                                        
      CALL BSTCNV( TIME(I,ISTEP), IH, MN, IS1, IS2 )                    
                                                                        
      RE  = RESLT(I,ISTEP)                                              
      AC  = ABS(ACSTD(I,ISTEP))                                         
      ARE = ABS(RE)                                                     
      IF( ARE .GE. AC) THEN                                             
          CALL BSORDR( ARE, F2, ORDER, IORDR)                           
      ELSE                                                              
          CALL BSORDR(  AC, F2, ORDER, IORDR )                          
      ENDIF                                                             
      RE  = RE/ORDER                                                    
      AC  = AC/ORDER                                                    
      IEFF = EFF(I,ISTEP)                                               
      WRITE(LU,9631) ITRAT(I,ISTEP),IEFF,WRONG(I,ISTEP),                
     .              TRSLT(I,ISTEP),TSTD(I,ISTEP),                       
     .              RE,AC,IORDR,PCNT(I,ISTEP),IH,MN,IS1,IS2             
 9631 FORMAT(I4,I4,F6.2,1P,E11.3, 0P,1X,F6.3,                           
     .              F10.6,'(+-',F8.6,')E',I3.2,1X,F6.3,                 
     .          1X,I3,':',I2,':',I2,'.',I2.2)                           
                                                                        
                                                                        
      RETURN                                                            
      END                                                               

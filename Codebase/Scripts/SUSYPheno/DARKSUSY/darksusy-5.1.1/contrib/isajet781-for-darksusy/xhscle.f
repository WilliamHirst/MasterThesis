************************************************************************
*    =========================================================         *
      SUBROUTINE XHSCLE( IFLG,VMIN,VMAX,VSTP,UNIT,SCALE,CHAR)           
*    =========================================================         *
* ((Function))                                                         *
*     Determine the vertical scale and make it's format                *
* ((Input))                                                            *
*     IFLG   : Flag which indicates whether logarithmic or linear      *
*              scale.  IFLG = ( 1 / any other ) = ( log / linear )     *
*     VMIN,VMAX : Minimum and maximum values of vertical window.       *
*     VSTEP  : Step of unit scale                                      *
*     UNIT   : Unit of one mark *(or o)                                *
* ((Output))                                                           *
*     NSCL   : Number of scale mark                                    *
*     NBLK   : Number of blanks between scale marks                    *
*     CHAR   : Format of scale                                         *
* ((Author))                                                           *
*     S.Kawabata    Oct '85  at KEK                                    *
*                                                                      *
************************************************************************
                                                                        
      CHARACTER*50 CHAR                                                 
      CHARACTER*52 SCALE                                                
      CHARACTER*1 PLUS,MINUS                                            
      DATA PLUS /'+'/, MINUS /'-'/                                      
                                                                        
C     IFLG =    1 : Log scale                                           
C            other: Linear scale                                        
      WRITE(SCALE,9000)                                                 
 9000 FORMAT(5('          '))                                           
      IF( IFLG .EQ. 1 ) THEN                                            
          SC  = 10.**VMIN                                               
      ELSE                                                              
          SC  = VMIN                                                    
      ENDIF                                                             
                                                                        
      WRITE(SCALE(1:8),9100) SC                                         
 9100 FORMAT(1P,E8.1)                                                   
      I2    = 8                                                         
      STV   = VSTP + VMIN                                               
      STV1  = STV                                                       
      VAL1  = VMIN                                                      
      CHAR(50:50) = PLUS                                                
      DO  100   I = 1, 49                                               
          VAL2    = VAL1 + UNIT                                         
          IF( STV .GE. VAL1 .AND. STV .LT. VAL2 ) THEN                  
              CHAR(I:I)  = PLUS                                         
              NSCL       = NSCL + 1                                     
              IF( IFLG .EQ. 1 ) THEN                                    
                 SC          = 10.0**STV                                
              ELSE                                                      
                 IF(     STV1 .EQ. 0.0 ) THEN                           
                         SC           = STV                             
                 ELSEIF( ABS(STV/STV1) .LT. 1.E-2 ) THEN                
                         SC           = 0.0                             
                 ELSE                                                   
                         SC          = STV                              
                 ENDIF                                                  
                 STV1       = STV                                       
              ENDIF                                                     
              STV  = STV + VSTP                                         
              IF( I2 .LT. I-1 ) THEN                                    
                  I2   = I + 8                                          
                  IF( I2 .LE. 52 ) THEN                                 
                      WRITE(SCALE(I+1:I2),9100) SC                      
                  ENDIF                                                 
              ENDIF                                                     
          ELSE                                                          
              CHAR(I:I) = MINUS                                         
          ENDIF                                                         
          VAL1      = VAL2                                              
  100 CONTINUE                                                          
C                                                                       
      IF( NSCL .EQ. 0 ) THEN                                            
          IF( IFLG .EQ. 1 ) THEN                                        
             SC       = 10.0**VMAX                                      
          ELSE                                                          
             SC       = VMAX                                            
          ENDIF                                                         
          WRITE(SCALE(44:52),9100) SC                                   
      ENDIF                                                             
C                                                                       
      RETURN                                                            
      END                                                               

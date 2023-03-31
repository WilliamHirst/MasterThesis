C***********************************************************************
C*                                                                     *
C*============================================================         *
C*  SUBROUTINE XHRNGE( IFLG, VMIN, VMAX, VTMIN, VTMAX, STEP)           *
C*============================================================         *
C*((Function))                                                         *
C*    Determine the vertical range of the histogram.                   *
C*((Input))                                                            *
C*    IFLG   : Flag which indicates whether logarithmic or linear      *
C*             scale.  IFLG = ( 1 / any other ) = ( log / linear )     *
C*    VMIN,VMAX : Minimum and maximum values of vertical window.       *
C*((Output))                                                           *
C*    VTMIN,VTMAX : Minimum and maxmum values of optimized vertical    *
C*                  window.                                            *
C*    STEP   : step of scale for the optimized vertical window         *
C*((Author))                                                           *
C*    S.Kawabata    Oct '85  at KEK                                    *
C*                                                                     *
C***********************************************************************
C                                                                       
      SUBROUTINE XHRNGE( IFLG, VMIN, VMAX, VTMIN, VTMAX, STEP)          
C                                                                       
C     IFLG =    1 : Log scale                                           
C            other: Linear scale                                        
C                                                                       
      PARAMETER ( NBIN  = 25 )                                          
      REAL    WIND(NBIN),STP1(NBIN),STP2(NBIN)                          
C                                                                       
      DATA WIND/                                                        
     .   1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.80, 2.00,  2.20,   
     .   2.50, 2.70, 3.00, 3.30, 3.60, 4.00, 4.50, 5.00, 5.50,  6.00,   
     .   6.50, 7.00, 8.00, 9.00, 10.0/                                  
*     DATA STP1/                                                        
*    .   0.20, 0.22, 0.30, 0.26, 0.28, 0.30, 0.32, 0.36, 0.40,  0.44,   
*    .   0.50, 0.54, 0.60, 0.66, 0.60, 0.80, 0.90, 1.00, 1.10,  1.00,   
*    .   1.30, 1.00, 1.60, 1.80, 2.00/                                  
      DATA STP1/                                                        
     .   0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.450,0.500,0.550,   
     .   0.625,0.675,0.750,0.825,0.900,1.000,1.125,1.250,1.375,1.500,   
     .   1.625,1.750,2.000,2.250,2.500/                                 
      DATA STP2/                                                        
     .   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  1.00,   
     .   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  2.00,   
     .   2.00, 2.00, 2.00, 2.00, 2.00/                                  
C                                                                       
          XMAX   = VMAX                                                 
          XMIN   = VMIN                                                 
          IFLAG  = IFLG                                                 
          IF( IFLG .NE. 1 .AND. VMIN .LT. 0.0 ) THEN                    
              IF( VMAX .LE. 0.0 )THEN                                   
                  IFLAG  = 2                                            
                  XMAX  = - VMIN                                        
                  XMIN  =  0.0                                          
              ELSE                                                      
                  AVMIN  = - VMIN                                       
                  XMIN  =0.0                                            
                  IF( VMAX .GE. AVMIN ) THEN                            
                      IFLAG  = 3                                        
                      XMAX  = VMAX                                      
                      XMIN1 = AVMIN                                     
                  ELSE                                                  
                      IFLAG  = 4                                        
                      XMAX  = AVMIN                                     
                      XMIN1 = VMAX                                      
                  ENDIF                                                 
              ENDIF                                                     
          ENDIF                                                         
          DSCALE = XMAX - XMIN                                          
          CALL XHORDR( DSCALE, DSF2, DSORDR, IORD)                      
                                                                        
          DO 100 I = 2, 25                                              
             IF( DSF2 .GE. WIND(I-1) .AND.                              
     .           DSF2 .LE. WIND( I )       ) GO TO 200                  
 100      CONTINUE                                                      
          I = 25                                                        
C                                                                       
 200      CONTINUE                                                      
                                                                        
          XMAX = WIND(I)*DSORDR + XMIN                                  
          IF(     DSORDR .GE. 10.0 .OR. IFLG .NE. 1 ) THEN              
                  STEP1  = STP1(I)                                      
                  STEP   = STEP1*DSORDR                                 
          ELSE                                                          
                  STEP1  = STP2(I)                                      
                  STEP   = STEP1                                        
          ENDIF                                                         
                                                                        
          IF(     IFLAG .LE. 1 ) THEN                                   
                  VTMAX  = XMAX                                         
                  VTMIN  = XMIN                                         
          ELSEIF( IFLAG .EQ. 2 ) THEN                                   
                  VTMAX  = XMIN                                         
                  VTMIN  = -XMAX                                        
          ELSE                                                          
                                                                        
                  XPLUS   = 0.0                                         
                  DO 300 J = 1, 10                                      
                     XPLUS = XPLUS + STEP                               
                     IF( XPLUS .GT. XMIN1 ) GO TO 400                   
 300              CONTINUE                                              
 400              XMIN = XPLUS                                          
                  XMAX = XMAX                                           
                  IF( IFIX((WIND(I)+0.1)/STEP1)+J .GT. 7 ) THEN         
                      STEP = 2.0*STEP                                   
                  ENDIF                                                 
                  IF( IFLAG .EQ. 3 ) THEN                               
                      VTMAX  = XMAX                                     
                      VTMIN  = -XMIN                                    
                  ELSE                                                  
                      VTMAX  = XMIN                                     
                      VTMIN  = -XMAX                                    
                  ENDIF                                                 
          ENDIF                                                         
C                                                                       
      RETURN                                                            
      END                                                               

************************************************************************
*    ===================                                               *
      SUBROUTINE BSINIT                                                 
*    ===================                                               *
* ((Purpose))                                                          *
*     Initialization of BASE50/SPRING50.                               *
*     Function of this routine is                                      *
*       (0) Set the size of histogram and scatter plot buffers         *
*       (1) Set the parameters INTV and IPNT                           *
*             INTV = ( 0 / 1 / any )                                   *
*                  = ( Batch / Batch(Unix) / Interactive )             *
*             IPNT = ( 0 / any )                                       *
*                  = ( IBM Type / Ascii printer )                      *
*       (2) Set the acceleration factor ALPHA by 1.5                   *
*            The range of this value is from 0.0 to 2.0.               *
*            ALPHA = 0.0 results in no grid-optimization.              *
*       (3) Set the grid-optimization flag IGOPT ( Default value 0 )   *
*             IGOPT = 0  :  The grid is optimized by VEGAS algorithm   *
*             IGOPT = 1  :  The grid is optimized so that the accuracy *
*                           of each iteration be minimized.            *
*       (4) Set Node-ID number NODEID and the number of nodes NUMNOD   *
*       (5) Set seed of radom number                                   *
*       (6) Set the values of BASES paremeters with default ones.      *
*       (7) Set the values of parameters with non-sense values,        *
*            which should be set again with the true values by User    *
*            before running BASES.                                     *
*                                                                      *
*        Coded by S.Kawabata         March '94                         *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (MXDIM = 50, NDMX = 50 )                                
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                   
     .               IG(MXDIM),NCALL                                    
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2                             
                                                                        
      COMMON /BASE0/ JFLAG,IBASES                                       
      COMMON /BASE6/ D(NDMX,MXDIM),                                     
     .               ALPH,XSAVE(NDMX,MXDIM),XTI,XTSI,XACC,IGOPT         
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
       COMMON/NINFO/ NODEID, NUMNOD                                     
       COMMON /BDATE/ IDATE(3),ITIME(2)                                 
*            IDATE(1) : year        ITIME(1) : hour                     
*            IDATE(2) : month       ITIME(2) : minute                   
*            IDATE(3) : day                                             
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1      
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1                  
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)                    
*=========================================================              
* (0) Initialization of timer and Histogram buffer                      
*     Timer initialization                                              
       CALL BSTIME( TIME0, 0 )                                          
       TIMEB1 = TIME0                                                   
       TIMINT = 0                                                       
                                                                        
*     Histogram buffer initialization                                   
       LU  = 6                                                          
       CALL BHINIT( LU )                                                
                                                                        
*=========================================================              
                                                                        
* (1) Set the parameters INTV and IPNT                                  
       INTV  = 0                                                        
       IPNT  = 1                                                        
* (2) Set the acceleration factor ALPHA by 1.5                          
       ALPH  = 1.5D0                                                    
* (3) Set the grid-optimization flag IGOPT                              
       IGOPT = 0                                                        
* (4) Set Node-ID number NODEID and the number of nodes NUMNOD          
*      IF( INTV .EQ. 0 ) THEN                                           
           NODEID = 0                                                   
           NUMNOD = 1                                                   
*      ELSE                                                             
*          NODEID = 0                                                   
*          NUMNOD = 1                                                   
*      ENDIF                                                            
                                                                        
C---------------------------------------------------------------        
C (5)  Set initial seeds of random number generator                     
C---------------------------------------------------------------        
       ISEED = 12345                                                    
C                                                                       
       CALL DRNSET( ISEED )                                             
C ---------------------------------------------------------------       
C (6),(7)  Set BASES parameters equal to default values                 
C ---------------------------------------------------------------       
C                                                                       
cccAB       NDIM   = -1                                                      
cccAB        NWILD  =  1                                                      
cccAB        ITMX1  = 15                                                      
cccAB        ITMX2  = 100                                                     
cccAB        NCALL  = 1000                                                    
cccAB        ACC1   = 0.2D0                                                   
cccAB        ACC2   = 0.01D0                                                  
       DO 100 I = 1,MXDIM                                               
          IG(I) = 1                                                     
          XU(I)  = -1.0D37                                              
  100  CONTINUE                                                         
                                                                        
*    Initialization of computing time table of BASES                    
       DO 200 I = 0, 2                                                  
          TIMEBS(I) = 0.0                                               
  200  CONTINUE                                                         
                                                                        
*-------------------------------------------                            
*      Don't change IBASES from this value                              
*-------------------------------------------                            
       IBASES =  1                                                      
                                                                        
       RETURN                                                           
       END                                                              

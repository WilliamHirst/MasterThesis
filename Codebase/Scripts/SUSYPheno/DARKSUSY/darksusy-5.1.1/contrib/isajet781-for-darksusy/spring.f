************************************************************************
*    ==================================                                *
      SUBROUTINE SPRING(FUNC, MXTRY )                                   
*    ==================================                                *
*         Main Program for the Event generation program SPRING.        *
*                                                                      *
*        Coded by S.Kawabata        September '84                      *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      EXTERNAL FUNC                                                     
      COMMON /BASE0/ NDUM,IBASES                                        
      PARAMETER (MXDIM = 50, NDMX = 50, LENG = 32768)                   
      COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                    
     .               IG(MXDIM),NCALL                                    
      COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),      
     .               ND,NG,NPG,MA(MXDIM)                                
      COMMON /BDATE/ IDATE(3),ITIME(2)                                  
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
                                                                        
      COMMON /SPRNG1/ XND, DXG, XJAC, DXMAX, NSP                        
      COMMON /SPRNG2/ MXTRYP,NEVENT, NTRIAL,MISS                        
                                                                        
      REAL*4 TIMEBS,TIMINT,TIMESP,TIME0,RTIME,TIMEB1,TIMEB2,TIMES1      
      COMMON /BTIME1/ TIME0,RTIME,TIMEB1,TIMEB2,TIMES1                  
      COMMON /BTIME2/ TIMEBS(0:2),TIMINT,TIMESP(0:2)                    
*                                                                      *
*----------------------------- Entry point ----------------------------*
*                                                                      *
*======================================================================*
*                  Initialization of the program                       *
*======================================================================*
*----------------------------------------------------------------------*
*                     initialize timer etc.                            *
*----------------------------------------------------------------------*
*                                                                      *
       IF( IBASES .GT. 0 ) THEN                                         
                                                                        
           CALL SPCHCK                                                  
                                                                        
           CALL BSTIME( TIME0, 0 )                                      
           TIMES1 = TIME0                                               
                                                                        
           MXTRYP = MXTRY                                               
           INTV   = 0                                                   
           IBASES = 0                                                   
           MISFLG = 0                                                   
                                                                        
           CALL BSDATE                                                  
                                                                        
           DO 10 I = 0,2                                                
              TIMESP(I) = 0.0                                           
   10      CONTINUE                                                     
*                                                                      *
            IF( MXTRY .LT. 10 ) MXTRY = 50                              
            NBIN    = MXTRY                                             
            IF( MXTRY .GT. 50) NBIN = 50                                
            MXTRY1  = MXTRY + 1                                         
            MISS    = 0                                                 
            NEVENT  = 0                                                 
            NTRIAL  = 0                                                 
                                                                        
            CALL SHINIT( MXTRY1 )                                       
                                                                        
*           -------------                                               
             CALL SHRSET                                                
*            -------------                                              
*----------------------------------------------------------------------*
*             Make the cumulative probability distribution             *
*----------------------------------------------------------------------*
*                                                                      *
            XND     = ND                                                
            DXG     = XND/NG                                            
            NSP     = NG**NWILD                                         
                                                                        
*///// DEBUG                                                            
*       MCALL   = NSP*NPG                                               
*       CALL BSPRNT( 4, MCALL, IDUM2, IDUM3, IDUM4 )                    
*                                                                       
            XJAC    = 1.0                                               
            DO 50 I = 1, NDIM                                           
               XJAC = XJAC*DX(I)                                        
   50       CONTINUE                                                    
            DXMAX   = 0.0D0                                             
            DO 100  I = 1,NSP                                           
               IF( DXD( I ) .LT. 0.0D0 ) THEN                           
                   WRITE(6,9100) I                                      
 9100              FORMAT(                                              
     .             /5X,'********** FATAL ERROR IN SPRING **********',   
     .             /5X,'*     Negative probability was found      *',   
     .             /5X,'*        in the ',I6,'-th Hypercube.      *',   
     .             /5X,'*******************************************')   
                   STOP                                                 
               ENDIF                                                    
                                                                        
               DXMAX    = DXMAX + DXD( I )                              
               DXD(I)   = DXMAX                                         
  100       CONTINUE                                                    
*        =====================                                          
          CALL BSUTIM( 1, 1 )                                           
*        =====================                                          
      ENDIF                                                             
*     =====================                                             
       CALL BSUTIM( 1, 2 )                                              
*     =====================                                             
      IF( IBASES .EQ. 1 ) THEN                                          
          WRITE(6,9000)                                                 
 9000     FORMAT(                                                       
     .      1X,'**************************************************',    
     .     /1X,'*    Flag IBASES was not equal to "0".           *',    
     .     /1X,'*                                                *',    
     .     /1X,'*   Process was terminated by this error.        *',    
     .     /1X,'*   Call S.Kawabata.                             *',    
     .     /1X,'**************************************************')    
           STOP                                                         
       ENDIF                                                            
*                                                                      *
*======================================================================*
*                       Event generation                               *
*======================================================================*
*     =====================                                             
  500  CALL BSUTIM( 1, 1 )                                              
*     =====================                                             
                                                                        
*     ==================================                                
        CALL SPRGEN( FUNC, MXTRY, IRET)                                 
*     ==================================                                
                                                                        
*     =====================                                             
       CALL BSUTIM( 1, 0 )                                              
*     =====================                                             
                                                                        
      CALL SHFILL( IRET )                                               
                                                                        
      IF( IRET .LE. MXTRY ) THEN                                        
          NTRIAL =NTRIAL + IRET                                         
          NEVENT = NEVENT + 1                                           
          CALL SHUPDT                                                   
      ELSE                                                              
          NTRIAL =NTRIAL + IRET - 1                                     
          MISS = MISS + 1                                               
          IF( MISFLG .EQ. 0 .AND. MISS .GT. MXTRY ) THEN                
              WRITE(6,9600) MXTRY                                       
 9600         FORMAT(1X,'****************************************',     
     .                  '****************************************',     
     .              /1X,'* (((( Warning ))))                     ',     
     .                  '                                       *',     
     .              /1X,'*                                       ',     
     .                  '                                       *',     
     .              /1X,'*  The number of mis-generations is foun',     
     .                  'd more than',I3,' times.                  *')  
              WRITE(6,9610)                                             
 9610         FORMAT(1X,'*                                       ',     
     .                  '                                       *',     
     .              /1X,'*(( Suggestion ))                       ',     
     .                  '                                       *',     
     .              /1X,'* (1) Try integration again with larger ',     
     .                  'number of sample points than this job. *',     
     .              /1X,'* or                                    ',     
     .                  '                                       *',     
     .              /1X,'* (2) The integral variables are not sui',     
     .                  'ted for the function.                  *',     
     .              /1X,'*     Take another integral variables !!',     
     .                  '                                       *',     
     .              /1X,'*                                       ',     
     .                  '                                       *',     
     .              /1X,'****************************************',     
     .                  '****************************************')     
            MISFLG = 1                                                  
          ENDIF                                                         
          GO TO 500                                                     
      ENDIF                                                             
*     =====================                                             
  600  CALL BSUTIM( 1, 1 )                                              
*     =====================                                             
                                                                        
      RETURN                                                            
      END                                                               

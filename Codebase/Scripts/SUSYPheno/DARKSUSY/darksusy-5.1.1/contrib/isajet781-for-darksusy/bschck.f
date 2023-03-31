************************************************************************
*    ===================                                               *
      SUBROUTINE BSCHCK                                                 
*    ===================                                               *
* ((Purpose))                                                          *
*     To check user's initialization parameters.                       *
*                                                                      *
*        Coded by S.Kawabata        Oct. '85                           *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER ( MXDIM = 50)                                           
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                   
     .               IG(MXDIM),NCALL                                    
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2                             
                                                                        
      COMMON /BASE0/ JFLAG,IBASES                                       
      COMMON /BASE1/ XLT(MXDIM),XUT(MXDIM),NDIMT,NWILDT,                
     .               IGT(MXDIM),NCALLT                                  
      COMMON /BASE2/ ACC1T,ACC2T,ITMX1T,ITMX2T                          
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
      COMMON /XHCNTL/ LOCK                                              
                                                                        
      LOCK  = 1                                                         
                                                                        
      IF( IBASES .NE.  1 ) THEN                                         
          WRITE(6,9000)                                                 
 9000     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   BSINIT was not called before calling BASES  *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
          STOP                                                          
      ENDIF                                                             
                                                                        
                                                                        
      IF( NDIM .LT. 1) THEN                                             
          WRITE(6,9100)                                                 
 9100     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   NDIM was not set before calling BASES.      *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
          STOP                                                          
      ENDIF                                                             
                                                                        
      NDIMT = NDIM                                                      
                                                                        
      DO 200 I = 1,NDIM     
         IF( XU(I) .LE. -1.0D37) THEN                                   
             WRITE(6,9200) I,I                                          
 9200        FORMAT(                                                    
     .        5X,'*************************************************',   
     .       /5X,'*                                               *',   
     .       /5X,'*   XL(',I6,' ).  XU(',I6,' ) were not set      *',   
     .       /5X,'*    before calling BASES.                      *',   
     .       /5X,'*   Process was terminated due to this error.   *',   
     .       /5X,'*                                               *',   
     .       /5X,'*************************************************')   
             STOP                                                       
         ENDIF                                                          
                                                                        
         IGT(I)  = IG(I)                                                
         XLT(I)  = XL(I)                                                
         XUT(I)  = XU(I)                                                
                                                                        
  200 CONTINUE                                                          
C                                                                       
C  Change the maximum number of the wild variables                      
C 10 ===> 15                                                            
      IF( NWILD .LT.  0) THEN                                           
          NWILD = MIN( NDIM, 15)                                        
          WRITE(6,9300) NWILD                                           
 9300     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   NWILD was not set before calling BASES.     *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   NWILD is set equal to the value(',I6,' ).   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
      ELSE                                                              
     .IF( NWILD .GT. 15) THEN                                           
          NWILDO = NWILD                                                
          NWILD  = MIN( NDIM, 15)                                       
          WRITE(6,9400) NWILDO, NWILD                                   
 9400     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   NWILD(',I6,' ) was too large number.        *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   NWILD is set equal to the value(',I6,' ).   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
      ENDIF                                                             
                                                                        
      NWILDT = NWILD                                                    
      NCALLT = NCALL                                                    
                                                                        
      ITMX1T = ITMX1                                                    
      ITMX2T = ITMX2                                                    
      ACC1T  = ACC1                                                     
      ACC2T  = ACC2                                                     
C                                                                       
      RETURN                                                            
      END                                                               

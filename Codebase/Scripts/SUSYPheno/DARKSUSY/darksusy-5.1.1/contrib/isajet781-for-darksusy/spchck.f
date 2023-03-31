************************************************************************
*    ===================                                               *
      SUBROUTINE SPCHCK                                                 
*    ===================                                               *
* ((Purpose))                                                          *
*     To check user's initialization parameters.                       *
*                                                                      *
*        Coded by S.Kawabata      April  '94                           *
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
                                                                        
      IF( NDIM .NE. NDIMT ) THEN                                        
          WRITE(6,9100) NDIM,NDIMT                                      
 9100     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Given NDIM(',I6,' ) does not match          *',      
     .    /5X,'*      to NDIM(',I6,' ) in BASES.               *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
          STOP                                                          
      ENDIF                                                             
                                                                        
      IF( NWILD .NE. NWILDT ) THEN                                      
          WRITE(6,9110) NWILD,NWILDT                                    
 9110     FORMAT(                                                       
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Given NWILD(',I6,' ) does not match         *',      
     .    /5X,'*      to NWILD(',I6,' ) in BASES.              *',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
          STOP                                                          
      ENDIF                                                             
                                                                        
      DO 200 I = 1,NDIM                                                 
         IF( XL(I) .NE. XLT(I) ) THEN                                   
             WRITE(6,9200) I,XL(I),I,XLT(I)                             
 9200        FORMAT(                                                    
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Given XL(',I3,' ) = ',D15.8,'            *',         
     .    /5X,'*      does not match to                        *',      
     .    /5X,'*      to XL(',I3,' ) = ',D15.8,' in BASES   *',         
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
             STOP                                                       
         ENDIF                                                          
         IF( XU(I) .NE. XUT(I) ) THEN                                   
             WRITE(6,9210) I,XU(I),I,XUT(I)                             
 9210        FORMAT(                                                    
     .     5X,'*************************************************',      
     .    /5X,'*                                               *',      
     .    /5X,'*   Given XU(',I3,' ) = ',D15.8,'            *',         
     .    /5X,'*      does not match to                        *',      
     .    /5X,'*      to XU(',I3,' ) = ',D15.8,' in BASES   *',         
     .    /5X,'*                                               *',      
     .    /5X,'*   Process was terminated due to this error.   *',      
     .    /5X,'*                                               *',      
     .    /5X,'*************************************************')      
             STOP                                                       
         ENDIF                                                          
  200 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               

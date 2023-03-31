************************************************************************
*    ================================================                  *
      SUBROUTINE BSPARM( MCALL, AC1, AC2, IT1, IT2 )                    
*    ================================================                  *
* ((Purpose))                                                          *
*     To set the BASES parameters.                                     *
* ((Input))                                                            *
*     MCALL  : The number of sample points per iteration.              *
*              This actual number is to be determined by taking the    *
*              number of dimensions into account.                      *
*     AC1 %  : The required accuracy at the grid optimization step     *
*     AC2 %  : The required accuracy at the integration step.          *
*     IT1    : The max. number of iteration at the grid opt. step.     *
*     IT2    : Thr max. number of iteration at the integration step.   *
* ((Output))                                                           *
*     These parameters are set in the labeled common /BPARM1/ and      *
*     /BPARM2/.                                                         
*                                                                      *
*        Coded by S.Kawabata         August '94                        *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (MXDIM = 50, NDMX = 50 )                                
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                   
     .               IG(MXDIM),NCALL                                    
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2                             
                                                                        
      INTEGER MCALL, IT1, IT2                                           
      REAL*8 AC1, AC2                                                   
                                                                        
      NCALL = MCALL                                                     
      ACC1  = AC1                                                       
      ACC2  = AC2                                                       
      ITMX1 = IT1                                                       
      ITMX2 = IT2                                                       
                                                                        
      RETURN                                                            
      END                                                               

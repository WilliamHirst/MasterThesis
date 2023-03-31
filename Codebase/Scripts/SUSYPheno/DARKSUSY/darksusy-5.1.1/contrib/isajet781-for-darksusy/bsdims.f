************************************************************************
*    ===========================================                       *
      SUBROUTINE BSDIMS( MDIM, MWILD, XLL, XUU )                        
*    ===========================================                       *
* ((Purpose))                                                          *
*     To set the BASES parameters.                                     *
* ((Input))                                                            *
*     MDIM   : The number of dimension of integral                     *
*     MWILD  : The number of wild variables                            *
*     XLL(i) : The lower value of the i-th integral variable           *
*     XUU(i) : The upper value of the i-th integral variable           *
* ((Output))                                                           *
*     These parameters are to be set in the labeled common /BPARM1/    *
* (Caution)                                                            *
*     The parameter IG(i) is not able to set by this routine.          *
*     If some of parameters IG(i) are required to be changed,          *
*     it is done by calling the subroutine BSGRID.                     *
*                                                                      *
*        Coded by S.Kawabata         August '94                        *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (MXDIM = 50, NDMX = 50 )                                
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                   
     .               IG(MXDIM),NCALL                                    
                                                                        
      DIMENSION XLL(MDIM), XUU(MDIM)                                    
                                                                        
*=========================================================              
                                                                        
      NDIM   = MDIM                                                     
      NWILD  = MWILD                                                    
      DO 100 I= 1, NDIM                                                 
         XL(I) = XLL(I)                                                 
         XU(I) = XUU(I)                                                 
  100 CONTINUE                                                          
                                                                        
       RETURN                                                           
       END                                                              

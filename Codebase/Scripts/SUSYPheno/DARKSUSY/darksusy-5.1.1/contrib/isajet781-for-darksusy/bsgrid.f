************************************************************************
*    ==================================                                *
      SUBROUTINE BSGRID( MDIM, IGG )                                    
*    ==================================                                *
* ((Purpose))                                                          *
*     To change the grid optimizing flag.                              *
* ((Input))                                                            *
*     MDIM   : The number of dimension of integral                     *
*     IGG(i) : The flag switches whether the grid of i-th variable     *
*              is to be optimized ( 1 ) or kept uniform ( 0 ).         *
* ((Output))                                                           *
*     These parameters are to be set in the labeled common /BPARM1/    *
*                                                                      *
*        Coded by S.Kawabata         August '94                        *
*                                                                      *
************************************************************************
                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (MXDIM = 50, NDMX = 50 )                                
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,                   
     .               IG(MXDIM),NCALL                                    
                                                                        
      DIMENSION IGG(MDIM)                                               
                                                                        
*=========================================================              
                                                                        
      NDIM   = MDIM                                                     
      DO 100 I= 1, NDIM                                                 
         IG(I) = IGG(I)                                                 
  100 CONTINUE                                                          
                                                                        
       RETURN                                                           
       END                                                              

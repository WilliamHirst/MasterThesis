

C-------------------------------------------------------
C#######################################################
C-------------------------------------------------------

      FUNCTION FUNC_INT(IPRINT)
      IMPLICIT NONE
C------- BASES COMMON BLOCKS ---------------------
      
      EXTERNAL FUNC
      REAL*8  FUNC_INT,FUNC
      
      INTEGER NDIM,IG,MXDIM,NWILD,NCALL
      REAL*8 XL,XU
      PARAMETER (MXDIM = 50 )                                           
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL    
c      DATA NWILD/2/,NCALL/1000/      
      INTEGER INTV, IPNT, NLOOP, MLOOP,ITMX1,ITMX2
      REAL*8 ACC1,ACC2
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2                             
c      DATA  ITMX1/5/,ITMX2/5/,ACC1/0.1d0/,ACC2/0.1d0/
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP      
c---------------------------------------------------
      REAL*8 ESTIM, SIGMA, CTIME

     
      COMMON /CTRL/
     &AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX,
     &NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP

      REAL*8 AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX
      INTEGER NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP
     
      INTEGER I,IPRINT,IT1,IT2
      
      CALL BSINIT 
      
      ITMX1=5
      ITMX2=5
      NDIM =NDIMUSER_EFF
      NWILD=NDIM
      ACC1=1.
      ACC2=1.
      
      DO 1 I=1,NDIMUSER_EFF
        XL(I)=0d0
        XU(I)=1d0
    1 CONTINUE
      XL(3)=0.000001
      XU(3)=XF
   
c-------------------------------------------------

      INTV=IPRINT
      
c       print *,'NCALLS=', NCALL,IPRINT,INTV,NPROC_MIN,NPROC_MAX
      
c      CALL BSINIT                                                  
      CALL BASES( FUNC, ESTIM, SIGMA, CTIME, IT1, IT2 )             
      
      FUNC_INT=ESTIM
      
     
      
      RETURN
      END

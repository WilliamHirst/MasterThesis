cccccccccccccccccccccccc===================================cccccccccccc

      Function FA12(as,I,x)
      IMPLICIT NONE
      REAL*8 FA12
      
      REAL*8 A
      COMMON/VARS/A(1800)

      COMMON /CTRL/
     &AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX,
     &NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP

      REAL*8 AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX
      INTEGER NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP
 
      REAL*8 PI,BesK1,BesK2,bj,bi,sqrlam,as,x,gj,gi,BF
      INTEGER ip1,ip2,JSUM,I

      REAL*8 b,g,pm
      CHARACTER*6 names
      INTEGER J,idof,IPTOT
      COMMON /INPART/ b(29),g(29),pm(29),idof(29),IPTOT
      COMMON /NAMES/ names(29)
      
      CHARACTER*6 PINF1,PINF2,PINF
     
     
      Data Pi / 3.14159265358979323846D0 /

      DO J=1,IPTOT
        IF(PINF(I,1).eq.NAMES(J)) IP1=J
        IF(PINF(I,2).eq.NAMES(J)) IP2=J
        b(J)=PM(J)/PM(1)
 	g(J)=IDOF(J)
      ENDDO
     
  
        BF=0d0
      
       JSUM=IPTOT
       if(NPROC_MIN.ge.1.and.NPROC_MAX.le.26.and.ISUM.eq.0) JSUM=1   
       if(NDIMUSER_EFF.eq.2) JSUM=1
       DO J=1,JSUM
 	BF=BF+b(J)**2*g(J)*BesK2(b(J)/x)*Exp((-b(J)+b(1))/x)
      ENDDO
      
C Changed from /Exp((-2d0*b(1)+as)/x) to avoid floating overflow. FEP

        FA12= 
     _ sqrlam(as,b(ip1),b(ip2))**2*
     _ g(ip1)*g(ip2)*BesK1(as/x) /4d0*Exp(-(-2d0*b(1)+as)/x)/(x*(BF)**2)

       End   

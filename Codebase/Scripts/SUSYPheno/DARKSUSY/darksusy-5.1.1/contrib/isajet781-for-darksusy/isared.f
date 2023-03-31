CDECK  ID>, ISARED.
cccccccccccccccccccccccccccccccccccccccccccccc
c      AUTHORS: Baer,Balazs,Belyaev
c
c      Last modification -> 10/27/2005 A.Belyaev
c      Last modification -> 10/04/2007 A.Belyaev
c      Last modification -> 12/06/2007 A.Belyaev -> sigma*v (v->0) added
c
cccccccccccccccccccccccccccccccccccccccccccccc
     
      SUBROUTINE ISARED(IPRT)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      
	
c----------USER--------------------------------      
      COMMON /CTRL/
     &AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX,
     &NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP

      REAL*8 AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX
      INTEGER NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP
     

      INTEGER IPRINT
      
      DATA   AS_MAX / 3D0 /, XFINI/0.05D0/, NDIMUSER/3/, NST_MAX/10/,
     &      NPROC_MIN/1/,   NPROC_MAX/1820/, NPROC_STEP/1/, IPRINT/1/,
     &      COS_MIN / -0.999999d0 / , COS_MAX / 0.999999d0 /,
     &      CS_MIN  /0d0/, ISUM/1/, SUPPEXPMAX/2d0/
     

      common/printlevel/NPRINT
      INTEGER NPRINT
c-----------------------------------------------      
      COMMON   /GOOD/    NNOGOOD,IALLOW 
      INTEGER  NNOGOOD,IALLOW                                                                                                       
c------BASES------------------------------------\
     
      INTEGER NDIM,IG,MXDIM,NWILD,NCALL
      REAL*8 XL,XU
   
      PARAMETER (MXDIM = 50 )                                           
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,IG(MXDIM),NCALL    
      DATA NWILD/2/,NCALL/1000/      
      
      INTEGER INTV, IPNT, NLOOP, MLOOP,ITMX1,ITMX2
      REAL*8 ACC1,ACC2
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2                             
      DATA  ITMX1/5/,ITMX2/5/,ACC1/1d0/,ACC2/1d0/
      COMMON /BSCNTL/ INTV, IPNT, NLOOP, MLOOP                          
c-----------------------------------------------
      REAL*8 xb,g,pm
      CHARACTER*6 names
      INTEGER J,idof,IPTOT
      COMMON /INPART/ xb(29),g(29),pm(29),idof(29),IPTOT
      COMMON /NAMES/ names(29)

      DATA IPTOT/29/
      DATA  (IDOF(J),J=1,29)
     _/  2,     2 ,     16,
     _   2 ,	1 ,     1,      1,      3,     3, 
     _  	1 ,     1 ,     1,      3,     3,    3,    3,
     _   2 ,    1 ,     1,      1,      3,     3,
     _  	1 ,     1 ,     1,      3,     3,    3,    3/

     
      DATA  (NAMES(J),J=1,29)
     _/'~o1',  '~o2' ,'~g',
     _  '~1-' ,'~e1' , '~e2' , '~e3' ,'~t1','~b1', 
     _         '~n1' , '~n2' , '~n3' ,'~u1','~d1','~c1','~s1',
     _  '~1+' ,'~E1' , '~E2' , '~E3' ,'~T1','~B1', 
     _         '~N1' , '~N2' , '~N3' ,'~U1','~D1','~C1','~S1'/
c-----------------------------------------------
      REAL*4 OMGH2,SIGMA,XFREEZ,FFF_V
      INTEGER NSTEPS
      COMMON/SUGRED/ OMGH2,SIGMA,XFREEZ,NSTEPS,FFF_V
c-----ISAJET-----------------------------------   
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS 
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_hd^2     GSS(14) = M_hu^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = vdq
C     GSS(31) = vuq
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,IGUTST
      SAVE /SUGPAS/
c----------COMPHEP------------------------------      
      COMMON/VARS/A(1800)
      REAL*8 A
c-----------------------------------------------      
c-----------------------------------------------      
      CHARACTER*6 PINF
      REAL*8 PBNTOGEV,PI,SQRT_GN,FFF,XFI,XF_OLD,SUMM,GEFF,
     &FUNC_INT,FFF_TMP,OMEGA,FFF_OLD
      
      INTEGER IPRT,NPROC_MAX_S,NPROC_MIN_S,IPM,IDEL,
     &NDIMUSER_SAVE,NST,NDIMUSER_EFF_SAVE, NCALL_SAV,I
    
      
c-----------------------------------------------      
      REAL*8 del(29)
    
      REAL*8 vvv,AS_MAX_SAVE,conv1,conv2,CS_V
      
      
      IPRINT=IPRT
      NPRINT=IPRINT
      
      XF=XFINI
      FFF=0d0
	
      SQRT_GN=dsqrt(67.07d0)/1d+20
      PBNTOGEV=1d0/0.3893796623/1d+09	
       
      Pi=DACOS(-1d0)
      
      FFF=0
      NNOGOOD=0
      
      call vini
***************************

      call isachp

*****************************
! Z1,Z2,W1, e1, e2, e3, n1,n2,n3,u1,d1,c1,s1,t1,b1,gss
****************************
	PM(1)=abs(A(93))   	!neutr1
        PM(2)=abs(A(95))   	!neutr2
        PM(3)=abs(A(131))  	!gluino
	
        PM(4)=abs(A(89))   	!charg
        PM(5)=abs(A(132))  	!sel
        PM(6)=abs(A(134))  	!smu
        PM(7)=abs(A(106))  	!sta
        PM(8)=abs(A(122))   	!stop
        PM(9)=abs(A(126))  	!sbot
	
   	
        PM(10)=abs(A(136))  	!snu_e
        PM(11)=abs(A(137))  	!snu_l
        PM(12)=abs(A(112))  	!snu_tau
        PM(13)=abs(A(138))  	!sup1
        PM(14)=abs(A(140))  	!sd_1
        PM(15)=abs(A(142))  	!sc1
        PM(16)=abs(A(144))  	!ss_1
	
        PM(17)=	PM(4)		!charg
        PM(18)=	PM(5)  		!sel
        PM(19)=	PM(6)  		!smu
        PM(20)=	PM(7)  		!sta
        PM(21)=	PM(8)		!stop
        PM(22)=	PM(9)  		!sbot
	
 
        PM(23)=	PM(10)  	!snu_e
        PM(24)=	PM(11)  	!snu_l
        PM(25)=	PM(12)  	!snu_tau
        PM(26)=	PM(13)  	!sup1
        PM(27)=	PM(14)  	!sd_1
        PM(28)=	PM(15)  	!sc1
        PM(29)=	PM(17)  	!ss_1
***********************************************      

 
      SUPPEXP=1E+20
 
      DO IPM=2,IPTOT
        if(min(pm(IPM),PM(1)).ne.pm(1).and.iallow.eq.0) then
          iallow=-11 !NEUTR IS NOT THE LIGHTEST
          if(iprt.gt.2) print *,'Z1 IS NOT LSP',IPM,PM(IPM),PM(1)
          NNOGOOD=5
          goto 555
          ENDIF 
          SUPPEXP=min(SUPPEXP,pm(IPM)/pm(1))
      ENDDO

      geff=0d0
      DO IDEL=1,IPTOT
        del(idel)=(pm(idel)-pm(1))/pm(1)
        geff=geff+ 2. * (1.+del(IDEL))**(3./2.) * Exp(-del(idel)/XF)
        if(iprint.gt.2) 
     &print *,  (1.+del(IDEL))**(3./2.) * Exp(-del(idel)/XF), del(idel)
      ENDDO
      IF(IPRINT.gt.1) print *,'GEFF=', GEFF,'   ','SUPPEXP=',SUPPEXP
      NDIMUSER_EFF=NDIMUSER
      IF(SUPPEXP.ge.SUPPEXPMAX.and.NDIMUSER.eq.3)  NDIMUSER_EFF=2 
      if(iprint.gt.2) 
     & print *, 'NDIMUSER_EFF=', NDIMUSER_EFF,SUPPEXP
     

cccccccccccccccccc    
c      goto 111

      NST=0
      NDIMUSER_SAVE=NDIMUSER
      NDIMUSER_EFF_SAVE=NDIMUSER_EFF
c       print *,'ISARED:', NPROC_MIN,NPROC_MAX

      IF(NDIMUSER.ge.2) NDIMUSER=2
      NDIMUSER_EFF=NDIMUSER
      NCALL_SAV=NCALL
      NCALL=NCALL/2.
888   continue
      NST=NST+1
      IF(NST.GT.NST_MAX) goto 666
      
      IPRINT=IPRT
c--------------------	
      FFF_OLD=FFF  
      FFF=FUNC_INT(IPRINT)
      IF(FFF.lt.CS_MIN) GOTO  666
c--------------------	  

      IF(NDIMUSER_SAVE.eq.2) GOTO 777
      if(IPRT.gt.1) print *,'FFF0=', FFF

        
      IF(FFF.le.0) then
         FFF=0.
         goto 666
      ENDIF


       XFI=LOG(
     &     PM(1)/(2d0*Pi**3)*geff/2d0*sqrt(45d0/(2d0*81d0))/SQRT_GN
     &     *FFF*PBNTOGEV*SQRT(XF)
     &     )
       

       XF_OLD=XF
       XF=1D0/XFI
       


       IF(iprint.gt.1) then
       print *,'======================='
       print *,'XF = ',XF,1d0/XF
       print *,'CS  = ',FFF
       print *,'NST = ',NST
       print *,'======================='
       endif

       IF(XF.le.0.) then
         XF=XFINI
         FFF=0
         goto 666
       ENDIF

ccc       IF(abs(XF-XF_OLD)/XF.gt.0.01) goto 888
       IF(abs(FFF-FFF_OLD)/FFF.gt.0.03) goto 888
ccc    print *,'xxxxxxxx 3d Integration xxxxxxxxxx'
 777   CONTINUE
 666   NCALL=NCALL*2.
       NDIMUSER=NDIMUSER_SAVE
       NDIMUSER_EFF=NDIMUSER_EFF_SAVE
       IF(NDIMUSER.eq.2    ) goto 999
       IF(NDIMUSER_EFF.eq.2) NCALL=NCALL/2.

       IF(FFF.lt.CS_MIN) THEN
          FFF=1.E-20
          GOTO 999
       ENDIF


       FFF=FUNC_INT(IPRINT)
 999   continue
  
c      print *,'IPRINT=',IPRINT
       IF(IPRINT.ge.3) then
         NPROC_MAX_S=NPROC_MAX
         NPROC_MIN_S=NPROC_MIN
         SUMM=0.
         
         
         DO I =NPROC_MIN_S,NPROC_MAX_S,NPROC_STEP
           if(SUPPEXP.ge.SUPPEXPMAX.and.I.gt.26) GOTO 444
           NPROC_MIN=I
           NPROC_MAX=I
           FFF_TMP=FUNC_INT(-1)
           IF(FFF_TMP/FFF.gt.0.01) THEN
             print '(I6,A4,F6.2,A2,4A8)',
     &       I,'  ',FFF_TMP/FFF*100,' %',
     &       (PINF(I,J),J=1,4)
             SUMM=SUMM+FFF_TMP/FFF
           ENDIF
         ENDDO
 444     continue	 

         NPROC_MAX=NPROC_MAX_S
         NPROC_MIN=NPROC_MIN_S
       endif

      
 
 555   continue
       if(iprint.gt.1) then
       print *,'===========FINAL======='
       print *, 'freez out temp=',1d0/XF
       print *, 'n steps       =',NST
       print *, 'CS (fb)       =',FFF*1000d0
       print *, 'OMEGA H^2     =',OMEGA(FFF)
       print *,'======================='
       endif

       OMGH2 =OMEGA(FFF)
       SIGMA =FFF
       XFREEZ=XF
       NSTEPS =NST

       NCALL=NCALL_SAV
      
cccccccccccccccccccccccccccccccccccccccccccc
 111   continue
       NDIMUSER_SAVE=NDIMUSER
       AS_MAX_SAVE  =AS_MAX
       NDIMUSER=1
       NDIMUSER_EFF=1

       
       VVV=1.E-03
       AS_MAX=2.*sqrt(1.+(VVV/2.)**2)  
       CONV1 =  2.998E+10  ! speed of light cm/sec
       CONV2 =  1.000E-36  ! 1pb ==> cm^2

       CS_V  =  FUNC_INT(IPRINT)
       FFF_V =  CS_V*VVV*CONV1*CONV2

       NDIMUSER=NDIMUSER_SAVE
       AS_MAX=AS_MAX_SAVE
       

cccccccccccccccccccccccccccccccccccccccccccc       

      RETURN
      END

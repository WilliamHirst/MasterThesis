CDECK  ID>, ISAWIG.
C-----------------------------------------------------------------------
C
C                          HERWIG-ISAJET interface
C
C-----------------------------------------------------------------------
C
C                    VERSION 1.200 of 24th September 2002.
C
C-----------------------------------------------------------------------
C  Subroutines to produce an output file from ISAJET with SUSY masses 
C  and decay modes in a form HERWIG can read in.
C-----------------------------------------------------------------------
C  We also include the code to calculate the R-parity violating decay
C  modes which ISAJET does not include
C-----------------------------------------------------------------------
C  We now include an interface to HDECAY to allow the use of 
C  next-to-leading order Higgs decay rates
C-----------------------------------------------------------------------
C
C  We have made changes to allow the use of the code with either ISAJET
C  7.63/64 or 7.58 this is to enable the Snowmass accord points to be 
C  still be generated with the defined version of ISAJET.
C
C  The default is ISAJET7.64 in order to obtain code compatible with
C  7.63/7.58 the compiler option
C
C  -DISAJET758   for ISAJET7.58
C  -DISAJET763   for ISAJET7.63
C
C  should be used.
C
C-----------------------------------------------------------------------
C
C  We have also made changes to make using HDECAY easier and to allow 
C  the use of the new version of HDECAY.
C
C  The default is to compile a dummy subroutine. If you wish to use
C  HDECAY you should use the following compiler options
C
C -DHDECAY2   for version 2.0 of HDECAY
C -DHDECAY3   for version 3.0 of HDECAY
C
C-----------------------------------------------------------------------
C
C       F. Paige, March 2005:
C  Patchy format; Isajet common blocks included with +CDE statements.
C  OPEN statement removed.
C  Unit number for output passed as first argument, IH, and DATA
C  statement removed.
C  HDECAY version number (0,2,3) passed as second argument, IHDCY.
C  Ignored unless corresponding version selected with Patchy flag.
C
C-----------------------------------------------------------------------
C
C  This block contains the following code
C
C  ISAWIG to output the HERWIG decay table
C
C  Subroutines for R parity
C
C  RPDECY calculates all the decay rates
C  RPINF1 is a function for the amplitude squared terms in 3-body ME's
C  RPINF2 is a function for the interference terms in 3-body ME's
C  RPINT1 is the integrand for RPINF1
C  RPINT2 is the integrand for LR interference terms
C  RPINT3 is the integrand for the remaining terms
C  RPMAIN routine for the user to enter the couplings
C  RPNORM adds R-parity modes to tables, and removes modes < MINBR
C  RPMODA adds an R-parity violating mode to the table
C  RPRATE is a function for the 2-body rates
C  RPRTCH routine to test for negative square roots
C
C  HDECAY interface routines
C
C  HDCYAD is a routine used by the HDECAY interface to add decay modes
C         to the ISAJET decay tables
C  HDCYSY is the interface routine to HDECAY
C
C--13/04/99 Modified to work with ISAJET 7.42  by Peter Richardson 
C--02/07/01 Modified to work with AMSB models in ISAJET 7.48 and 7.51
C                       by Bryan Webber
C--09/04/02 Fixes to various RPV bugs and top decays to charged Higgs
C                       by Peter Richardson
C--19/09/02 Changes to allow use with either ISAJET7.58 or ISAJET7.63
C                       by Peter Richardson
C--19/09/02 Incorperated the HDECAY interface and made changes to allow
C           the use of HDECAY2 or 3 by Peter Richardson
C--24/09/02 Checked the code worked with ISAJET7.64 Peter Richardson
C
C-----------------------------------------------------------------------
C All the RPV rates from JHEP 0004:008,2000 
C                        by H. Dreiner, P. Richardson and M.H. Seymour
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CDECK  ID>, ISAWIG
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author : Bryan Webber, Kosuke Odagiri & Peter Richardson
C-----------------------------------------------------------------------       
      SUBROUTINE ISAWIG(IH,IHDCY,MT,M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                      MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR)
C-----------------------------------------------------------------------
C
C  Writes a table of SUSY particle (and top quark) properties
C  and decays in a format suitable for reading by HERWIG.
C  Call at end of SSRUN or SUGRUN.
C  MT = top mass used in ISASUSY or ISASUGRA.
C  Other parameters passed via common /SUGMG/ & /SSMODE/
C
C--18/08/98 modified by BRW to include gravitino (for GMSB option)
C--11/01/99 modified by KO to output mixing matrices in correct format
C--18/01/99 modified by KO to output correct sign of ALFAH
C--02/04/99 modified by PR to include R-parity violation
C--13/04/99 modified by PR to work with ISAJET 7.42
C--23/04/99 modified by PR to inculde KO A and mu terms
C--25/05/99 modified by PR to change R-parity ME code to 300
C--18/06/99 modified by PR to move data statements and rewrite a format
C                             statement caused problems with LINUX 
C--14/08/99 modified by PR to change file format and remove modes < MINBR
C--16/11/99 modified by PR to interface with HDECAY
C--31/03/00 modified by PR to fix a bug in the chargino mixing,
C                             the sfermion mixing angles
C--25/05/00 modified by BRW to include 1st and 2nd generation pseudoscalar
C                              mesons as possible decay products (for AMSB)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C          ISAJET common blocks and EQUIVALENCE
C
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
      REAL MSS(72)
      EQUIVALENCE (AMGLSS,MSS(1))
C
C  Common Block for the R-parity violating couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C--Inputs from ISAJET
      REAL MT,M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,
     &                        MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR
C
      REAL GEV2S,WIDTH,LTIM,MASS,THX,THY
      INTEGER IDEC(5),I,J,K,L,M,ID,JD,IH,NSS,NPP
      INTEGER IHDCY
C  SUSY HERWIG names, charge & PDG codes, ISAJET ID & mass codes
      PARAMETER (NSS=65,NPP=NSS+43)
      INTEGER IMISA(NSS),NDEC(NSS),IDHW(NPP),IDISA(NPP)
      DATA (IDISA(I),IMISA(I),I=1,NSS-2)/
     & 22, 4, 21, 2, 23, 6, 24, 8, 25,12, 26,16,-22, 4,-21, 2,
     &-23, 6,-24, 8,-25,12,-26,16, 42, 5, 41, 3, 43, 7, 44, 9,
     & 45,13, 46,17,-42, 5,-41, 3,-43, 7,-44, 9,-45,13,-46,17,
     & 32,18, 31,26, 34,20, 33,27, 36,24, 35,28,-32,18,-31,26,
     &-34,20,-33,27,-36,24,-35,28, 52,19, 51, 0, 54,21, 53, 0,
     & 56,25, 55, 0,-52,19,-51, 0,-54,21,-53, 0,-56,25,-55, 0,
     & 29, 1, 30,31, 40,32, 50,33, 60,34, 39,51, 49,52,-39,51,
     &-49,52, 82,55, 83,56, 84,57, 86,58,-86,58, 91,66/
C  non-SUSY HERWIG names & codes, PDG & ISAJET codes
      DATA (IDHW(I),IDISA(I),I=NSS-1,NPP)/
     &  6,  6, 12, -6,  1,  2,  2,  1,  3,  3,  4,  4,  5,  5,  7, -2,
     &  8, -1,  9, -3, 10, -4, 11, -5,121, 12,122, 11,123, 14,124, 13,
     &125, 16,126, 15,127,-12,128,-11,129,-14,130,-13,131,-16,132,-15,
     & 13,  9, 59, 10,198, 80,199,-80,200, 90,
     & 21, 110, 38, 120, 30,-120, 22, 220, 46, 130, 34,-130,
     & 50, 230, 42,-230, 25, 330,175, 140,140,-140,171, 240,
     &136,-240,179, 340,144,-340,163, 440/
C--constants etc.
      DATA GEV2S/6.582122E-25/
      DATA NDEC/NSS*0/
C Input variable
      CHARACTER *1 HYORN
C  Decide whether or not to use HDECAY
C  Valid version (2 or 3) must be defined with Patchy flag.
C  Calculate the R-parity violating modes
      CALL RPMAIN
C  Add R-parity violating modes if needed and remove all modes less MINBR
      CALL RPNORM
C  Tell HERWIG where to store SUSY particles
      DO I=1,57
        IDHW(I)=400+I
      ENDDO
      DO I=58,62
        IDHW(I)=145+I
      ENDDO
      IDHW(63)=458
C
C  Output SUSY particle + top quark table
C
      WRITE (IH,'(I4)') NSS
      DO I=1,NSS
        IF (I.LT.63) THEN
          MASS=0.
          IF (IMISA(I).NE.0) MASS=MSS(IMISA(I))
          IF (I.GT.49.AND.I.LT.58) MASS=-MASS
        ELSEIF (I.EQ.63) THEN
          MASS=0.
        ELSE
C t and tbar
          MASS=MT
        ENDIF
C  Compute lifetime in ps
        ID=ABS(IDISA(I))
        WIDTH=0.
        DO J=1,NSSMOD
          IF (ISSMOD(J).EQ.ID) THEN
            NDEC(I)=NDEC(I)+1
            WIDTH=WIDTH+GSSMOD(J)
          ENDIF
        ENDDO
        IF (WIDTH.NE.0.) THEN
          LTIM=GEV2S/WIDTH
        ELSE
          LTIM=1E30
        ENDIF
        WRITE (IH,1) IDHW(I),MASS,LTIM
      ENDDO
    1 FORMAT(I5,F12.4,E15.5)
C 
C  Output decay modes
C
      DO I=1,NSS
        WRITE (IH,'(I4)') NDEC(I)
        IF (NDEC(I).NE.0) THEN
          ID=ABS(IDISA(I))
          DO J=1,NSSMOD
            IF (ISSMOD(J).EQ.ID) THEN
C  Translate decay products
              DO K=1,5
                L=JSSMOD(K,J)
                IF (L.EQ.0) THEN
                  IDEC(K)=0
                ELSE
                  IF (IDISA(I).LT.0) L=-L
                  DO M=1,NPP
                    IF (IDISA(M).EQ.L) THEN
                      IDEC(K)=IDHW(M)
                      GO TO 10
                    ENDIF
                  ENDDO
C   Antiparticle=particle
                  L=-L
                  DO M=1,NPP
                    IF (IDISA(M).EQ.L) THEN
                      IDEC(K)=IDHW(M)
                      GO TO 10
                    ENDIF
                  ENDDO
                  PRINT *,' Unknown ISAJET ID =',L
                  IDEC(K)=20
   10             CONTINUE
                ENDIF
              ENDDO
              JD=IDHW(I)
              IF (JD.EQ.6.OR.JD.EQ.12) THEN
C--bug fix 9/04/02 by P.R. for t --> H b
C   Special for t and tbar
                IF (IDEC(1).GT.120.AND.IDEC(1).LT.133.AND.
     &              IDEC(3).NE.0) THEN
C   Leptonic decay: ISAJET order is wrong for M.E. calcn
                  WRITE (IH,11) JD,BSSMOD(J),100,IDEC(2),IDEC(1),
     &                  IDEC(3),0,0
                ELSEIF(IDEC(3).NE.0) THEN
C   Nonleptonic decay
                  WRITE (IH,11) JD,BSSMOD(J),100,(IDEC(K),K=1,5)
C--Higgs decay
                ELSE
                  WRITE (IH,11) JD,BSSMOD(J),0,(IDEC(K),K=1,5)
                ENDIF
C--end of fix
C   RPARITY 3-body matrix elements, all code 200, HERWIG decides what to do
              ELSEIF(JD.GE.449.AND.JD.LE.457.AND.IDEC(1).LE.140.
     &          AND.IDEC(2).LE.140.AND.IDEC(3).LE.140) THEN
                WRITE (IH,11) JD,BSSMOD(J),300,(IDEC(K),K=1,5)
              ELSE
                WRITE (IH,11) JD,BSSMOD(J),0,(IDEC(K),K=1,5)
              ENDIF
   11         FORMAT(I6,E16.8,6I6)
   12         FORMAT(1X,A8,'  -->   ',3(A8,2X),1P,E12.4)
            ENDIF
          ENDDO
        ENDIF
      ENDDO  
C  tan(beta) and alpha, neutralinos, charginos, and sfermion mixings
      WRITE (IH,'(2F16.8)') 1./MSS(30),-MSS(59)
C  neutralino mixing matrix - in ascending mass order (rows)
C             and in the order (bino, w3ino, higgs1, higgs2) (columns)
      WRITE (IH,13) MSS(38),MSS(37),-MSS(36),-MSS(35)
      WRITE (IH,13) MSS(42),MSS(41),-MSS(40),-MSS(39)
      WRITE (IH,13) MSS(46),MSS(45),-MSS(44),-MSS(43)
      WRITE (IH,13) MSS(50),MSS(49),-MSS(48),-MSS(47)
C--Bug fix 22/03/00 by Peter Richardson
      THX=SIGN(1.0D0,1.0D0/TAN(MSS(53)))
      THY=SIGN(1.0D0,1.0D0/TAN(MSS(54)))
c                   WMXVSS(1,1) WMXVSS(1,2)  WMXVSS(2,1) WMXVSS(2,2)
      WRITE (IH,13)     -SIN(MSS(54)),    -COS(MSS(54)),
     &              -THY*COS(MSS(54)),THY*SIN(MSS(54))
c                   WMXUSS(1,1) WMXUSS(1,2)  WMXUSS(2,1) WMXUSS(2,2)
      WRITE (IH,13)      -SIN(MSS(53)),  -COS(MSS(53)),
     &               -THX*COS(MSS(53)),THX*SIN(MSS(53))
C  sfermion mixing A terms and mu for 
      WRITE (IH,'(3F16.8)') -MSS(61),-MSS(63),-MSS(65)
      WRITE (IH,'(3F16.8)')  MSS(60),MSS(62),MSS(64)
      WRITE (IH,'( F16.8)') -MSS(29)
C  R-parity violating couplings
      WRITE(IH,'(L5)') RPARTY
      IF(.NOT.RPARTY) THEN
        WRITE(IH,20) (((LAMDA1(I,J,K),K=1,3),J=1,3),I=1,3)
        WRITE(IH,20) (((LAMDA2(I,J,K),K=1,3),J=1,3),I=1,3)
        WRITE(IH,20) (((LAMDA3(I,J,K),K=1,3),J=1,3),I=1,3)
      ENDIF
      RETURN
 13   FORMAT(4F16.8)
 19   FORMAT(I6)
 20   FORMAT(27E16.8)
      END

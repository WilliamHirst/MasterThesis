CDECK  ID>, SETHSS.
      SUBROUTINE SETHSS
C
C          Set the MSSM Higgs parameters in /HCON/.
C          HMASS  = Higgs mass for HTYPE
C          HGAM   = Higgs width
C          HGAMSS = Higgs partial widths. Note HGAMSS is not
C                   necessarily diagonal for SUSY decays.
C          ZSTARS = minimum allowed mass for Z*
C
C          Note LISTSS(78) => W+, LISTSS(79) => W-, LISTSS(80) => Z0
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
C          LISTSS IDENT and JETTYPE codes
C       ISGL  ISUPL -ISUPL  ISDNL -ISDNL  ISSTL -ISSTL  ISCHL -ISCHL
C          1      2      3      4      5      6      7      8      9
C      ISBT1 -ISBT1  ISTP1 -ISTP1  ISUPR -ISUPR  ISDNR -ISDNR  ISSTR
C         10     11     12     13     14     15     16     17     18
C     -ISSTR  ISCHR -ISCHR  ISBT2 -ISBT2  ISTP2 -ISTP2   ISW1  -ISW1
C         19     20     21     22     23     24     25     26     27
C       ISW2  -ISW2   ISZ1   ISZ2   ISZ3   ISZ4  ISNEL -ISNEL   ISEL
C         28     29     30     31     32     33     34     35     36
C      -ISEL  ISNML -ISNML  ISMUL -ISMUL  ISNTL -ISNTL ISTAU1-ISTAU1
C         37     38     39     40     41     42     43     44     45
C       ISER  -ISER  ISMUR -ISMUR ISTAU2-ISTAU2      9      1     -1
C         46     47     48     49     50     51     52     53     54
C          2     -2      3     -3      4     -4      5     -5      6
C         55     56     57     58     59     60     61     62     63
C         -6     11    -11     12    -12     13    -13     14    -14
C         64     65     66     67     68     69     70     71     72
C         15    -15     16    -16     10     80    -80     90   ISHL
C         73     74     75     76     77     78     79     80     81
C       ISHH   ISHA   ISHC  -ISHC
C         82     83     84     85
      COMMON/LISTSS/LISTSS(85)
      INTEGER LISTSS
      SAVE /LISTSS/
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
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
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
      DOUBLE PRECISION AQDP,BQDP,EZDP
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR
C
      REAL AMASS
      REAL AM12
      INTEGER I,J,N,IQ1,IQ2,IW,K
      INTEGER LISTJ(25),LISTW(4)
C
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
      DATA LISTW/10,80,-80,90/
C
C          Initialize
C
      IF(IHTYPE.EQ.0) THEN
        WRITE(ITLIS,*) ' YOU MUST SELECT AN HTYPE FOR SUSY HIGGS'
        WRITE(ITLIS,*) ' JOB TERMINATED'
        STOP99
      ENDIF
      HMASS=AMASS(IHTYPE)
      HGAM=0.
      DO 100 I=1,85
        DO 110 J=1,85
          HGAMSS(I,J)=0
110     CONTINUE
100   CONTINUE
C
C          Extract widths from SSMODE common block
C          Note the only 3-body modes are Zff or Wff
C          These are added to the ZZ and WW entries in HCONSS,
C          and the Z* or W* decay is generated later, as for SM Higgs
C
      DO 200 N=1,NSSMOD
        IF(ISSMOD(N).NE.IHTYPE) GO TO 200
        HGAM=HGAM+GSSMOD(N)
        IF(JSSMOD(3,N).NE.0) THEN
C          3-body modes
          IF(IABS(JSSMOD(1,N)).EQ.80) THEN
            HGAMSS(78,79)=HGAMSS(78,79)+0.5*GSSMOD(N)
            HGAMSS(79,78)=HGAMSS(79,78)+0.5*GSSMOD(N)
          ELSEIF(JSSMOD(1,N).EQ.90) THEN
            HGAMSS(80,80)=HGAMSS(80,80)+GSSMOD(N)
          ELSE
            WRITE(ITLIS,1000) ISSMOD(N),(JSSMOD(K,N),K=1,5)
1000        FORMAT(' SETHSS: UNEXPECTED MODE ',I8,' --> ',5I8)
            STOP 99
          ENDIF
          GO TO 200
        ELSE
C          2-body modes
          DO 210 I=1,85
            IF(JSSMOD(1,N).NE.LISTSS(I)) GO TO 210
            DO 220 J=1,85
              IF(JSSMOD(2,N).NE.LISTSS(J)) GO TO 220
              HGAMSS(I,J)=HGAMSS(I,J)+.5*GSSMOD(N)
              HGAMSS(J,I)=HGAMSS(J,I)+.5*GSSMOD(N)
              GO TO 200
220         CONTINUE
210       CONTINUE
        ENDIF
        WRITE(ITLIS,1000) ISSMOD(N),(JSSMOD(K,N),K=1,5)
        STOP99
200   CONTINUE
C
C          W* and Z* mass limits
C
      DO 300 I=1,2
        ZSTARS(1,I)=0.
        DO 310 IW=2,4
          ZSTARS(IW,I)=AMASS(LISTW(IW))
          DO 320 IQ1=2,25
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 320
            IF(GOWW(IQ1,I).AND.GOWW(IQ2,I)) THEN
              AM12=AMASS(LISTJ(IQ1))+AMASS(LISTJ(IQ2))+1.0
              ZSTARS(IW,I)=MIN(ZSTARS(IW,I),AM12)
            ENDIF
320       CONTINUE
310     CONTINUE
300   CONTINUE
      RETURN
      END

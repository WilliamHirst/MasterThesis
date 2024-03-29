CDECK  ID>, ISALHE.
      SUBROUTINE ISALHE
C
C     USING NOEVOL AND NOHADR, DECAY SUBPROCESS PARTICLES TO FILL
C     PARTCL COMMON BLOCK. THEN WRITE TO A .lhe FILE,
C     SO EVENT CAN BE PASSED TO OTHER GENERATORS FOR
C     SHOWERING, HADRONIZATION AND UNDERLYING EVENT
C
      IMPLICIT NONE
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
      INTEGER MXJETS
      PARAMETER (MXJETS=10)
      COMMON/PJETS/PJETS(5,MXJETS),IDJETS(MXJETS),QWJET(5),IDENTW 
     $,PPAIR(5,4),IDPAIR(4),JPAIR(4),NPAIR,IFRAME(MXJETS)
      SAVE /PJETS/
      INTEGER   IDJETS,IDENTW,IDPAIR,JPAIR,NPAIR,IFRAME
      REAL      PJETS,QWJET,PPAIR
      COMMON/PINITS/PINITS(5,2),IDINIT(2)
      SAVE /PINITS/
      INTEGER   IDINIT
      REAL      PINITS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
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
      COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
     $,XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      SAVE /XMSSM/
      REAL XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS
     $,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      LOGICAL GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
C
      INTEGER I,IFL1,IFL2,IP1,JET,NFIRST,IP
      INTEGER LISTJ(17),LISTW(4),LISTSM(30),IPAK,ID
      INTEGER ICOLOR(2,100),ISTAT,ITRANS,I1
      INTEGER IF1,IF2,IF3,JSPIN,INDEX,indx1,indx2,indx3,indx4
      INTEGER ND,N1,N2,N3,L1,IMO1,IMO2
      REAL AMASS
C
       INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
      DATA LISTSM/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,81/
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8/
      DATA LISTW/10,80,-80,90/
      DATA IPAK/100/
C     FILL PARTCL FROM JETPAR: FINAL PARTONS
        NPTCL=NJET
        DO 100 I=1,NJET
        PPTCL(1,I)=PT(I)*COS(PHI(I))
        PPTCL(2,I)=PT(I)*SIN(PHI(I))
        PPTCL(3,I)=P(I)*CTH(I)
        IF(KEYS(1)) THEN
          IDENT(I)=LISTJ(JETTYP(I))
        ELSEIF(KEYS(2)) THEN
          IDENT(I)=IDJETS(I)
        ELSEIF(KEYS(5).OR.(KEYS(10).AND.GOMSSM)) THEN
          IDENT(I)=LISTSS(JETTYP(I))
        ELSEIF(KEYS(6)) THEN
          IDENT(I)=LISTW(JETTYP(I))
        ELSEIF(KEYS(8)) THEN
          IF(JETTYP(1).LE.13) THEN
            IFL1=LISTJ(JETTYP(1))
          ELSE
            IFL1=10
          ENDIF
          IF(JETTYP(2).LE.13) THEN
            IFL2=LISTJ(JETTYP(2))
          ELSE
            IFL2=10
          ENDIF
          IDENT(1)=IFL1
          IDENT(2)=IFL2
        ELSEIF(KEYS(10)) THEN
          IDENT(I)=LISTSM(JETTYP(I))
        ENDIF
        PPTCL(5,I)=AMASS(IDENT(I))
        PPTCL(4,I)=SQRT(P(I)**2+PPTCL(5,I)**2)
        IORIG(I)=-(IPACK*I)
        IDCAY(I)=0
100   CONTINUE
C     Implement color connection for 2-> 2 subprocess
      IF (NPTCL.EQ.2) THEN
        CALL COLR22(IDINIT(1),IDINIT(2),IDENT(1),IDENT(2),ICOLOR)
      END IF
C     NOW DECAY FINAL STATE PARTONS
      DO 610 IP=1,NJET
        NFIRST=NPTCL+1
        JET=IP
        CALL DECAY(IP)
c
        ND=NPTCL-(NFIRST-1)
        N1=NFIRST
        N2=NFIRST+1
        L1=IP
        IF (ND.EQ.2) THEN
          CALL COLR12(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,ICOLOR)
        END IF
        IF (ND.EQ.3) THEN
          N3=NFIRST+2
          CALL COLR13(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,
     $                IDENT(N3),N3,ICOLOR)
        END IF
C
        DO 620 IP1=NFIRST,NPTCL
620     IORIG(IP1)=ISIGN(IABS(IORIG(IP1))+IPACK*JET,IORIG(IP1))
610   CONTINUE
C     NOW DECAY THE DECAY PRODUCTS
      IP=NJET+1
700   NFIRST=NPTCL+1
      JET=IABS(IORIG(IP))/IPACK
      CALL DECAY(IP)
c
        ND=NPTCL-(NFIRST-1)
        N1=NFIRST
        N2=NFIRST+1
        L1=IP
        IF (ND.EQ.2) THEN
          CALL COLR12(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,ICOLOR)
        END IF
        IF (ND.EQ.3) THEN
          N3=NFIRST+2
          CALL COLR13(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,
     $                IDENT(N3),N3,ICOLOR)
        END IF
C
      DO 720 IP1=NFIRST,NPTCL
720     IORIG(IP1)=ISIGN(IABS(IORIG(IP1))+IPACK*JET,IORIG(IP1))
      IP=IP+1
      IF (IP.LE.NPTCL) GO TO 700
C
C     Now output to isajet.lhe
      WRITE(LHEOUT,1001) 
C     Here one needs to invert particle IDs using LISTJ 
C     for idinit or LISTSS for IDENT if running SUSY
C     in order to match up with INOUT reaction code.
C
        ID=IPAK**3*JETTYP(2)+IPAK**2*JETTYP(1)+IPAK*INITYP(2)+INITYP(1)
c     If we have SUSY production, just dump out one type of subprocess, 
C     since Pythia can only handle 500 or less
        IF (GOMSSM) THEN
          ID=2160
        END IF
       WRITE(LHEOUT,1002) NPTCL+2,ID,1.,QSQ,ALFA,ALFQSQ
C     Write out initial state particles
        DO I=1,2
          WRITE(LHEOUT,1003) ITRANS(IDINIT(I),1),-1,0,0,
     $ICOLOR(1,I),ICOLOR(2,I),PINITS(1,I),PINITS(2,I),PINITS(3,I),
     $PINITS(4,I),PINITS(5,I),0.,9.
        END DO
        DO I=1,NPTCL
          IF (IDCAY(I).EQ.0.) THEN 
            ISTAT=1
          ELSE
            ISTAT=2
          END IF
          IF (IORIG(I).EQ.0) ISTAT=-1
          I1=IABS(IORIG(I))
          JET=I1/IPACK
          I1=I1-IPACK*JET
          I1=ISIGN(I1,IORIG(I))
          IF (I.LE.2) THEN
            IMO1=1
            IMO2=2
          ELSE
            IMO1=I1+2
            IMO2=0
          END IF
          WRITE(LHEOUT,1003) ITRANS(IDENT(I),1),ISTAT,IMO1,IMO2,
     $ICOLOR(1,I+2),ICOLOR(2,I+2),PPTCL(1,I),PPTCL(2,I),PPTCL(3,I),
     $PPTCL(4,I),PPTCL(5,I),0.,9.
        END DO
      WRITE(LHEOUT,1004) 
1001  FORMAT('<event>')
1002  FORMAT(4X,I3,4X,I8,3X,F12.5,3X,E12.6,3X,F12.6,3X,F12.6)
1003  FORMAT(6X,I8,3(2X,I4),2(2X,I3),5(2X,E12.6),2(1X,F2.0))
1004  FORMAT('</event>')
      RETURN
      END

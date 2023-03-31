CDECK  ID>, FLAVOR.
      SUBROUTINE FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
C
C          This subroutine unpacks the IDENT code ID=+/-IJKL
C
C          Mesons--
C          I=0, J<=K, +/- is sign for J
C          ID=110 for PI0, ID=220 for ETA, etc.
C
C          Baryons--
C          I<=J<=K in general
C          J<I<K for second state antisymmetric in (I,J), eg. L = 2130
C
C          Other--
C          ID=1,...,6 for quarks
C          ID=9 for gluon
C          ID=10 for photon
C          ID=11,...,16 for leptons
C          ID=20 for KS, ID=-20 for KL
C
C          I=21...26 for left scalar quarks
C          I=29 for gluino
C          I=30 for Z1SS
C          I=31...36 for left scalar leptons
C          I=39 for W1SS
C          I=40 for Z2SS
C          I=41...46 for right scalar quarks
C          I=49 for W2SS
C          I=50 for Z3SS
C          I=51...56 for right scalar leptons
C          I=60 for Z4SS
C
C          ID=80 for W+
C          ID=81,...,89 for Higgs mesons
C          ID=90 for Z0
C          ID=91 for gravitino
C          ID=92 for graviton
C
C          Incomplete meson multiplets used in b decays:
C          ID=10121     A1+(1260)
C          ID=10111     A10(1260)
C          ID=10131     K1+(1270)
C          ID=10231     K10(1270)
C          ID=30131     K1*+(1400)
C          ID=30231     K1*0(1400)
C          ID=132       K2*+(1430)
C          ID=232       K2*0(1430)
C          ID=10110     F0(980)     (mass = 1000 to allow K+K- decay)
C          ID=112       F2(1270)
C          ID=10441     PSI(2S)
C          ID=20440     CHI0
C          ID=20441     CHI1
C          ID=20442     CHI2
C
C          Diquarks--
C          ID=+/-IJ00, I<J for diquark composed of I,J.
C
C          INDEX is a sequence number used internally
C
C          Ver. 7.03: Make more robust by returning INDEX = 0 for
C          bad ID codes. Does not check for valid baryons, e.g.,
C          uuu with J = 1/2. Test on LABEL(1:3) = 'ERR' for this.
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QLMASS/AMLEP(100),NQLEP,NMES,NBARY
      SAVE /QLMASS/
      INTEGER   NQLEP,NMES,NBARY
      REAL      AMLEP
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
C
      INTEGER ID,IFL1,IFL2,IFL3,JSPIN,INDEX
      INTEGER I,J,K,IDABS,INDXSP
C
      IDABS=IABS(ID)
C
C          Select case
C
      IF(IDABS.GT.NQLEP-1.AND.IDABS.LT.80) GO TO 400
      IF(IDABS.GT.92.AND.IDABS.LE.100) GO TO 400
C          Quarks: ID < 100
      IF(IDABS.LT.100) GO TO 200
      I=MOD(IDABS/1000,10)
      J=MOD(IDABS/100,10)
      K=MOD(IDABS/10,10)
      JSPIN=MOD(IDABS,10)
C          Special hadrons
      IF(IDABS.GE.10000.OR.JSPIN.GT.1) GO TO 500
      IF(I.EQ.9.OR.J.EQ.9.OR.K.EQ.9) GO TO 400
C          Mesons: 100 < ID < 1000
      IF(IDABS.LT.1000) GO TO 100
C          Diquarks: ID > 1000 but K = 0
      IF(K.EQ.0.AND.JSPIN.EQ.0) GO TO 300
C
C          Baryons
C          Only X,Y baryons are QQX, QQY, Q=U,D,S.
C
      IF(I.GT.K.OR.J.GT.K.OR.J.EQ.0) GO TO 400
      IF(K.GT.6.AND.(I.GT.3.OR.J.GT.3)) GO TO 400
      IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,ID)
      IF(K.LE.6) THEN
        INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)+(K-1)*K*(2*K-1)/6
     1  +109*JSPIN+36*NMES+NQLEP+13
      ELSE
        INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)+9*(K-7)+91
     1  +109*JSPIN+36*NMES+NQLEP+13
      ENDIF
      RETURN
C
C          Mesons
C
100   CONTINUE
      IF(J.GT.K) GO TO 400
      IF(J.EQ.K.AND.ID.LT.0) GO TO 400
      IFL1=0
      IFL2=ISIGN(J,ID)
      IFL3=ISIGN(K,-ID)
      INDEX=J+K*(K-1)/2+36*JSPIN+NQLEP
      INDEX=INDEX+13
      RETURN
C
C          Quarks, leptons, etc
C
200   CONTINUE
      IFL1=0
      IFL2=0
      IFL3=0
      JSPIN=0
      INDEX=IDABS
      IF(IDABS.LT.20) RETURN
C          Define INDEX=20 for KS, INDEX=21 for KL
      INDEX=IDABS+1
      IF(ID.EQ.20) INDEX=20
C          INDEX=NQLEP+1,...,NQLEP+13 for W+, Higgs, Z0, GVSS, GRAV
      IF(IDABS.LT.80) RETURN
      INDEX=NQLEP+IDABS-79
      RETURN
C
C          Diquarks
C
300   IF(JSPIN.GT.0.OR.I.GT.J) GO TO 400
      IF(I.GT.6.OR.J.GT.6) GO TO 400
      IFL1=ISIGN(I,ID)
      IFL2=ISIGN(J,ID)
      IFL3=0
      JSPIN=0
      INDEX=109*NBARY+36*NMES+NQLEP+13+I+J*(J-1)/2
      RETURN
C
C          Error
C
400   CONTINUE
      IFL1=0
      IFL2=0
      IFL3=0
      JSPIN=0
      INDEX=0
      RETURN
C
C          Special mesons - used only for B decays
C
500   INDXSP=400
      IF(IDABS.EQ.10121) THEN
        INDEX=INDXSP+1
      ELSEIF(IDABS.EQ.10111) THEN
        INDEX=INDXSP+2
      ELSEIF(IDABS.EQ.10131) THEN
        INDEX=INDXSP+3
      ELSEIF(IDABS.EQ.10231) THEN
        INDEX=INDXSP+4
      ELSEIF(IDABS.EQ.30131) THEN
        INDEX=INDXSP+5
      ELSEIF(IDABS.EQ.30231) THEN
        INDEX=INDXSP+6
      ELSEIF(IDABS.EQ.132) THEN
        INDEX=INDXSP+7
      ELSEIF(IDABS.EQ.232) THEN
        INDEX=INDXSP+8
      ELSEIF(IDABS.EQ.10110) THEN
        INDEX=INDXSP+9
      ELSEIF(IDABS.EQ.112) THEN
        INDEX=INDXSP+10
      ELSEIF(IDABS.EQ.10441) THEN
        INDEX=INDXSP+11
      ELSEIF(IDABS.EQ.20440) THEN
        INDEX=INDXSP+12
      ELSEIF(IDABS.EQ.20441) THEN
        INDEX=INDXSP+13
      ELSEIF(IDABS.EQ.20442) THEN
        INDEX=INDXSP+14
      ELSEIF(IDABS.EQ.IDTAUL) THEN
        INDEX=INDXSP+15
      ELSEIF(IDABS.EQ.IDTAUR) THEN
        INDEX=INDXSP+16
      ELSE
        INDEX=0
      ENDIF
      IF(INDEX.GT.0) THEN
        IFL1=0
        IFL2=ISIGN(J,ID)
        IFL3=ISIGN(K,-ID)
      ELSE
        IFL1=0
        IFL2=0
        IFL3=0
      ENDIF
C
      RETURN
      END

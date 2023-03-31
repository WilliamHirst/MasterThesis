CDECK  ID>, AMASS.
      FUNCTION AMASS(ID)
C
C          Returns the mass of the particle with IDENT code ID.
C          Quark-based IDENT code.
C          Ver 7.10: Update masses and split B baryon degeneracy.
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
      INTEGER ID
      REAL AMASS
      REAL AMMES0(10),AMMES1(10),AMBAR0(30),AMBAR1(30)
      INTEGER IFL1,IFL2,IFL3,JSPIN,INDEX,IFL1A,IFL2A,IFL3A,IDA
C
C          0- meson mass table
C          pi0, pi+, eta, k+, k0, etap, ad0, d-, ds-, etac
C
      DATA AMMES0/.13496,.13957,.54745,.49367,.49767,.95775,1.8645
     $,1.8693,1.9688,2.9788/
C
C          1- meson mass table
C          rho0, rho+, omega, k*+, k*0, phi, ad*0, d*-, d*s-, jpsi
C
      DATA AMMES1/.7681,.7681,.78195,.89159,.89610,1.0194,2.0071
     $,2.0101,2.1103,3.0969/
C
C          1/2+ baryon mass table
C          x,p,n,-,-,s+,s0,s-,l,xi0,xi-,x,x,x
C          sc++,sc+,sc0,lc+,usc.,dsc.,ssc.,sdc.,suc.,ucc.,dcc.,scc.
C
      DATA AMBAR0/-1.,.93828,.93957,2*-1.,1.1894,1.1925,1.1974
     $,1.1156,1.3149,1.3213,3*-1.,2.4527,2.4529,2.4525,2.2849
     $,2.50,2.50,2.60,2.40,2.40,3.55,3.55,3.70,4*-1./
C
C          3/2+ baryon mass table
C          dl++,dl+,dl0,dl-,-,s*+,s*0,s*-,x,xi*0,xi*-,om-,x,x
C          uuc*,udc*,ddc*,x,usc*,dsc*,ssc*,x,x,,ucc*,dcc*,scc*,ccc*
C
      DATA AMBAR1/1.232,1.232,1.232,1.232,-1.,1.3823,1.3820
     $,1.3875,-1.,1.5318,1.5350,1.6722,2*-1.
     $,2.63,2.63,2.63,-1.,2.70,2.70,2.80,2*-1.,3.75,3.75
     $,3.90,4.80,3*-1./
C
C          Entry
C
      AMASS=-1.
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IDA=IABS(ID)
      IFL1A=IABS(IFL1)
      IFL2A=IABS(IFL2)
      IFL3A=IABS(IFL3)
      IF(IDA.GT.10000.OR.JSPIN.GT.1) GO TO 500
C
C          Diquarks
C
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) THEN
        AMASS=AMLEP(IFL1A)+AMLEP(IFL2A)
C
C          b and t particles. Only a few b masses are known, but we
C          guess a few others to make sure decays are allowed:
C
      ELSEIF(IFL3A.GT.4) THEN
        IF(IDA.EQ.150.OR.IDA.EQ.250) THEN
          AMASS=5.2786
        ELSEIF(IDA.EQ.151.OR.IDA.EQ.251) THEN
          AMASS=5.3246
        ELSEIF(IDA.EQ.350) THEN
          AMASS=5.3693
        ELSEIF(IDA.EQ.351) THEN
          AMASS=5.3693+0.04
        ELSEIF(IDA.EQ.2150) THEN
          AMASS=5.641
        ELSEIF(IDA.EQ.1150.OR.IDA.EQ.1250.OR.IDA.EQ.2250) THEN
          AMASS=5.641+0.171
        ELSEIF(IDA.EQ.2151) THEN
          AMASS=5.641+.04
        ELSEIF(IDA.EQ.1151.OR.IDA.EQ.1251.OR.IDA.EQ.2251) THEN
          AMASS=5.641+0.171+0.04
        ELSE
          AMASS=AMLEP(IFL2A)+AMLEP(IFL3A)-.03+.04*JSPIN
          IF(IFL1.NE.0) AMASS=AMASS+AMLEP(IFL1A)
        ENDIF
C
C          Quarks and leptons
C
      ELSEIF(IFL2.EQ.0) THEN
        AMASS=AMLEP(INDEX)
C
C          Mesons
C
      ELSEIF(IFL1.EQ.0) THEN
        INDEX=INDEX-36*JSPIN-NQLEP
        INDEX=INDEX-13
        AMASS=(1-JSPIN)*AMMES0(INDEX)+JSPIN*AMMES1(INDEX)
C
C          Baryons
C
      ELSE
        INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
        INDEX=INDEX-13
        AMASS=(1-JSPIN)*AMBAR0(INDEX)+JSPIN*AMBAR1(INDEX)
      ENDIF
      RETURN
C
C          Special hadrons - used only in B decays
C
500   IF(IDA.EQ.10121.OR.IDA.EQ.10111) THEN
        AMASS=1.230
      ELSEIF(IDA.EQ.10131.OR.IDA.EQ.10231) THEN
        AMASS=1.273
      ELSEIF(IDA.EQ.30131.OR.IDA.EQ.30231) THEN
        AMASS=1.412
      ELSEIF(IDA.EQ.132) THEN
        AMASS=1.4254
      ELSEIF(IDA.EQ.232) THEN
        AMASS=1.4324
      ELSEIF(IDA.EQ.10110) THEN
        AMASS=0.980+0.020
      ELSEIF(IDA.EQ.112) THEN
        AMASS=1.275
      ELSEIF(IDA.EQ.10441) THEN
        AMASS=3.686
      ELSEIF(IDA.EQ.20440) THEN
        AMASS=3.4151
      ELSEIF(IDA.EQ.20441) THEN
        AMASS=3.51053
      ELSEIF(IDA.EQ.20442) THEN
        AMASS=3.56617
      ELSEIF(IDA.EQ.IDTAUL.OR.IDA.EQ.IDTAUR) THEN
        AMASS=AMLEP(16)
      ELSE
        AMASS=0
      ENDIF
      RETURN
      END

CDECK  ID>, ZZSTAR.
      FUNCTION ZZSTAR(HM,IW)
C
C          Generate W* or Z* mass for H -> W W* or H -> Z Z* decay,
C          including the W or Z width in the propagator.
C          Ref: Marciano and Sirlin, Phys. Rev. D30, 248 (1984).
C
C          HM = generated Higgs mass, i.e. QMW**2
C          IW = 2   3   4
C               W+  W-  Z0
C
      IMPLICIT NONE
C
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
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
C
      REAL    HM,WM,WG,ZZSTAR,EPS,FBAR,R1,R2,RANF,X,F,DELTA,R,XM1
      INTEGER I,IW
C          WM and WG are the W or Z mass and width
      WM=WMASS(IW)
      WG=WGAM(IW)
      EPS=WM/HM
      DELTA=WM*WG/HM**2
C          FBAR is maximum of F below
      FBAR=12.*EPS**2*(1.-EPS)**2*(1.-EPS**2)
      R1=(2.*EPS-1.)/DELTA
      R2=EPS**2/DELTA
      R1=ATAN(R1)
      R2=ATAN(R2)
C          Generate Breit-Wigner and test remainder F against FBAR
      DO 100 I=1,NTRIES
        R=R1-RANF()*(R1-R2)
        XM1=DELTA*TAN(R)
        X=XM1+1.
        F=SQRT((X-2.*EPS)*(X+2.*EPS))
     $  *(X**2-12.*EPS**2*X+8.*EPS**2+12.*EPS**4)
        XM1=SQRT(XM1)
        ZZSTAR=HM*SQRT((EPS-XM1)*(EPS+XM1))
        IF(F.GT.FBAR*RANF()) RETURN
100   CONTINUE
C
      WRITE(ITLIS,9999) NTRIES
9999  FORMAT(' ERROR IN ZZSTAR ... NO MASS FOUND')
      STOP 99
      END

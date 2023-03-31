CDECK  ID>, SSMQCD.
      DOUBLE PRECISION FUNCTION SSMQCD(DM,DQ)
C-----------------------------------------------------------------------
C     Calculate leading-log running mass for quark with mass DM at 
C     scale Q, using alpha_s which is continuous across thresholds.
C     See Drees and Hikasa, Phys. Lett. B240: 455-464, Eq. 4.5.
C
C     Note the threshold is at Q = 2 m, not at Q = m as in MSbar.
C
C     Bisset's QCDRAD, WDHFFC
C     I think there was an error in the sense that some thresholds
C     were mixed up between m and 2m. I am changing everything to m
C     Modified Nov. 22, 2005, by XT.

C-----------------------------------------------------------------------
      IMPLICIT NONE
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C
      DOUBLE PRECISION DM,DQ,DLAM4,DLAM5,DLAM6,DNEFF,POW,RENORM
     $,DQBT,DQTP
C
C         Do nothing for light quarks
C
      IF(DM.LT.1.0) THEN
        SSMQCD=DM
        RETURN
      ENDIF
C
C          Calculate running mass
C
      DLAM4=DBLE(ALQCD4)
      DQBT=DBLE(AMBT)
      DQTP=DBLE(AMTP)
      SSMQCD=0
C          Q <  m(b)
      DNEFF=4
      POW=12.D0/(33.D0-2.*DNEFF)
      IF(DQ.LT.DQBT) THEN
        RENORM=(LOG(DM/DLAM4)/LOG(DQ/DLAM4))**POW
        SSMQCD=RENORM*DM
        RETURN
      ELSE
        RENORM=(LOG(DM/DLAM4)/LOG(DQBT/DLAM4))**POW
      ENDIF
C          m(b) < Q <  m(t)
      DNEFF=5
      POW=12.D0/(33.D0-2.*DNEFF)
      DLAM5=DEXP((25.D0*LOG(DLAM4)-LOG(DQBT**2))/23.D0)
      IF(DQ.GE.DQBT.AND.DQ.LT.DQTP) THEN
        RENORM=RENORM
     $  *(LOG(DQBT/DLAM5)/LOG(DQ/DLAM5))**POW
        SSMQCD=RENORM*DM
        RETURN
      ELSE
        RENORM=RENORM
     $  *(LOG(DQBT/DLAM5)/LOG(DQTP/DLAM5))**POW
      ENDIF
C           m(t) < Q
      DNEFF=6
      POW=12.D0/(33.D0-2.*DNEFF)
      DLAM6=DEXP((25.D0*LOG(DLAM4)-LOG(DQBT**2)
     $-LOG(4*AMTP**2))/21.D0) 
      RENORM=RENORM
     $*(LOG(DQTP/DLAM6)/LOG(DQ/DLAM6))**POW
      SSMQCD=RENORM*DM
      RETURN
      END

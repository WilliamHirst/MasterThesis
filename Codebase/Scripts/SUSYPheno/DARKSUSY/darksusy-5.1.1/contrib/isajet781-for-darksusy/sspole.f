CDECK  ID>, SSPOLE.
      REAL FUNCTION SSPOLE(MGMS,MUSQ,AS)
C*********************************************************************
C* Computes the on-shell (pole) gluino mass for given running (MSbar)*
C* gluino mass, defined at scale MUSQ, and given alpha_s (AS). The   *
C* squark masses are stored in the SQUARK COMMON block.              *
C* This function needs the complex functions B0 and B1.              *
C* Contributed by M. Drees; modified by H. Baer                      *
C  B0 contributions from Pierce et al.                               *
C  included on 9/23/05 by J.Ferrandis and H. Baer                    *
C                                                                    *
C  Version 7.30: Cast COMPLEX*16 to REAL*8 in standard way. :-(      *
C*********************************************************************
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
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
      REAL MGMS,MUSQ,AS,MGSQ,FAC
      DOUBLE PRECISION DMUSQ,DFAC
      COMPLEX*16 SSB1,SSB0
      DMUSQ=MUSQ
      XLAM = LOG(DMUSQ)
      MGSQ = MGMS*MGMS
C
C          Cast COMPLEX*16 to REAL*8:
C
      DFAC = -SSB1(MGSQ,0.,AMULSS) -SSB1(MGSQ,0.,AMURSS)
     $-SSB1(MGSQ,0.,AMDLSS) -SSB1(MGSQ,0.,AMDRSS)
     $-SSB1(MGSQ,0.,AMSLSS) -SSB1(MGSQ,0.,AMSRSS)
     $-SSB1(MGSQ,0.,AMCLSS) -SSB1(MGSQ,0.,AMCRSS)-
     $ ( SSB1(MGSQ,AMTP,AMT1SS)+ SSB1(MGSQ,AMTP,AMT2SS)+
     $  SSB1(MGSQ,4.0,AMB1SS) + SSB1(MGSQ,4.0,AMB2SS) )
     $ - AMTP*SIN(2.*THETAT)*(SSB0(MGSQ,AMTP,AMT1SS)-
     $   SSB0(MGSQ,AMTP,AMT2SS))/MGMS
     $ - 4.0*SIN(2.*THETAB)*(SSB0(MGSQ,4.0,AMB1SS)-
     $   SSB0(MGSQ,4.0,AMB2SS))/MGMS
      DFAC = DFAC + 15.D0 + 9.D0*LOG(DMUSQ/MGSQ)
      FAC=DFAC
      SSPOLE = MGMS*(1.0 + .0796*AS*FAC )
      RETURN
      END

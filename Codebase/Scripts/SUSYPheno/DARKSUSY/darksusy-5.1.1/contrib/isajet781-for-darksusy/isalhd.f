CDECK  ID>, ISALHD.
C--------------------------------------------------------------------
      SUBROUTINE ISALHD(LOUT,ID,J,JMAX)
C--------------------------------------------------------------------
C     
C     Output decay modes for ID in 'Les Houches accord 3' (LHA3) format
C     C. Balazs, May 21 2005, v0.1
C
C     Note these need not be contiguous, so the loop is over all modes 
C     in /SSMODE/.

CsB   ISAJET common blocks from SSPRT ...
CsB   ... explicitly included from v7.71
      IMPLICIT NONE
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
CsB   End of ISAJET common blocks
CsB   Local ISAJET related variables
      INTEGER ID,J,JMAX, I,K,NOUT,LOUT
      CHARACTER*5 SSID,LBLIN,LBLOUT(4)
      CHARACTER*40 VERSN,VISAJE
C
CsB   Local LHA3 related variables
      Integer PDGID,PDGIN,PDGOUT(4),iCnt
      Real Width
C
      If (J.EQ.1) then
C      	LOUT = 92
CsB     Open output file
      	Open(LOUT,FILE='ISALHD.out',FORM='FORMATTED')
CsB     Write header
        WRITE(LOUT,7000) 
     .  ' ISAJET decay tables in SUSY Les Houches accord format'
        WRITE(LOUT,7000)
     .  ' Created by ISALHD. Last revision: C. Balazs, 2005 May 25'
        VERSN=VISAJE()
        VERSN=VERSN(14:)
        WRITE(LOUT,7001)    'DCINFO', 
     ,                      'Program information'
        WRITE(LOUT,7012) 1, 'ISASUGRA from ISAJET       ',
     ,                      'Spectrum Calculator'
        WRITE(LOUT,7012) 2,  VERSN, 
     ,                      'Version number'

      End If
C
      Width=0.d0
      DO 90 I=1,NSSMOD
        IF(ISSMOD(I).NE.ID) GO TO 90
        LBLIN=SSID(ISSMOD(I))
        Width = Width + GSSMOD(I)
 90   CONTINUE
C
      iCnt=0
      NOUT=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).NE.ID) GO TO 100
        NOUT=NOUT+1
        LBLIN=SSID(ISSMOD(I))
        PDGIN=PDGID(LOUT,ISSMOD(I))
        If (iCnt.Eq.0) then
          Write(LOUT,'(A)') '#         PDG         Width'
          Write(LOUT,7500) PDGIN,Width,LBLIN//' decays'
          Write(LOUT,'(A)') '#          BR          NDA       ID1       
     .ID2       ID3       ID4'
          iCnt=1
        End If
        DO 110 K=1,4
        PDGOUT(K)=PDGID(LOUT,JSSMOD(K,I))
110     LBLOUT(K)=SSID(JSSMOD(K,I))
        If (PDGOUT(4).Eq.0) then
          If (PDGOUT(3).Eq.0) then
            WRITE(LOUT,7502) BSSMOD(I),2,
     ,                       PDGOUT(1),PDGOUT(2),
     ,                       LBLIN,(LBLOUT(K),K=1,4)
          Else
            WRITE(LOUT,7503) BSSMOD(I),3,
     ,                       PDGOUT(1),PDGOUT(2),PDGOUT(3),
     ,                       LBLIN,(LBLOUT(K),K=1,4)
          End If
        Else
            WRITE(LOUT,7504) BSSMOD(I),4,
     ,                       PDGOUT(1),PDGOUT(2),PDGOUT(3),PDGOUT(4),
     ,                       LBLIN,(LBLOUT(K),K=1,4)
        End If
100   CONTINUE
C
      If (J.EQ.JMAX) Close(LOUT)
C
CsB LHA3 format statements
C
 7000 FORMAT('# ',A)
C     Format to use for block statements
 7001 FORMAT('Block',1x,A,27x,'#',1x,A)
C     Indexed Char(12)
 7012 FORMAT(1x,I5,3x,A27,3x,'#',1x,A)
C     Write Decay Table
 7500 FORMAT('DECAY',1x,I9,1P,E16.8,0P,3x,'#',1x,A)
 7502 FORMAT(4x,1P,E16.8, 0P,3x,I2, 
     ,  0P,4x,I9, 0P,1x,I9,  23x,'#',1x,
     ,  A5,'  -->  ',4(A5,2X))
 7503 FORMAT(4x,1P,E16.8, 0P,3x,I2, 
     ,  0P,4x,I9, 0P,1x,I9, 0P,1x,I9,  13x,'#',1x,
     ,  A5,'  -->  ',4(A5,2X))
 7504 FORMAT(4x,1P,E16.8, 0P,3x,I2, 
     ,  0P,4x,I9, 0P,1x,I9, 0P,1x,I9, 0P,1x,I9,  3x,'#',1x,
     ,  A5,'  -->  ',4(A5,2X))
      RETURN
      END

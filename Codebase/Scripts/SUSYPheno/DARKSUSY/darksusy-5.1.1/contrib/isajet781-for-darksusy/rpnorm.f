CDECK  ID>, RPNORM
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson 
C----------------------------------------------------------------------- 
      SUBROUTINE RPNORM  
C----------------------------------------------------------------------- 
C     SUBROUTINE TO REMOVE ALL MODES WITH BRANCHING RATIO LESS THAN 1E-5
C     AND TO ADD THE RPARITY VIOLATING MODES IF NEEDED
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C--Common block containing the SUSY parameters
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
C
C--Common Block to contain R-parity violating decay rates
C
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      REAL GSSMD2,BSSMD2
      INTEGER ISSMD2,JSSMD2,NSSMD2
      SAVE /RPARRT/
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C  Common Block for the R-parity violating couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C--Local varaibles
      INTEGER SUSYLP,SUSYMN(9),SUSYMX(9),SUSYSP(9),I,J,NRMMDS,NRMMDR,
     &        ANUMRM,L,K
      REAL RATE,ZERO,RRATE,RPBR,EPS,RMBRRT,MINBR
      PARAMETER(ZERO=0.,MINBR=1E-5,EPS=1E-50)  
      DATA SUSYMN /31,51,21,41,29,82,86,30,39/
      DATA SUSYMX /36,56,26,46,29,84,86,60,49/
      DATA SUSYSP / 1, 1, 1, 1, 1, 1, 1,10,10/
C--Now we need to recalculate rates,etc.
      DO SUSYLP = 1,9
        DO 50 I=SUSYMN(SUSYLP),SUSYMX(SUSYLP),SUSYSP(SUSYLP) 
C--If Rparity conserved just remove modes less then MINBR
          IF(RPARTY) THEN 
            RMBRRT = ZERO
            NRMMDS = 0
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I.AND.BSSMOD(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMOD(J)
                ISSMOD(J) = 100000
                NRMMDS = NRMMDS+1 
              ENDIF
            ENDDO
C--Now remove the modes we have removed, and renormalise BR's to 1
            IF(NRMMDS.GT.0.OR.NRMMDR.GT.0) THEN
              ANUMRM = 0
              RMBRRT = 1/(1-RMBRRT)
              DO J=1,NSSMOD
 10             IF(ISSMOD(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                ISSMOD(J) = ISSMOD(J+ANUMRM)
                DO L=1,5
                  JSSMOD(L,J)=JSSMOD(L,J+ANUMRM)
                ENDDO
                GSSMOD(J)=GSSMOD(J+ANUMRM)
                BSSMOD(J)=BSSMOD(J+ANUMRM)
                IF(ISSMOD(J).EQ.100000) GOTO 10
                IF(ISSMOD(J).EQ.I) THEN
                  BSSMOD(J) = BSSMOD(J)*RMBRRT
                  GSSMOD(J) = GSSMOD(J)*RMBRRT
                ENDIF
              ENDDO
              ANUMRM = 0
              NSSMOD = NSSMOD-NRMMDS
            ENDIF
          ELSE
C--If R-paroty violated need to calculate rate and renormalise
C--Obtain the MSSM rate
            RATE  = ZERO
            RRATE = ZERO
            RPBR  = ZERO
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I) THEN
                RATE = RATE + GSSMOD(J)
              ENDIF
            ENDDO
C--Calculate the Rparity violating rate      
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I) THEN
                RRATE = RRATE + GSSMD2(J)
              ENDIF
            ENDDO        
            IF(RRATE.LT.EPS) GOTO 50
            RPBR = RATE/(RATE+RRATE)
            RATE = RATE+RRATE
C--Reset MSSM rates and branching ratios
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I) BSSMOD(J) = BSSMOD(J)*RPBR
            ENDDO
C--Calculate Rparity violating branching ratios
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I) BSSMD2(J) = GSSMD2(J)/RATE         
            ENDDO
C--Now remove any modes of the particle with branching ratio
C--less than MINBR
            RMBRRT = ZERO
            NRMMDS = 0
            NRMMDR = 0
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I.AND.BSSMOD(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMOD(J)
                ISSMOD(J) = 100000
                NRMMDS = NRMMDS+1 
              ENDIF
            ENDDO
            IF(RRATE.LT.EPS) GOTO 50
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I.AND.BSSMD2(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMD2(J)
                ISSMD2(J) = 100000
                NRMMDR = NRMMDR +1
              ENDIF
            ENDDO  
C--Now remove the modes we have removed, and renormalise BR's to 1
            IF(NRMMDS.GT.0.OR.NRMMDR.GT.0) THEN
              ANUMRM = 0
              RMBRRT = 1/(1-RMBRRT)
              DO J=1,NSSMOD
 20             IF(ISSMOD(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                ISSMOD(J) = ISSMOD(J+ANUMRM)
                DO L=1,5
                  JSSMOD(L,J)=JSSMOD(L,J+ANUMRM)
                ENDDO
                GSSMOD(J)=GSSMOD(J+ANUMRM)
                BSSMOD(J)=BSSMOD(J+ANUMRM)
                IF(ISSMOD(J).EQ.100000) GOTO 20
                IF(ISSMOD(J).EQ.I) THEN
                  BSSMOD(J) = BSSMOD(J)*RMBRRT
                  GSSMOD(J) = GSSMOD(J)*RMBRRT
                ENDIF
              ENDDO
              ANUMRM = 0
              DO J=1,NSSMD2
 30             IF(ISSMD2(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                  ISSMD2(J) = ISSMD2(J+ANUMRM)
                  DO L=1,5
                    JSSMD2(L,J)=JSSMD2(L,J+ANUMRM)
                  ENDDO
                  GSSMD2(J)=GSSMD2(J+ANUMRM)
                  BSSMD2(J)=BSSMD2(J+ANUMRM)
                  IF(ISSMD2(J).EQ.100000) GOTO 30
                  IF(ISSMD2(J).EQ.I) THEN
                  BSSMD2(J) = BSSMD2(J)*RMBRRT
                  GSSMD2(J) = GSSMD2(J)*RMBRRT
               ENDIF
              ENDDO
              NSSMOD = NSSMOD-NRMMDS
              NSSMD2 = NSSMD2-NRMMDR
            ENDIF
          ENDIF
 50     CONTINUE
      ENDDO
      IF(RPARTY) RETURN
C--NOW WE NEED TO COMBINE THE RPARITY AND MSSM MODES
      IF((NSSMOD+NSSMD2).GT.MXSS) THEN      
        WRITE(LOUT,*) 'WARNING REMOVING',NSSMOD+NSSMD2-MXSS,'MODES'
        print *,'WARNING EXCEEDS ISAJET NO OF MODES'
        print *,'REMOVING THE',NSSMOD+NSSMD2-MXSS,
     &          'MODES WITH LOWEST BRANCHING RATIO'
        print *,'RECOMMEND YOU RERUN WITH HIGHER MXSS'
      ENDIF
 110   IF((NSSMOD+NSSMD2).GT.MXSS) THEN
        RMBRRT = 1
        NRMMDS = 0
        DO J=1,NSSMOD
          IF(BSSMOD(J).LT.RMBRRT) THEN
            RMBRRT = BSSMOD(J)
            NRMMDS = J
          ENDIF
        ENDDO
        DO J=1,NSSMD2
          IF(BSSMD2(J).LT.RMBRRT) THEN
            RMBRRT = BSSMD2(J)
            NRMMDS = J+MXSS
          ENDIF
        ENDDO
C--remove mode with lowest branching ratio and rescale BR's
        RMBRRT = 1/(1-RMBRRT)
        IF(NRMMDS.LE.MXSS) NRMMDR = ISSMOD(NRMMDS)
        IF(NRMMDS.GT.MXSS) NRMMDR = ISSMD2(NRMMDS-MXSS)
        DO J=1,NSSMOD
          IF(ISSMOD(J).EQ.NRMMDR) THEN
            GSSMOD(J) = RMBRRT*GSSMOD(J)
            BSSMOD(J) = RMBRRT*BSSMOD(J)
          ENDIF
        ENDDO
        DO J=1,NSSMD2
          IF(ISSMD2(J).EQ.NRMMDR) THEN
            GSSMD2(J) = RMBRRT*GSSMD2(J)
            BSSMD2(J) = RMBRRT*BSSMD2(J)
          ENDIF
        ENDDO
        IF(NRMMDS.LE.MXSS) THEN
          DO J=NRMMDS,NSSMOD-1
            ISSMOD(J) = ISSMOD(J+1)
            DO L=1,5
              JSSMOD(L,J) = JSSMOD(L,J+1)
            ENDDO
            GSSMOD(J)=GSSMOD(J+1)
            BSSMOD(J)=BSSMOD(J+1)
          ENDDO
          NSSMOD = NSSMOD-1
        ELSE
          DO J=(NRMMDS-MXSS),NSSMD2-1
            ISSMD2(J) = ISSMD2(J+1)
            DO L=1,5
              JSSMD2(L,J) = JSSMD2(L,J+1)
            ENDDO
            GSSMD2(J) = GSSMD2(J+1)
            BSSMD2(J) = BSSMD2(J+1)
          ENDDO
          NSSMD2 = NSSMD2-1
        ENDIF   
      ENDIF   
      IF((NSSMD2+NSSMOD).GT.MXSS) GOTO 110
C--Now less than maximum so add the R parity violating modes
      DO J=1,NSSMD2
        NSSMOD=NSSMOD+1
        ISSMOD(NSSMOD) = ISSMD2(J)
        DO K=1,5
          JSSMOD(K,NSSMOD) = JSSMD2(K,J)
        ENDDO
        GSSMOD(NSSMOD) = GSSMD2(J)
        BSSMOD(NSSMOD) = BSSMD2(J)
        MSSMOD(NSSMOD) = 0
      ENDDO
      END

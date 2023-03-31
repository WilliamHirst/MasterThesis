CDECK  ID>, ST3MAT.
      FUNCTION ST3MAT(N,X)
C
C Purpose: Contains the equations for the SQUARED matrix element taken
C         mainly from Porod and Worhmann, PRD 55 2907 (1997).
C         There is an erratum to the paper at PRD 67 059902 (2003)
C
C         Requires: couplings, masses, widths
C
C         Passed: N - required by integration routine, but never
C                     used in the following function
C                 X - two values which run from 0 to 1 and are
C                     substituted for E_b and E_W
C
C         NB: The mt used in this function is the pole mass.
C
C         All references to the stop1 and lightest neutralino are
C         written as single entry arrays. This may make it easier
C         to generalise the function at a later date...
C
C         8/6/6
C
      IMPLICIT NONE
C
      COMMON/ST3MAS/MSB,MST,MCH,MNEU,MW,MB,MT
      DOUBLE PRECISION MSB(2),MST(1),MCH(2),MNEU(1),MW,MB,MT
      SAVE/ST3MAS/
C
      COMMON/ST3WID/SBWIDTH,CHWIDTH,TWIDTH
      DOUBLE PRECISION SBWIDTH(2),CHWIDTH(2),TWIDTH
      SAVE/ST3WID/
C
      COMMON/ST3COUP/CA1,CB1,CA2,CB2,CC2,CD2,CE2,CF2,CA3,CB3,CC3,
     $               CA12,CB12,CC12,CD12,CA13,CB13,CC13,CD13,
     $               CA23,CB23,CC23,CD23,CE23,CF23,CG23,CH23,
     $               A11,A32
      DOUBLE PRECISION A11(2),A32
      DOUBLE COMPLEX CA1(2,2),CB1(2,2),CA2(2,2),CB2(2,2),CC2(2,2),
     $               CD2(2,2),CE2(2,2),CF2(2,2),CA3,CB3,CC3,
     $               CA12(2,2),CB12(2,2),CC12(2,2),CD12(2,2),
     $               CA13(2),CB13(2),CC13(2),CD13(2),
     $               CA23(2),CB23(2),CC23(2),CD23(2),CE23(2),
     $               CF23(2),CG23(2),CH23(2)
      SAVE/ST3COUP/
C
      DOUBLE PRECISION ST3MAT,EB,EW,AP,BP,CP,X(2),
     $                 ENEU,PBPW,PBPNEU,PWPNEU,
     $                 M1SQ,PSB,M2SQ,PCH,M3SQ,PT,M1M2,M1M3,M2M3
      INTEGER I,J,N
C
CDefine the energies that we are integrating over. EB comes from
Ckinematics and EW from the phase space detla-fn. The integral is
Cbeing carried out over x(1) and x(2) which take values from 0 to 1
Cand have been substituted for E_b and E_W respectively
C
      EB=((MST(1)**2+MB**2-(MW+DSQRT(MNEU(1)**2))**2)/(2.D0*MST(1))-MB)
     .   *X(1)+MB
C
      AP=-4.D0*MB**2-4.D0*MST(1)**2+8.D0*MST(1)*EB
      BP=2.D0*(MST(1)**2+MB**2+MW**2-MNEU(1)**2-2.D0*MST(1)*EB)*(2.D0*
     .   MST(1)-2.D0*EB)
      CP=-(MST(1)**2+MB**2+MW**2-MNEU(1)**2-2.D0*MST(1)*EB)**2-4.D0*
     .   (EB**2-MB**2)*MW**2
C
CThe following relies on the fact that BP is positive and AP is negative.
CError message if not.
C
      EW=-DSQRT(BP**2-4.D0*AP*CP)/AP*X(2)+1.D0/(2.D0*AP)*(-BP+DSQRT
     .   (BP**2-4.D0*AP*CP))
      IF(BP.LT.0.D0.OR.AP.GT.0.D0)
     .            WRITE(*,*)'WARNING: KINEMATICS MAY BE INCORRECT',EW,EB
C
CEnergy, Momenta and Dot products which are used to simplify the algebra
C
      ENEU=MST(1)-EB-EW
      PSB=DSQRT(MST(1)**2+MW**2-2.D0*EW*MST(1))
      PCH=DSQRT(MST(1)**2+MB**2-2.D0*EB*MST(1))
      PT=DSQRT(MST(1)**2+MNEU(1)**2-2.D0*ENEU*MST(1))
      PBPW=(MST(1)**2+MNEU(1)**2-MB**2-MW**2)/2.D0-ENEU*MST(1)
      PBPNEU=(MST(1)**2+MW**2-MB**2-MNEU(1)**2)/2.D0-EW*MST(1)
      PWPNEU=(MST(1)**2+MB**2-MW**2-MNEU(1)**2)/2.D0-EB*MST(1)
C
CNow we are ready to define the different parts of the matrix element
C
CM1**2
C
      M1SQ=0.D0
      DO I=1,2
        DO J=1,2
C
          M1SQ=M1SQ+16.D0*A11(I)*A11(J)
C
CThe next 4 lines are the propagators
C
     .         *((PSB**2-MSB(I)**2)*(PSB**2-MSB(J)**2)
     .         +MSB(I)*MSB(J)*SBWIDTH(I)*SBWIDTH(J))
     .         /(((PSB**2-MSB(I)**2)**2+MSB(I)**2*SBWIDTH(I)**2)
     .         *((PSB**2-MSB(J)**2)**2+MSB(J)**2*SBWIDTH(J)**2))
C
     .         *((1.D0/MW**2*(PBPW**2+PWPNEU**2+2.D0*PBPW*PWPNEU)
     .         -MB**2-MNEU(1)**2-2.D0*PBPNEU)*(CA1(I,J)*PBPNEU-CB1(I,J)
     .         *MB
     .         *MNEU(1)))

        END DO
      END DO
C
CM2**2
C
      M2SQ=0.D0
      DO I=1,2
        DO J=1,2
C
          M2SQ=M2SQ+
C
CThe next 4 lines are the propagators
C
     .         ((PCH**2-MCH(I)**2)*(PCH**2-MCH(J)**2)
     .         +MCH(I)*MCH(J)*CHWIDTH(I)*CHWIDTH(J))
     .         /(((PCH**2-MCH(I)**2)**2+MCH(I)**2*CHWIDTH(I)**2)
     .         *((PCH**2-MCH(J)**2)**2+MCH(J)**2*CHWIDTH(J)**2))
C
     .         *(4.D0*CA2(I,J)*(PBPNEU*(MNEU(1)**2-MW**2)+4.D0*PWPNEU
     .         *(PBPW+PBPNEU)+2.D0*MNEU(1)**2*PBPW+2.D0/MW**2
     .         *PWPNEU*(2.D0*PBPNEU*PWPNEU-PBPW*MNEU(1)**2))
C
     .         +12.D0*CB2(I,J)*MNEU(1)*MB*MCH(I)*MCH(J)
C
     .         -12.D0*CC2(I,J)*MNEU(1)*(PBPNEU+PBPW)
C
     .         +CD2(I,J)*MCH(I)*MCH(J)*(8.D0/MW**2*PWPNEU*PBPW
     .         +4.D0*PBPNEU)
C
     .         +12.D0*CE2(I,J)*MNEU(1)*MB*(MW**2+MNEU(1)**2
     .         +2.D0*PWPNEU)
C
     .         -CF2(I,J)*MB*(12.D0*PWPNEU+8.D0/MW**2*PWPNEU**2
     .         +4.D0*MNEU(1)**2))
        END DO
      END DO
C
CM3**2
C
      M3SQ=1.D0/((PT**2-MT**2)**2+MT**2*TWIDTH**2)
C
     .     *(8.D0*CA3*(MB**2*(2.D0*PWPNEU+PBPNEU)
     .     +4.D0/MW**2*PBPNEU*PBPW**2-2.D0/MW**2*MB**2*PWPNEU*PBPW
     .     -MW**2*PBPNEU+4.D0*PBPW*(PBPNEU+PWPNEU))
C
     .     +8.D0*MT**2*CB3*(2.D0/MW**2*PWPNEU*PBPW+PBPNEU)
C
     .     -16.D0*MT*MNEU(1)*CC3*(3.D0*PBPW+2.D0/MW**2*PBPW**2+MB**2))

C
C2*M1M2
C
      M1M2=0.D0
      DO I=1,2
        DO J=1,2
C
          M1M2=M1M2+DBLE(8.D0*A11(J)
C
CThe next 4 lines are the propagators
C
     .         *((PCH**2-MCH(I)**2)*(PSB**2-MSB(J)**2)
     .         +MCH(I)*MSB(J)*CHWIDTH(I)*SBWIDTH(J))
     .         /(((PCH**2-MCH(I)**2)**2+MCH(I)**2*CHWIDTH(I)**2)
     .         *((PSB**2-MSB(J)**2)**2+MSB(J)**2*SBWIDTH(J)**2))
C
     .         *(CA12(I,J)*(MNEU(1)**2/MW**2*(2.D0*PBPW**2+2.D0*PWPNEU
     .         *PBPW)-2.D0*PBPNEU*PWPNEU-2.D0/MW**2*PBPNEU*(2.D0
     .         *PWPNEU**2+2.D0*PWPNEU*PBPW)-2.D0*MB**2*(MNEU(1)**2
     .         +PWPNEU)+2.D0*MNEU(1)**2*(PBPNEU+PBPW)+2.D0*PBPNEU
     .         *PBPW+4.D0*PBPNEU**2)
C
     .         -MNEU(1)*MCH(I)*CB12(I,J)*(2.D0*EB*MST(1)/MW**2
     .         *(EW*MST(1)
     .         -MW**2)+EW*MST(1)/MW**2*(MNEU(1)**2-MST(1)**2-MB**2
     .         -MW**2+2.D0*EW*MST(1)))
C
     .         +(MB*MNEU(1)*CC12(I,J)+MB*MCH(I)*CD12(I,J))
     .         *(1.D0/MW**2
     .         *(2.D0
     .         *PWPNEU**2+2.D0*PWPNEU*PBPW)-2.D0*MNEU(1)**2
     .         -2.D0*PBPNEU)))
C
        END DO
      END DO
C
C2*M1M3
C
      M1M3=0.D0
      DO I=1,2
C
        M1M3=M1M3+DBLE(8.D0*A32*A11(I)
C
CThe next 4 lines are the propagators
C
     .        *((PT**2-MT**2)*(PSB**2-MSB(I)**2)+MT*MSB(I)
     .        *TWIDTH*SBWIDTH(I))
     .        /(((PT**2-MT**2)**2+MT**2*TWIDTH**2)
     .        *((PSB**2-MSB(I)**2)**2+MSB(I)**2*SBWIDTH(I)**2))
C
     .        *(CA13(I)*(-MB**2/MW**2*(2.D0*PWPNEU**2+2.D0*PWPNEU
     .        *PBPW)+2.D0*PBPNEU*PBPW+2.D0/MW**2*PBPNEU*(2.D0
     .        *PBPW**2+2.D0*PWPNEU*PBPW)+2.D0*MNEU(1)**2*(MB**2
     .        +PBPW)-2.D0*MB**2*(PWPNEU+PBPNEU)-2.D0*PBPNEU*PWPNEU
     .        -4.D0*PBPNEU**2)
C
     .        +MB*MT*CB13(I)*(1.D0/MW**2*(2.D0*PWPNEU**2+2.D0*PWPNEU
     .        *PBPW)-2.D0*MNEU(1)**2-2.D0*PBPNEU)
C
     .        -(MT*MNEU(1)*CC13(I)+MB*MNEU(1)*CD13(I))
     .        *(1.D0/MW**2*(2.D0*PBPW**2+2.D0*PWPNEU*PBPW)-2.D0*MB**2
     .        -2.D0*PBPNEU)))
      END DO
C
C2*M2M3
C
      M2M3=0.D0
      DO I=1,2
C
        M2M3=M2M3+DBLE(8.D0*A32
C
CThe next 4 lines are the propagators
C
     .        *((PT**2-MT**2)*(PCH**2-MCH(I)**2)+MT*MCH(I)
     .        *TWIDTH*CHWIDTH(I))
     .        /(((PT**2-MT**2)**2+MT**2*TWIDTH**2)
     .        *((PCH**2-MCH(I)**2)**2+MCH(I)**2*CHWIDTH(I)**2))
C
     .        *(CA23(I)*(2.D0*MB**2/MW**2*PWPNEU**2-MNEU(1)**2*MB**2
     .        +2.D0*PBPNEU*PWPNEU-PBPNEU*MW**2-2.D0/MW**2*PBPW
     .        *(2.D0*PBPNEU*PWPNEU-PBPW*MNEU(1)**2)+PBPW*MNEU(1)**2
     .        +2.D0*PBPW*PBPNEU+PWPNEU*MB**2+4.D0*PBPNEU**2+4.D0
     .        *PBPW*PWPNEU)
C
     .        +MCH(I)*MT*CB23(I)*(2.D0/MW**2*PBPW*PWPNEU+PBPNEU)
C
     .        -3.D0*MNEU(1)*MT*CC23(I)*(PBPNEU+PBPW)
C
     .        -MNEU(1)*MCH(I)*CD23(I)*(2.D0/MW**2*PBPW**2+MB**2
     .        +3.D0*PBPW)
C
     .        -MB*MT*CE23(I)*(2.D0/MW**2*PWPNEU**2+MNEU(1)**2
     .        +3.D0*PWPNEU)
C
     .        -3.D0*MCH(I)*MB*CF23(I)*(PBPNEU+PWPNEU)
C
     .        +MB*MNEU(1)*CG23(I)*(3.D0*MW**2+3.D0*PWPNEU+3.D0*PBPW
     .        +PBPNEU+2.D0/MW**2*PWPNEU*PBPW)
C
     .        +3.D0*MCH(I)*MB*MNEU(1)*MT*CH23(I)))
      END DO
C
CThis is the final answer:
C
      ST3MAT=(M1SQ+M2SQ+M3SQ+M1M2+M1M3+M2M3)*((MST(1)**2+MB**2-
     .        (MW+DSQRT(MNEU(1)**2))**2)/(2.D0*MST(1))-MB)*(-DSQRT(BP**2
     .        -4.D0*AP*CP))/AP
C
      RETURN
      END

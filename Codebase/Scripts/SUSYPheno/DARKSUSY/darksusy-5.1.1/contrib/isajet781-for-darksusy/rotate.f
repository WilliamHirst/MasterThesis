!
      SUBROUTINE ROTATE(SMASS)
!
!Purpose: Rotate the parameters to the original current
!         basis
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3),TRIUP(3,3),TRIDP(3,3)
      DOUBLE COMPLEX MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3),TRIUDUM(3,3),TRIDDUM(3,3)
      DOUBLE COMPLEX MTSFUDUM(3,3),MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL      
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate from the quark mass basis
!Otherwise we rotate between different current bases.
!
      DO I=1,3
        DO J=1,3
          IF(SMASS.EQ.1)THEN
            VLUQ(I,J)=VLU(I,J)
            VLDQ(I,J)=VLD(I,J)
          ELSE
            IF(SVLQ.EQ.1)THEN
              VLUQ(I,J)=VLU(I,J)
              VLDQ(I,J)=VLU(I,J)
            ELSE
              VLUQ(I,J)=VLD(I,J)
              VLDQ(I,J)=VLD(I,J)
            END IF
          END IF
        END DO
      END DO
!
!Convert the G's to matrices
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=G(3+(I-1)*3+J)
          FDP(I,J)=G(12+(I-1)*3+J)
          MQP(I,J)=G(62+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MDP(I,J)=G(89+(I-1)*3+J)
          AUP(I,J)=G(33+(I-1)*3+J)
          ADP(I,J)=G(42+(I-1)*3+J)
!
          LUP(I,J)=G(111+(I-1)*3+J)
          LDP(I,J)=G(120+(I-1)*3+J)
!
          GTPQP(I,J)=G(138+(I-1)*3+J)
          GTPUP(I,J)=G(156+(I-1)*3+J)
          GTPDP(I,J)=G(165+(I-1)*3+J)
          GTQP(I,J)=G(185+(I-1)*3+J)
          GTSQP(I,J)=G(205+(I-1)*3+J)
          GTSUP(I,J)=G(214+(I-1)*3+J)
          GTSDP(I,J)=G(223+(I-1)*3+J)
          FTUQP(I,J)=G(232+(I-1)*3+J)
          FTDQP(I,J)=G(241+(I-1)*3+J)
          FTUUP(I,J)=G(259+(I-1)*3+J)
          FTDDP(I,J)=G(268+(I-1)*3+J)
!
          FMUP(I,J)=G(293+(I-1)*3+J)
          FMDP(I,J)=G(302+(I-1)*3+J)
          MQMP(I,J)=G(352+(I-1)*3+J)
          MUMP(I,J)=G(370+(I-1)*3+J)
          MDMP(I,J)=G(379+(I-1)*3+J)
          AUMP(I,J)=G(323+(I-1)*3+J)
          ADMP(I,J)=G(332+(I-1)*3+J)
!
          TRIUP(I,J)=G(399+(I-1)*3+J)
          TRIDP(I,J)=G(408+(I-1)*3+J)
          MTSFUP(I,J)=G(429+(I-1)*3+J)
          MTSFDP(I,J)=G(438+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(2,MQP,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(2,MUP,VRU,I,J)
          MDDUM(I,J)=CMATMUL(2,MDP,VRD,I,J)
          AUDUM(I,J)=CMATMUL(0,AUP,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(0,ADP,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(0,LUP,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(0,LDP,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(2,GTPQP,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(2,GTPUP,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(2,GTPDP,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(2,GTQP,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(2,GTSQP,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(2,GTSUP,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(2,GTSDP,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(0,FTUQP,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(0,FTDQP,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(0,FTUUP,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(0,FTDDP,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(0,FMUP,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(0,FMDP,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(2,MQMP,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(2,MUMP,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(2,MDMP,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(0,AUMP,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(0,ADMP,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(0,TRIUP,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(0,TRIDP,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(0,MTSFUP,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(0,MTSFDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=CMATMUL(1,VLUQT,FUDUM,I,J)
          G(12+(I-1)*3+J)=CMATMUL(1,VLDQT,FDDUM,I,J)
          G(62+(I-1)*3+J)=CMATMUL(0,VLUQ,MQDUM,I,J)
          G(80+(I-1)*3+J)=CMATMUL(0,VRU,MUDUM,I,J)
          G(89+(I-1)*3+J)=CMATMUL(0,VRD,MDDUM,I,J)
          G(33+(I-1)*3+J)=CMATMUL(1,VLUQT,AUDUM,I,J)
          G(42+(I-1)*3+J)=CMATMUL(1,VLDQT,ADDUM,I,J)
!
          G(111+(I-1)*3+J)=CMATMUL(1,VLUQT,LUDUM,I,J)
          G(120+(I-1)*3+J)=CMATMUL(1,VLDQT,LDDUM,I,J)
          G(138+(I-1)*3+J)=CMATMUL(0,VLUQ,GTPQDUM,I,J)
          G(156+(I-1)*3+J)=CMATMUL(0,VRU,GTPUDUM,I,J)
          G(165+(I-1)*3+J)=CMATMUL(0,VRD,GTPDDUM,I,J)
          G(185+(I-1)*3+J)=CMATMUL(0,VLUQ,GTQDUM,I,J)
          G(205+(I-1)*3+J)=CMATMUL(0,VLUQ,GTSQDUM,I,J)
          G(214+(I-1)*3+J)=CMATMUL(0,VRU,GTSUDUM,I,J)
          G(223+(I-1)*3+J)=CMATMUL(0,VRD,GTSDDUM,I,J)
          G(232+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUQDUM,I,J)
          G(241+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDQDUM,I,J)
          G(259+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUUDUM,I,J)
          G(268+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDDDUM,I,J)
!
          G(293+(I-1)*3+J)=CMATMUL(1,VLUQT,FMUDUM,I,J)
          G(302+(I-1)*3+J)=CMATMUL(1,VLDQT,FMDDUM,I,J)
          G(352+(I-1)*3+J)=CMATMUL(0,VLUQ,MQMDUM,I,J)
          G(370+(I-1)*3+J)=CMATMUL(0,VRU,MUMDUM,I,J)
          G(379+(I-1)*3+J)=CMATMUL(0,VRD,MDMDUM,I,J)
          G(323+(I-1)*3+J)=CMATMUL(1,VLUQT,AUMDUM,I,J)
          G(332+(I-1)*3+J)=CMATMUL(1,VLDQT,ADMDUM,I,J)
!
          G(399+(I-1)*3+J)=CMATMUL(1,VLUQT,TRIUDUM,I,J)
          G(408+(I-1)*3+J)=CMATMUL(1,VLDQT,TRIDDUM,I,J)
          G(429+(I-1)*3+J)=CMATMUL(1,VLUQT,MTSFUDUM,I,J)
          G(438+(I-1)*3+J)=CMATMUL(1,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END

!
      SUBROUTINE ROTBACK(SMASS)
!
!Purpose: Rotate all the matrices from the original current basis.
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
      DOUBLE COMPLEX FU(3,3),FD(3,3),MQ(3,3),MU(3,3),MD(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),LU(3,3),LD(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPU(3,3),GTPD(3,3),GTQ(3,3)
      DOUBLE COMPLEX GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTUU(3,3),FTDD(3,3)
      DOUBLE COMPLEX FMU(3,3),FMD(3,3),MQM(3,3),MUM(3,3),MDM(3,3)
      DOUBLE COMPLEX AUM(3,3),ADM(3,3)
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),MTSFU(3,3),MTSFD(3,3)
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3)
      DOUBLE COMPLEX TRIUP(3,3),TRIDP(3,3),MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3)
      DOUBLE COMPLEX TRIUDUM(3,3),TRIDDUM(3,3),MTSFUDUM(3,3)
      DOUBLE COMPLEX MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate to the quark mass basis
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
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO   
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          AU(I,J)=G(33+(I-1)*3+J)
          AD(I,J)=G(42+(I-1)*3+J)
          MQ(I,J)=G(62+(I-1)*3+J)
          MU(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
!
          LU(I,J)=G(111+(I-1)*3+J)
          LD(I,J)=G(120+(I-1)*3+J)
!
          GTPQ(I,J)=G(138+(I-1)*3+J)
          GTPU(I,J)=G(156+(I-1)*3+J)
          GTPD(I,J)=G(165+(I-1)*3+J)
          GTQ(I,J)=G(185+(I-1)*3+J)
          GTSQ(I,J)=G(205+(I-1)*3+J)
          GTSU(I,J)=G(214+(I-1)*3+J)
          GTSD(I,J)=G(223+(I-1)*3+J)
          FTUQ(I,J)=G(232+(I-1)*3+J)
          FTDQ(I,J)=G(241+(I-1)*3+J)
          FTUU(I,J)=G(259+(I-1)*3+J)
          FTDD(I,J)=G(268+(I-1)*3+J)
!
          FMU(I,J)=G(293+(I-1)*3+J)
          FMD(I,J)=G(302+(I-1)*3+J)
          AUM(I,J)=G(323+(I-1)*3+J)
          ADM(I,J)=G(332+(I-1)*3+J)
          MQM(I,J)=G(352+(I-1)*3+J)
          MUM(I,J)=G(370+(I-1)*3+J)
          MDM(I,J)=G(379+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          TRID(I,J)=G(408+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(2,FU,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(2,FD,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(0,MQ,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(0,MU,VRU,I,J)
          MDDUM(I,J)=CMATMUL(0,MD,VRD,I,J)
          AUDUM(I,J)=CMATMUL(2,AU,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(2,AD,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(2,LU,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(2,LD,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(0,GTPQ,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(0,GTPU,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(0,GTPD,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(0,GTQ,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(0,GTSQ,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(0,GTSU,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(0,GTSD,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(2,FTUQ,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(2,FTDQ,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(2,FTUU,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(2,FTDD,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(2,FMU,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(2,FMD,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(0,MQM,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(0,MUM,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(0,MDM,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(2,AUM,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(2,ADM,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(2,TRIU,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(2,TRID,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(2,MTSFU,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(2,MTSFD,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUQT,FUDUM,I,J)
          FDP(I,J)=CMATMUL(0,VLDQT,FDDUM,I,J)
          MQP(I,J)=CMATMUL(1,VLUQ,MQDUM,I,J)
          MUP(I,J)=CMATMUL(1,VRU,MUDUM,I,J)
          MDP(I,J)=CMATMUL(1,VRD,MDDUM,I,J)
          AUP(I,J)=CMATMUL(0,VLUQT,AUDUM,I,J)
          ADP(I,J)=CMATMUL(0,VLDQT,ADDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUQT,LUDUM,I,J)
          LDP(I,J)=CMATMUL(0,VLDQT,LDDUM,I,J)
          GTPQP(I,J)=CMATMUL(1,VLUQ,GTPQDUM,I,J)
          GTPUP(I,J)=CMATMUL(1,VRU,GTPUDUM,I,J)
          GTPDP(I,J)=CMATMUL(1,VRD,GTPDDUM,I,J)
          GTQP(I,J)=CMATMUL(1,VLUQ,GTQDUM,I,J)
          GTSQP(I,J)=CMATMUL(1,VLUQ,GTSQDUM,I,J)
          GTSUP(I,J)=CMATMUL(1,VRU,GTSUDUM,I,J)
          GTSDP(I,J)=CMATMUL(1,VRD,GTSDDUM,I,J)
          FTUQP(I,J)=CMATMUL(0,VLUQT,FTUQDUM,I,J)
          FTDQP(I,J)=CMATMUL(0,VLDQT,FTDQDUM,I,J)
          FTUUP(I,J)=CMATMUL(0,VLUQT,FTUUDUM,I,J)
          FTDDP(I,J)=CMATMUL(0,VLDQT,FTDDDUM,I,J)
          FMUP(I,J)=CMATMUL(0,VLUQT,FMUDUM,I,J)
          FMDP(I,J)=CMATMUL(0,VLDQT,FMDDUM,I,J)
          MQMP(I,J)=CMATMUL(1,VLUQ,MQMDUM,I,J)
          MUMP(I,J)=CMATMUL(1,VRU,MUMDUM,I,J)
          MDMP(I,J)=CMATMUL(1,VRD,MDMDUM,I,J)
          AUMP(I,J)=CMATMUL(0,VLUQT,AUMDUM,I,J)
          ADMP(I,J)=CMATMUL(0,VLDQT,ADMDUM,I,J)
          TRIUP(I,J)=CMATMUL(0,VLUQT,TRIUDUM,I,J)
          TRIDP(I,J)=CMATMUL(0,VLDQT,TRIDDUM,I,J)
          MTSFUP(I,J)=CMATMUL(0,VLUQT,MTSFUDUM,I,J)
          MTSFDP(I,J)=CMATMUL(0,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=FUP(I,J)
          G(12+(I-1)*3+J)=FDP(I,J)
          G(62+(I-1)*3+J)=MQP(I,J)
          G(80+(I-1)*3+J)=MUP(I,J)
          G(89+(I-1)*3+J)=MDP(I,J)
          G(33+(I-1)*3+J)=AUP(I,J)
          G(42+(I-1)*3+J)=ADP(I,J)
!
          G(111+(I-1)*3+J)=LUP(I,J)
          G(120+(I-1)*3+J)=LDP(I,J)
!
          G(138+(I-1)*3+J)=GTPQP(I,J)
          G(156+(I-1)*3+J)=GTPUP(I,J)
          G(165+(I-1)*3+J)=GTPDP(I,J)
          G(185+(I-1)*3+J)=GTQP(I,J)
          G(205+(I-1)*3+J)=GTSQP(I,J)
          G(214+(I-1)*3+J)=GTSUP(I,J)
          G(223+(I-1)*3+J)=GTSDP(I,J)
          G(232+(I-1)*3+J)=FTUQP(I,J)
          G(241+(I-1)*3+J)=FTDQP(I,J)
          G(259+(I-1)*3+J)=FTUUP(I,J)
          G(268+(I-1)*3+J)=FTDDP(I,J)
!
          G(293+(I-1)*3+J)=FMUP(I,J)
          G(302+(I-1)*3+J)=FMDP(I,J)
          G(352+(I-1)*3+J)=MQMP(I,J)
          G(370+(I-1)*3+J)=MUMP(I,J)
          G(379+(I-1)*3+J)=MDMP(I,J)
          G(323+(I-1)*3+J)=AUMP(I,J)
          G(332+(I-1)*3+J)=ADMP(I,J)
!
          G(399+(I-1)*3+J)=TRIUP(I,J)
          G(408+(I-1)*3+J)=TRIDP(I,J)
          G(429+(I-1)*3+J)=MTSFUP(I,J)
          G(438+(I-1)*3+J)=MTSFDP(I,J)
        END DO
      END DO
!
      RETURN
      END

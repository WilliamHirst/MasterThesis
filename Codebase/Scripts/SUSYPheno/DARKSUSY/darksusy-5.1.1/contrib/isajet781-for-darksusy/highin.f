!
      SUBROUTINE HIGHIN(TIME)
!
!Purpose: To put the high scale conditions into G(601) ready
!         for running back down
!
!The switch 'TIME' is used to tell if this is the first call to
!HIGHIN. If so, the MSSM mu parameter is also set here.
!In future calls it will be set along with b at m_H.
!
!NOTE ABOUT INPUT OF SOFT MATRICES:
!
!All soft matrices are input in the current basis.
!To ensure that it is possible to have good control over additional
!flavour violating inputs, it is possible to set the coefficients
!for three independent matrices which have been carefully chosen 
!so that they do not introduce new flavour structure. In addition
!to the coefficients of these three matrices, there is also the
!possibility for introducing an entirely general matrix, which should
!be chosen with care.
!
!a_{ude}(GUT)=f_{ude}[A_{{ude}}\1+W_{ude}f^\dagger_{ude}f_{ude}
!         +X_{ude}f^\dagger_{ude}f_{ude}f^\dagger_{ude}f_{ude}]+Z_{ude}
!
!m^2_{UDE}(GUT)=m^2_{{UDE}0}[\1+R_{UDE}f^T_{ude}f^*_{ude}
!                 +S_{UDE}f^T_{ude}f^*_{ude}f^T_{ude}f^*_{ude}]+T_{UDE}
!
!m^2_{QL}(GUT)=m^2_{{QL}0}\1+T_{QL}
!
!where \1 is the unit matrix, denoted in the code by CID(I,J).
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/HISUSY/CM1,CM2,CM3,CM1P,CM2P,CM3P
     $             ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $             ,CTQ,CTU,CTD,CTL,CTE
     $             ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $             ,CZU,CZD,CZE,CMHU,CMHD,CU,CD,CE
      DOUBLE COMPLEX CM1,CM2,CM3,CM1P,CM2P,CM3P
     $              ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $              ,CTQ(3,3),CTU(3,3),CTD(3,3),CTL(3,3),CTE(3,3)
     $              ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $              ,CZU(3,3),CZD(3,3),CZE(3,3),CMHU,CMHD
      INTEGER CU,CD,CE
      SAVE/HISUSY/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
      COMMON /RGEFNM/ FNRGE
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3),FUS(3,3),FDS(3,3),FES(3,3)
     $              ,AU(3,3),AD(3,3),AE(3,3)
     $              ,MQSQ(3,3),MLSQ(3,3),MUPSQ(3,3),MDSQ(3,3),MESQ(3,3)
     $              ,FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
     $              ,FUDFUFUDFU(3,3),FDDFDFDDFD(3,3),FEDFEFEDFE(3,3)
     $              ,FUTFUS(3,3),FDTFDS(3,3),FETFES(3,3)
     $              ,FUTFUSFUTFUS(3,3),FDTFDSFDTFDS(3,3)
     $              ,FETFESFETFES(3,3)
      DOUBLE COMPLEX DUMAU(3,3),DUMAD(3,3),DUMAE(3,3),FUSFUT(3,3)
     $              ,FDSFDT(3,3)
      DOUBLE COMPLEX CID(3,3),CMATMUL,GTMP(601)
!
      DOUBLE PRECISION MUSQ,TU,TD
      INTEGER I,J,K,TIME
      CHARACTER*128 FILENAME,STRADD,FNRGE
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Set the initial values for the old number of matter sfermions
!and their rotations.
!
      OLDNSQ=3
      OLDNSU=3
      OLDNSD=3
      OLDNSL=3
      OLDNSE=3
      DO I=1,3
        THSQ(I)=1
        THSU(I)=1
        THSD(I)=1
        THSL(I)=1
        THSE(I)=1
      END DO
      DO I=1,3
        DO J=1,3
          RQTOT(I,J)=CID(I,J)
          RUPTOT(I,J)=CID(I,J)
          RDTOT(I,J)=CID(I,J)
          RLTOT(I,J)=CID(I,J)
          RETOT(I,J)=CID(I,J)
          DO K=1,2
            RQSAV(K,I,J)=CID(I,J)
            RUPSAV(K,I,J)=CID(I,J)
            RDSAV(K,I,J)=CID(I,J)
            RLSAV(K,I,J)=CID(I,J)
            RESAV(K,I,J)=CID(I,J)            
          END DO
        END DO
      END DO
!
!Convert the G's to Matrices
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
        END DO
      END DO
!
!If we are using mSUGRA boundary conditions
!
      IF(SUG.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            AU(I,J)=DCMPLX(A0)*FU(I,J)
            AD(I,J)=DCMPLX(A0)*FD(I,J)
            AE(I,J)=DCMPLX(A0)*FE(I,J)
            MQSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MLSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MUPSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MDSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MESQ(I,J)=DCMPLX(M0**2)*CID(I,J)
          END DO
        END DO
!
        G(31)=DCMPLX(M12)
        G(32)=DCMPLX(M12)
        G(33)=DCMPLX(M12)
        G(599)=(0.D0,0.D0)
        G(600)=(0.D0,0.D0)
        G(601)=(0.D0,0.D0)
        IF(TIME.EQ.1)THEN
          G(61)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
          G(62)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
        ELSE
!
!(G(61)+G(62)) WAS FIXED BY EWSB AT M_SUSY
!ABS(G(108)) WAS FIXED AT HIGHEST SUSY THRESHOLD
!
          G(108)=ABS(G(108))*EXP((0.D0,1.D0)*PHASEMU)
          G(61)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
          G(62)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
        END IF
        G(351)=DCMPLX(M0**2)
        G(352)=DCMPLX(M0**2)
!
!Or use the non-mSUGRA matrices entered by the user
!
      ELSE
!
        DO I=1,3
          DO J=1,3
            FUS(I,J)=CONJG(FU(I,J))
            FDS(I,J)=CONJG(FD(I,J))
            FES(I,J)=CONJG(FE(I,J))
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            FUDFU(I,J)=CMATMUL(1,FU,FU,I,J)
            FDDFD(I,J)=CMATMUL(1,FD,FD,I,J)
            FEDFE(I,J)=CMATMUL(1,FE,FE,I,J)
            FUTFUS(I,J)=CMATMUL(1,FUS,FUS,I,J)
            FDTFDS(I,J)=CMATMUL(1,FDS,FDS,I,J)
            FETFES(I,J)=CMATMUL(1,FES,FES,I,J)
            FUSFUT(I,J)=CMATMUL(2,FUS,FUS,I,J)
            FDSFDT(I,J)=CMATMUL(2,FDS,FDS,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            FUTFUSFUTFUS(I,J)=CMATMUL(0,FUTFUS,FUTFUS,I,J)
            FDTFDSFDTFDS(I,J)=CMATMUL(0,FDTFDS,FDTFDS,I,J)
            FETFESFETFES(I,J)=CMATMUL(0,FETFES,FETFES,I,J)
            FUDFUFUDFU(I,J)=CMATMUL(0,FUDFU,FUDFU,I,J)
            FDDFDFDDFD(I,J)=CMATMUL(0,FDDFD,FDDFD,I,J)
            FEDFEFEDFE(I,J)=CMATMUL(0,FEDFE,FEDFE,I,J)
          END DO
        END DO
!
!Here we use the relations in the comment at the start of this
!subroutine
!
        DO I=1,3
          DO J=1,3
            DUMAU(I,J)=CAU0*CID(I,J)+CWU*FUDFU(I,J)+CXU*FUDFUFUDFU(I,J)
            DUMAD(I,J)=CAD0*CID(I,J)+CWD*FDDFD(I,J)+CXD*FDDFDFDDFD(I,J)
            DUMAE(I,J)=CAE0*CID(I,J)+CWE*FEDFE(I,J)+CXE*FEDFEFEDFE(I,J)
          END DO
        END DO
!
!The next two lines can be hardwired to allow for entry of a T_Q
!which is proportional to the Yukawas, a la Minimal Flavour Violation
!See arXiv:0810.5765.
!
        TU=0.D0
        TD=0.D0
        IF(TU.NE.0.D0.OR.TD.NE.0.D0)THEN
          WRITE(*,*)'NON-ZERO MFV CONDITIONS HARDWIRED AT GUT'
          WRITE(*,*)'TU=',TU
          WRITE(*,*)'TD=',TD
        END IF
        DO I=1,3
          DO J=1,3
            AU(I,J)=CMATMUL(0,FU,DUMAU,I,J)+CZU(I,J)
            AD(I,J)=CMATMUL(0,FD,DUMAD,I,J)+CZD(I,J)
            AE(I,J)=CMATMUL(0,FE,DUMAE,I,J)+CZE(I,J)
!
            MQSQ(I,J)=CMQ0**2*CID(I,J)+CTQ(I,J)
            MQSQ(I,J)=CMQ0**2*(CID(I,J)+TU*FUSFUT(I,J)
     $                        +TD*FDSFDT(I,J))+CTQ(I,J)
            MUPSQ(I,J)=CMU0**2*(DBLE(CU)*CID(I,J)+CRU*FUTFUS(I,J)
     $                         +CSU*FUTFUSFUTFUS(I,J))+CTU(I,J)
            MDSQ(I,J)=CMD0**2*(DBLE(CD)*CID(I,J)+CRD*FDTFDS(I,J)
     $                         +CSD*FDTFDSFDTFDS(I,J))+CTD(I,J)
            MLSQ(I,J)=CML0**2*CID(I,J)+CTL(I,J)
            MESQ(I,J)=CME0**2*(DBLE(CE)*CID(I,J)+CRE*FETFES(I,J)
     $                         +CSE*FETFESFETFES(I,J))+CTE(I,J)
          END DO
        END DO
!
        G(31)=CM1
        G(32)=CM2
        G(33)=CM3
        G(599)=CM1P
        G(600)=CM2P
        G(601)=CM3P
        IF(TIME.EQ.1)THEN
          G(61)=CMHU**2+G(108)*CONJG(G(108))
          G(62)=CMHD**2+G(108)*CONJG(G(108))
        ELSE
!
!(G(61)+G(62)) WAS FIXED BY EWSB AT M_SUSY
!ABS(G(108)) WAS FIXED AT HIGHEST SUSY THRESHOLD
!
          G(108)=ABS(G(108))*EXP((0.D0,1.D0)*PHASEMU)
          G(61)=CMHU**2+G(108)*CONJG(G(108))
          G(62)=CMHD**2+G(108)*CONJG(G(108))
        END IF
        G(351)=CMHU**2
        G(352)=CMHD**2
      END IF
!
!Now put back into the g's
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=FU(I,J)
          G(12+(I-1)*3+J)=FD(I,J)
          G(21+(I-1)*3+J)=FE(I,J)
          G(33+(I-1)*3+J)=AU(I,J)
          G(42+(I-1)*3+J)=AD(I,J)
          G(51+(I-1)*3+J)=AE(I,J)
          G(62+(I-1)*3+J)=MQSQ(I,J)
          G(71+(I-1)*3+J)=MLSQ(I,J)
          G(80+(I-1)*3+J)=MUPSQ(I,J)
          G(89+(I-1)*3+J)=MDSQ(I,J)
          G(98+(I-1)*3+J)=MESQ(I,J)
        END DO
      END DO
!
!Next, I need to set the starting values for the tilde couplings and
!MSSM terms
!
      DO I=1,3
        DO J=1,3
          G(138+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(147+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(156+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(165+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(174+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(185+(I-1)*3+J)=G(2)*CID(I,J)
          G(194+(I-1)*3+J)=G(2)*CID(I,J)
          G(205+(I-1)*3+J)=G(3)*CID(I,J)
          G(214+(I-1)*3+J)=G(3)*CID(I,J)
          G(223+(I-1)*3+J)=G(3)*CID(I,J)
          G(232+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(241+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(250+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(259+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(268+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(277+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(293+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(302+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(311+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(323+(I-1)*3+J)=G(33+(I-1)*3+J)
          G(332+(I-1)*3+J)=G(42+(I-1)*3+J)
          G(341+(I-1)*3+J)=G(51+(I-1)*3+J)
          G(352+(I-1)*3+J)=G(62+(I-1)*3+J)
          G(361+(I-1)*3+J)=G(71+(I-1)*3+J)
          G(370+(I-1)*3+J)=G(80+(I-1)*3+J)
          G(379+(I-1)*3+J)=G(89+(I-1)*3+J)
          G(388+(I-1)*3+J)=G(98+(I-1)*3+J)
          G(429+(I-1)*3+J)=CONJG(G(108))*G(3+(I-1)*3+J)
          G(438+(I-1)*3+J)=CONJG(G(108))*G(12+(I-1)*3+J)
          G(447+(I-1)*3+J)=CONJG(G(108))*G(21+(I-1)*3+J)
        END DO
      END DO
      G(184)=DSQRT(3.D0/5.D0)*G(1)
      G(185)=DSQRT(3.D0/5.D0)*G(1)
      G(204)=G(2)
      G(205)=G(2)
      G(291)=G(1)
      G(292)=G(2)
      G(293)=G(3)
      G(321)=G(31)-(0.D0,1.D0)*G(599)
      G(322)=G(32)-(0.D0,1.D0)*G(600)
      G(323)=G(33)-(0.D0,1.D0)*G(601)
      IF(TIME.EQ.1)THEN
        G(398)=G(108)
      END IF
!
!Rotate the the basis where the Yukawas are diagonal at m_t and
!write the results to a file.
!
      DO I=1,601 !Save the original couplings
        GTMP(I)=G(I)
      END DO
      CALL ROTBACK(0)
      FILENAME=STRADD(FNRGE,'.gtout')
      CALL GOUT601(G,MHIGH,1,FILENAME)
!
!Replace the rotated couplings by the old saved ones
!This removes any inaccuracies introduced by the rotation.
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
      RETURN
      END

!
      SUBROUTINE RGEFLAV(RM0,RM12,RA0,RTANB,RSIGNMU,RMT)
!
!Purpose: To run the RGEs and derive values for all the parameters at
!         M_t/m_H
!
!Hard wired constants:
!        KM(3,3) - the mixing matrix
!        MWEAK(6) - the low scale masses
!
!Inputs - see inrge.dat:
!        High scale parameters
!        MHIGH - Value of high scale if not GUT scale
!        SUG - Whether to use mSUGRA
!        UNI - Whether to run up to unified gauge couplings 
!        ACC - Accuracy of the RGEs (ie. one or two loops)
!        COMP - Whether to compute the complex RGEs
!        PHASEMU - phase for \mu at the EWSB scale
!        SVLQ - Basis in which to output results
!        VLU,VRU,VLD,VRD - Quark rotation matrices
!        g(1-3),f_b,f_tau,lambda,vev (at M_z) and f_t(M_t)
!
!Some of the inputs are passed since RGEFLAV is used as a subroutine
!
!Requires:
!        EISPACK - cg.f and dependencies
!        RGEs - crkstp.f, crgeutils.f, csmrgedr.f, crge215.o, crge601.o
!               drkstp.f, drgeutils.f, dsmrgedr.f, drge215.o, drge601.o
!        INPUTS - rgeread.f,inrge.dat
!        OTHER - sugefffl.f,sqsix.f
!
!Output - Writes to a file (weakout.dat) the full complement of
!         couplings and masses at m_H. These are written
!         in the basis where one of either the up- or down-type
!         Yukawa matrices are diagonal, depending on user input.
!
!         If SQSIX is called, there is also output to the files sqm2u.dat
!         and sqm2d.dat which contian the (6x6) up-type and down-type
!         squark mass matrices and (when in the appropriate basis)
!         corresponding eigenvectors and eigenvalues.
!
!         In addition, there are commented lines in:
!            UPMZMHIGH
!            DOWNMHIGHMZ
!            UPMMHIGH2
!            DECRUN (located in sqsix.f)
!
!         that can be reinstated to print to separate files the running
!         of various parameters. The lines include OPEN and CLOSE
!         commands as well as calls to OUTCOUP, UOUTCOUP and DECOUTCOUP
!         NB: These commands require the presence of a subdirectory
!             called 'out', which must be created before running the
!             executable.
!
!HOW TO SEPARATE OUT THE GUT SCALE FROM THE FLAVOUR SCALE:
!
!   If at some later time it is required that the GUT scale can
!   be separated from the flavour scale, it may be easiest
!   to implement by slightly altering DOWNMSCOND. Simply
!   run all the way to the GUT scale, but at each run down,
!   check if the scale has passed the required flavour scale.
!   When the running passes this scale, call the new subroutine
!   that would contain the flavour BCs.
!
!
!NOTE: For the time being, RGEFLAV must be called with
!      equivalent mSUGRA inputs. This is not a necessity,
!      but at present mSUGRA is the default GUT scale input.
!
      IMPLICIT NONE
!
!     SUGRA parameters passed by ISAJET
!
      COMMON /RGEFNM/ FNRGE
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
!     All 601 RGE terms. Definitions are in *rge601.f
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
!     Used to tell RGEs that we want to run to
!     two loops
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
!     NU is the number of active up-like quarks and SMDR2LP tells
!     csmrgedr.f to run to two loops
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
!     Constants read in by RGEREAD from a file
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
!Warning flag from the accuracy subroutine: ORTH
!
      COMMON/ORTHWARNING/ORTHFLAG
      INTEGER ORTHFLAG
      SAVE/ORTHWARNING/
!
      INTEGER I,SWSQSIX,NSTEP,FINALIT,THCNG
      DOUBLE PRECISION QEND,PI
      REAL RM0,RM12,RA0,RTANB,RSIGNMU,RMT
      CHARACTER*128 FILENAME,STRADD,FNRGE
!
      PI=4.D0*DATAN(1.D0)
!
!Convert real inputs to double precision.
!I change some signs here since RGEFLAV uses book notation
!and ISAJET inputs are in MV notation.
!
      M0=DBLE(RM0)
      M12=DBLE(RM12)
      A0=-DBLE(RA0)
      TANB=DBLE(RTANB)
      SIGNMU=-DBLE(RSIGNMU)
      PHASEMU=PI*(1.D0-SIGNMU)/2.D0
      MT=RMT
!
      SWSQSIX=1
!
!May need to be increased to allow for more accurate convergence
!of iterations.
!
      FINALIT=16
      THCNG=0
!
!Next the rest of the inputs from INRGE.DAT
!AS PART OF RGEREAD, READ IN THE SWITCHES ACC AND COMP
!
      CALL RGEREAD
!
      SW2LP=ACC   !These two lines set the accuracy of the MSSM
      SMDR2LP=ACC !and SM running
!
      DO I=1,601
        G(I)=(0.D0,0.D0)
      END DO
!
!The light weak-scale quark masses are derived at Mz
!
      CALL MASS
!
!We are ready to run up.
!
      CALL UPMZMHIGH
!
!Now insert the high scale conditions read in earlier
!
      CALL HIGHIN(1)
!
!Run back down to the weak scale and iterate
!
      DO I=1,FINALIT
!
        ORTHFLAG=0
!
!May need to be increased to allow for more accurate running.
!
        NSTEP=INT(100.D0*1.6D0**I)
        WRITE(*,*)'ITERATION NUMBER ',I
!
!Sets when the iteration of the threshold locations ends. Could be
!set based upon a test of convergence. We continue iterating
!after fixing the thresholds since if the thresholds are still
!changing the Yukawas will not converge on a solution where
!they are diagonal at M_t.
!
        IF(I.GT.5)THEN
          NSTEP=INT(100.D0*1.6D0**5)
          IF(I.GT.10)THCNG=1
        END IF
        CALL DOWNMHIGHMZ(QEND,NSTEP,THCNG)
        CALL UPMZMHIGH2(NSTEP)
        CALL HIGHIN(0)
!
!If ORTHFLAG is too large, there is a high probability that there
!have been accuracy issues in the rotation between the squark
!mass basis and the output basis. ORTHFLAG obtains contributions
!to its value whenever ORTH is unlikely to have been able to shift
!the error in the diagonalisation into an entry that is much larger
!than the size of the error.
!Print a warning
!
        IF(ORTHFLAG.GT.INT(NSTEP/10))WRITE(*,*)'WARNING FROM ORTHFLAG'
!
      END DO
!
!One final run back down
!
      IF(SWSQSIX.EQ.0)THEN
        CALL DOWNMHIGHMZ(QEND,NSTEP,1)
        CALL ROTBACK(0)
!
!Write the results to a file
!
        FILENAME=STRADD(FNRGE,'.wkout')
        CALL GOUT601(G,QEND,3,FILENAME)
!
      ELSE
!
!If required, the next lines are the call to the decay calculation prog.
!
        IF(SVLQ.EQ.0)WRITE(*,*)'DOWN BASIS CHOSEN - NO STOP CALCULATION'
        CALL DOWNMHIGHMH(QEND,NSTEP)
        CALL ROTBACK(0)
!
!Write the results to a file
!
        FILENAME=STRADD(FNRGE,'.wkout')
        CALL GOUT601(G,QEND,3,FILENAME)
!
        CALL SQSIX(G,QEND)
      END IF
!
      RETURN
      END

!
      SUBROUTINE VGEN(T12,T13,T23,D13,MAT)
!
!Purpose: To generate an arbitrary unitary matrix
!         Is passed the values of alpha (T12), beta (T13),
!         gamma (T23) and delta (D13)
!
      IMPLICIT NONE
!
      INTEGER I,J
      DOUBLE PRECISION T12,C12,S12,T13,C13,S13,T23,C23,S23,D13
      DOUBLE COMPLEX MAT(3,3),MAT1(3,3),MAT2(3,3),MAT3(3,3)
     $              ,TMPMAT(3,3),CMATMUL
!
!Find the values of sin and cos of the angles
!
      C12=DCOS(T12)
      S12=DSIN(T12)
      C13=DCOS(T13)
      S13=DSIN(T13)
      C23=DCOS(T23)
      S23=DSIN(T23)
!
!These are the individual matrices
!
      MAT1(1,1)=(1.D0,0.D0)
      MAT1(1,2)=(0.D0,0.D0)
      MAT1(1,3)=(0.D0,0.D0)
      MAT1(2,1)=(0.D0,0.D0)
      MAT1(2,2)=DCMPLX(C23)
      MAT1(2,3)=DCMPLX(S23)
      MAT1(3,1)=(0.D0,0.D0)
      MAT1(3,2)=-DCMPLX(S23)
      MAT1(3,3)=DCMPLX(C23)
      MAT2(1,1)=DCMPLX(C13)
      MAT2(1,2)=(0.D0,0.D0)
      MAT2(1,3)=S13*DCMPLX(DCOS(D13),DSIN(D13))
      MAT2(2,1)=(0.D0,0.D0)
      MAT2(2,2)=(1.D0,0.D0)
      MAT2(2,3)=(0.D0,0.D0)
      MAT2(3,1)=-S13*DCMPLX(DCOS(D13),-DSIN(D13))
      MAT2(3,2)=(0.D0,0.D0)
      MAT2(3,3)=DCMPLX(C13)
      MAT3(1,1)=DCMPLX(C12)
      MAT3(1,2)=DCMPLX(S12)
      MAT3(1,3)=(0.D0,0.D0)
      MAT3(2,1)=-DCMPLX(S12)
      MAT3(2,2)=DCMPLX(C12)
      MAT3(2,3)=(0.D0,0.D0)
      MAT3(3,1)=(0.D0,0.D0)
      MAT3(3,2)=(0.D0,0.D0)
      MAT3(3,3)=(1.D0,0.D0)
!
      DO I=1,3
        DO J=1,3
          TMPMAT(I,J)=CMATMUL(0,MAT1,MAT2,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MAT(I,J)=CMATMUL(0,TMPMAT,MAT3,I,J)
        END DO
      END DO
!
      RETURN
      END

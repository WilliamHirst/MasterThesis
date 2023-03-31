!
      SUBROUTINE MASS
!
!Purpose: Input the low scale masses of the quarks and leptons
!         in the order u,c,d,s,e,mu and thus derive the yukawas.
!         From paper by Fusaoka and Koide.
!
      IMPLICIT NONE
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      DOUBLE PRECISION MATQ(2,3,3),MATL(3,3)
      INTEGER I,J
!
!Use same (Pierce) prescription as SUGRA.F for vev at MZ
!Eq. 19 NPB491 p.3 (1997)
!Remove log since it is presumably the bulk of the SUSY 
!contribution if all sparticles are decoupled at M_Z
!
!      VSMMZ=DCMPLX((248.6D0+0.9D0*LOG(RGEMS/MZ))/DSQRT(2.D0))
      VSMMZ=DCMPLX(248.6D0/DSQRT(2.D0))
!
      MWEAK(1)=.00233D0
      MWEAK(2)=.677D0
      MWEAK(3)=.00469D0
      MWEAK(4)=.0934D0
      MWEAK(5)=0.4873D-3
      MWEAK(6)=102.87D-3
!
!The values of f_b, f_t and f_tau are derived from RGEs
!
      RETURN
      END

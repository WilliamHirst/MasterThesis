C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP4X(SW2,ZMASS,FMASS , GCHF)
C
C This subroutine sets up the coupling constant for the fermion-fermion-
C Higgs vertex in the STANDARD MODEL.  The coupling is COMPLEX and the  
C array of the coupling specifies the chirality of the flowing-IN       
C fermion.  GCHF(1) denotes a left-handed coupling, and GCHF(2) a right-
C handed coupling.                                                      
C                                                                       
C INPUT:                                                                
C       real    SW2            : square of sine of the weak angle       
C       real    ZMASS          : Z       mass                           
C       real    FMASS          : fermion mass                           
C                                                                       
C OUTPUT:                                                               
C       complex GCHF(2)        : coupling of fermion and Higgs          
C
      IMPLICIT NONE
      COMPLEX*16 GCHF(2)
      REAL*8    SW2,ZMASS,FMASS,ALPHA,FOURPI,EZ,G
C
      ALPHA=1.D0/128.D0
C      ALPHA=1./REAL(137.0359895)
      FOURPI=4.D0*3.14159265358979323846D0
      EZ=SQRT(ALPHA*FOURPI)/SQRT(SW2*(1.D0-SW2))
      G=EZ*FMASS*0.5D0/ZMASS
C
      GCHF(1) = DCMPLX( -G )
      GCHF(2) = DCMPLX( -G )
C
      RETURN
      END

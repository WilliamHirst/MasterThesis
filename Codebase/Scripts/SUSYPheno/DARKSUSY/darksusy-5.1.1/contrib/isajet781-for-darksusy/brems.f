C-----------------------------------------------------------------
      SUBROUTINE BREMS(AEMB,ASMB,MS,MC,MB,Q,CI,GBREM)            
C-----------------------------------------------------------------
C
C    Computes bremsstrahlung corrections according to formulae
C    given in Appendix B of Greub.
C
      IMPLICIT NONE
      REAL*8 AEMB,ASMB,MB,MC,MS,Q,CI(8)
      REAL*8 GBREM
      EXTERNAL FNTG
      COMMON /TAUS/ FAC22,FAC27,FAC28,FAC78
      REAL*8 FAC22,FAC27,FAC28,FAC78
      SAVE /TAUS/ 
      REAL*8 PI,FACAL,FACBL,RINT,BREMSF,BREMS88
      REAL*8 DTRINT
      
      PI=4.d0*ATAN(1.d0)          
C     
      FACAL=AEMB*ASMB/768.d0/PI**5/4.d0
      FACBL=AEMB*ASMB/96.d0/PI**5
      FAC22= 2.*4./9.d0*CI(2)**2
      FAC27=-32.*2./3.d0*CI(2)*CI(7)
      FAC28=-32.*2./3.d0*(-1./3.d0)*CI(2)*CI(8)
      FAC78=-128.d0*(-1./3.d0)*CI(7)*CI(8)
C
      RINT=DTRINT(FNTG,0,64,1.d-3,0.d0,1.d0,1.d0,0.d0,1.d0,1.d0)
      BREMSF=FACAL*RINT
      BREMS88=FACBL*(-1./3.d0*CI(8))**2
     $        *(16./3.d0-4.*PI**2/3.d0+4.d0*LOG(MB/Q))
C
      GBREM=BREMSF+BREMS88
      
      RETURN
      END 

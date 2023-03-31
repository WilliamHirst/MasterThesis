C---------------------------------------------------------------------- 
      REAL*8 FUNCTION FNTG(X,Y)            
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,Y
      COMMON /TAUS/ FAC22,FAC27,FAC28,FAC78
      REAL*8 FAC22,FAC27,FAC28,FAC78
      SAVE /TAUS/ 
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      REAL*8 XKARE,XKAIM
      
      REAL*8 T,U,S,SC,RHO,O22,O27,O28,O78
      
      RHO=(MS/MB)**2
      T= X
      U= Y
      S=T+U-1.d0+RHO
      SC=S*(MB/MC)**2
C
      O22=(XKARE(SC)**2+XKAIM(SC)**2)*(1.d0-S)
      O27=1./2.d0*S*XKARE(SC)
      O28=1./2.d0*S*XKARE(SC)
      O78=1./2.d0*S*(1.d0+1./2.d0*T*U)/T/U
      
      FNTG=FAC22*O22+FAC27*O27+FAC28*O28+FAC78*O78

      RETURN 
      END

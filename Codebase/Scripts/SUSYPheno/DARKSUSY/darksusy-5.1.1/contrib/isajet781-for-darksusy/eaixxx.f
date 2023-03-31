C
C ======================================================================
C
      SUBROUTINE EAIXXX(EB,EA,SHLF,CHLF,PHI,NHE,NHA , EAI)
C
C This subroutine computes an off-shell electron wavefunction after     
C emitting a photon from the electron beam, with a special care for the 
C small angle region.  The momenta are measured in the laboratory frame,
C where the e- beam is along the positive z axis.                       
C                                                                       
C INPUT:                                                                
C       real    EB             : energy (GeV)    of beam  e-            
C       real    EA             : energy (GeV)    of final photon        
C       real    SHLF           : sin(theta/2)    of final photon        
C       real    CHLF           : cos(theta/2)    of final photon        
C       real    PHI            : azimuthal angle of final photon        
C       integer NHE  = -1 or 1 : helicity        of beam  e-            
C       integer NHA  = -1 or 1 : helicity        of final photon        
C                                                                       
C OUTPUT:                                                               
C       complex EAI(6)         : off-shell electron             |e',A,e>
C
      IMPLICIT NONE
      COMPLEX*16 EAI(6),PHS
      REAL*8  EB,EA,SHLF,CHLF,PHI,ME,ALPHA,GAL,RNHE,X,C,S,D,COEFF,
     &        XNNP,XNNM,SNP,CSP
      INTEGER NHE,NHA,NN
C
      ME   = 0.51099906D-3
      ALPHA=1./128.
      GAL  =SQRT(ALPHA*4.*3.14159265D0)
C
      NN=NHA*NHE
      RNHE=NHE
      X=EA/EB
      C=(CHLF+SHLF)*(CHLF-SHLF)
      S=2.*CHLF*SHLF
      D=-1./(EA*EB*(4.*SHLF**2+(ME/EB)**2*C))
      COEFF=-NN*GAL*SQRT(EB)*D
      XNNP=X*(1+NN)
      XNNM=X*(1-NN)
      SNP=SIN(PHI)
      CSP=COS(PHI)
      PHS=DCMPLX( CSP , RNHE*SNP )
C
      EAI((5-3*NHE)/2) = -RNHE*COEFF*ME*S*(1.+XNNP*.5)
      EAI((5-NHE)/2)   =  XNNP*COEFF*ME*CHLF**2*PHS
      EAI((5+NHE)/2)   =  RNHE*COEFF*EB*S*(-2.+XNNM)
      EAI((5+3*NHE)/2) =  XNNM*COEFF*EB*SHLF**2*PHS*2.
C
      EAI(5) =  EB*DCMPLX( 1.-X , 1.-X*C )
      EAI(6) = -EB*X*S*DCMPLX( CSP , SNP )
C
      RETURN
      END

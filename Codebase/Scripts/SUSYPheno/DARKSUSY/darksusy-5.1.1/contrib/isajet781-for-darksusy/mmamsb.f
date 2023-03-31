CDECK  ID>, MMAMSB.
      SUBROUTINE MMAMSB(M0,MHF,G)
C...In this MMAMSB model, M0 is alpha, MHF is the gravitino mass. 
c...  Formulas are from hep-ph/0507110 
c...( Lebedev et. al.)
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS 
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
      REAL MHF,M0
      REAL*8 G(31)
      REAL*8 MS,DPI,AL
      INTEGER I
      real*8 mwhu,mwhd,mwql,mwur,mwdr,mwll,mwer,gk1,gk2,gk3
     &,atop,abot,atau
C...Notice for gaugino masses, the signs are different 
C   from ph/0507110 (Lebedev). 
C
        mwql=XAMIN(1)
        mwdr=XAMIN(2)
        mwur=XAMIN(3)
        mwll=XAMIN(4)
        mwer=XAMIN(5)
        mwhd=XAMIN(6)
        mwhu=XAMIN(7)
        gk1=XAMIN(8)
        gk2=XAMIN(9)
        gk3=XAMIN(10)
        DPI=4.D0*DATAN(1.D0)
        AL=M0
        MS=DBLE(MHF)/16.D0/DPI/DPI        
        atop=3.-mwhu-mwql-mwur
        abot=3.-mwhd-mwql-mwdr
        atau=3.-mwhd-mwll-mwer
        G(7)=G(7)+MS*AL*gk1
        G(8)=G(8)+MS*AL*gk2
        G(9)=G(9)+MS*AL*gk3
        G(10)=G(10)-MS*AL*atau
        G(11)=G(11)-MS*AL*abot
        G(12)=G(12)-MS*AL*atop
        G(13)=G(13)+MS*MS*(AL*AL*(1.-mwhd)+2*AL*(-3.D0/2.D0*G(2)**2*gk2
     ,-3.D0/10.D0*G(1)**2*gk1+3.D0*G(5)**2*abot+1.D0*G(4)**2*atau))
        G(14)=G(14)+MS*MS*(AL*AL*(1.-mwhu)+2*AL*(-3.D0/2.D0*G(2)**2*gk2
     ,-3.D0/10.D0*G(1)**2*gk1+3.D0*G(6)**2*atop))
        G(15)=G(15)+MS*MS*(AL*AL*(1.-mwer)
     ,+2*AL*(-6.D0/5.D0*G(1)**2*gk1))
        G(16)=G(16)+MS*MS*(AL*AL*(1.-mwll)+2*AL*(-3.D0/2.D0*G(2)**2*gk2
     ,-3.D0/10.D0*G(1)**2*gk1))
        G(17)=G(17)+MS*MS*(AL*AL*(1.-mwdr)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-2.D0/15.D0*G(1)**2*gk1))
        G(18)=G(18)+MS*MS*(AL*AL*(1.-mwur)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-8.D0/15.D0*G(1)**2*gk1))
        G(19)=G(19)+MS*MS*(AL*AL*(1.-mwql)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-3.D0/2.D0*G(2)**2*gk2-1.D0/30.D0*G(1)**2*gk1))
        G(20)=G(20)+MS*MS*(AL*AL*(1.-mwer)+2*AL*(-6.D0/5.D0*G(1)**2*gk1
     ,+2.D0*G(4)**2*atau))
        G(21)=G(21)+MS*MS*(AL*AL*(1.-mwll)+2*AL*(-3.D0/2.D0*G(2)**2*gk2
     ,-3.D0/10.D0*G(1)**2*gk1+1.D0*G(4)**2*atau))
        G(22)=G(22)+MS*MS*(AL*AL*(1.-mwdr)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-2.D0/15.D0*G(1)**2*gk1+2.D0*G(5)**2*abot))
        G(23)=G(23)+MS*MS*(AL*AL*(1.-mwur)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-8.D0/15.D0*G(1)**2*gk1+2.D0*G(6)**2*atop))
        G(24)=G(24)+MS*MS*(AL*AL*(1.-mwql)+2*AL*(-8.D0/3.D0*G(3)**2*gk3
     ,-3.D0/2.D0*G(2)**2*gk2-1.D0/30.D0*G(1)**2*gk1
     ,+1.D0*G(6)**2*atop+1.D0*G(5)**2*abot))
c        WRITE(6,*) 'g7-12=',g(7),g(8),g(9),g(10),g(11),g(12)
c        WRITE(6,*) 'g13-24=',sign(1,g(13))*sqrt(abs(g(13))),
c     $sign(1,g(14))*sqrt(abs(g(14))),sign(1,g(15))*sqrt(abs(g(15))),
c     $sign(1,g(16))*sqrt(abs(g(16))),sign(1,g(17))*sqrt(abs(g(17))),
c     $sign(1,g(18))*sqrt(abs(g(18))),sign(1,g(19))*sqrt(abs(g(19))),
c     $sign(1,g(20))*sqrt(abs(g(20))),sign(1,g(21))*sqrt(abs(g(21))),
c     $sign(1,g(22))*sqrt(abs(g(22))),sign(1,g(23))*sqrt(abs(g(23))),
c     $sign(1,g(24))*sqrt(abs(g(24)))  
      RETURN
      END

* dse_h0HH.F
* generated 30-May-2012 12:06
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"
#include "types.h"


	ComplexType function dse_h0HH(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Ind3, Ind4
	integer Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup2
	ComplexType dup3, dup4, dup5, dup6, dup7, dup8, dup9, tmp1
	ComplexType tmp2

	dse_h0HH = 0

	LOOP(Gen3, gM,3,1)

        dse_h0HH = dse_h0HH - 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(tM,Gen3)*
     &        (B1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)) + 
     &          k2*DB1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)) + 
     &          2*DB0(k2,Mf2(tM,Gen3),Mf2(tM,Gen3))*Mf2(tM,Gen3)))/
     &      (MW2*SB2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp1 = DB0(k2,MASf2(All3,tM),MASf2(All4,tM))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_h0HH = dse_h0HH - 
     &    1/(48.D0*Pi)*(Alfa1L*tmp1*
     &        ((Delta(Ind1,Ind2)*
     &              (MW*MZ*SAB*SB*(-3 + 4*SW2) + 
     &                6*CA*CW*Mf2(tM,Ind1))*UASf(All4,Ind2,tM) + 
     &             3*CW*(CA*Kf(Ind1,Ind2,tM) + 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(tM,Ind1))*
     &              UASf(All4,3 + Ind2,tM))*UASfC(All3,Ind1,tM) + 
     &          (3*CW*(CA*KfC(Ind2,Ind1,tM) + 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(tM,Ind1))*
     &              UASf(All4,Ind2,tM) + 
     &             2*Delta(Ind1,Ind2)*
     &              (-2*MW*MZ*SAB*SB*SW2 + 3*CA*CW*Mf2(tM,Ind1))*
     &              UASf(All4,3 + Ind2,tM))*UASfC(All3,3 + Ind1,tM)
     &          )*((Delta(Ind3,Ind4)*
     &              (CAB*MW*MZ*SB*(-3 + 4*SW2) - 
     &                6*CW*SA*Mf2(tM,Ind3))*UASf(All3,Ind4,tM) + 
     &             3*CW*(-(SA*Kf(Ind3,Ind4,tM)) + 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(tM,Ind3))*
     &              UASf(All3,3 + Ind4,tM))*UASfC(All4,Ind3,tM) - 
     &          (3*CW*(SA*KfC(Ind4,Ind3,tM) - 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(tM,Ind3))*
     &              UASf(All3,Ind4,tM) + 
     &             2*Delta(Ind3,Ind4)*
     &              (2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM,Ind3))*
     &              UASf(All3,3 + Ind4,tM))*UASfC(All4,3 + Ind3,tM)
     &          ))/(CW2*MW2*SB2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_h0HH =", dse_h0HH ENDL
#endif

	if( mssmpart .eq. 1 ) return

	LOOP(Gen3, gM,3,1)

        dse_h0HH = dse_h0HH + 
     &    3/(4.D0*Pi)*(Alfa1L*CA*SA*Mf2(bX,Gen3)*
     &        (B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          k2*DB1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          2*DB0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CB2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = DB0(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_h0HH = dse_h0HH - 
     &    1/(48.D0*Pi)*(Alfa1L*tmp2*
     &        ((Delta(Ind3,Ind4)*
     &              (CB*MW*MZ*SAB*(-3 + 2*SW2) + 
     &                6*CW*SA*Mf2(bX,Ind3))*UASf(All4,Ind4,bX) + 
     &             3*CW*(SA*Kf(Ind3,Ind4,bX) + 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,Ind3,bX) + 
     &          (3*CW*(SA*KfC(Ind4,Ind3,bX) + 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,Ind4,bX) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CB*MW*MZ*SAB*SW2 - 3*CW*SA*Mf2(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,3 + Ind3,bX)
     &          )*((Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*(-3 + 2*SW2) + 
     &                6*CA*CW*Mf2(bX,Ind1))*UASf(All3,Ind2,bX) + 
     &             3*CW*(CA*Kf(Ind1,Ind2,bX) - 
     &                MUEC*SA*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All3,3 + Ind2,bX))*UASfC(All4,Ind1,bX) + 
     &          (3*CW*(CA*KfC(Ind2,Ind1,bX) - 
     &                MUE*SA*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &              UASf(All3,Ind2,bX) - 
     &             2*Delta(Ind1,Ind2)*
     &              (CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Ind1))*
     &              UASf(All3,3 + Ind2,bX))*UASfC(All4,3 + Ind1,bX)
     &          ))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_h0HH =", dse_h0HH ENDL
#endif

	if( mssmpart .eq. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_h0HH = dse_h0HH + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (-(CAB*CB2*MW2*MZ2*SAB*
     &             DB0(k2,MSf2(1,1,Gen3),MSf2(1,1,Gen3))) + 
     &          4*CA*CW2*SA*Mf2(2,Gen3)*
     &           (B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*DB1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             2*DB0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3))))
     &       /(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        dse_h0HH = dse_h0HH - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,1,2,Gen3)*
     &           ((CB*MW*MZ*SAB*(-1 + 2*SW2) + 
     &                2*CW*SA*Mf2(2,Gen3))*USfC(Sfe3,1,2,Gen3) + 
     &             CW*(SA*KfC(Gen3,Gen3,2) + CA*MUE*Mf(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)) + 
     &          USf(Sfe4,2,2,Gen3)*
     &           (CW*(SA*Kf(Gen3,Gen3,2) + CA*MUEC*Mf(2,Gen3))*
     &              USfC(Sfe3,1,2,Gen3) + 
     &             2*(-(CB*MW*MZ*SAB*SW2) + CW*SA*Mf2(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)))*
     &        (USf(Sfe3,1,2,Gen3)*
     &           ((CAB*CB*MW*MZ*(-1 + 2*SW2) + 
     &                2*CA*CW*Mf2(2,Gen3))*USfC(Sfe4,1,2,Gen3) + 
     &             CW*(CA*KfC(Gen3,Gen3,2) - MUE*SA*Mf(2,Gen3))*
     &              USfC(Sfe4,2,2,Gen3)) + 
     &          USf(Sfe3,2,2,Gen3)*
     &           (CW*(CA*Kf(Gen3,Gen3,2) - MUEC*SA*Mf(2,Gen3))*
     &              USfC(Sfe4,1,2,Gen3) + 
     &             2*(-(CAB*CB*MW*MZ*SW2) + CA*CW*Mf2(2,Gen3))*
     &              USfC(Sfe4,2,2,Gen3))))/(CB2*CW2*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_h0HH =", dse_h0HH ENDL
#endif

	if( mssmpart .eq. 3 ) return

        dse_h0HH = dse_h0HH - 
     &    1/(32.D0*Pi)*(Alfa1L*
     &        (-2*(-(CAB*MW2*S2B**2*SAB) + 
     &             CBA*CW2*(k2*SBA + MA02*SBA))*DB0(k2,MA02,MZ2) + 
     &          MW2*(SAB*(C2B**2*CAB*DB0(k2,MA02,MA02) + 
     &                3*C2A**2*CAB*DB0(k2,Mh02,Mh02) - 
     &                6*C2A*S2A*SAB*DB0(k2,Mh02,Mh02)) + 
     &             (2*CAB*S2A + C2A*SAB)*
     &              (4*S2A*SAB*DB0(k2,Mh02,MHH2) + 
     &                C2A*CAB*
     &                 (-2*DB0(k2,Mh02,MHH2) + 3*DB0(k2,MHH2,MHH2))
     &                ) - 2*(-(C2B*CAB) + 2*CBA*CW2)*
     &              (C2B*SAB + 2*CW2*SBA)*DB0(k2,MHp2,MHp2)) - 
     &          4*(CAB*MW2*S2B*(CBA*CW2 - S2B*SAB) + 
     &             CW2*(CBA*CW2*(k2 + MHp2 - MW2) + MW2*S2B*SAB)*
     &              SBA)*DB0(k2,MHp2,MW2) + 
     &          (C2B**2*CAB*MW2*SAB + 
     &             2*CBA*(-7*MW2 + CW2*(k2 + MZ2))*SBA)*
     &           DB0(k2,MZ2,MZ2) + 
     &          2*((C2B**2*CAB*MW2*SAB + 
     &                CBA*CW2**2*(2*k2*SBA - 12*MW2*SBA))*
     &              DB0(k2,MW2,MW2) + 
     &             CBA*CW2*SBA*
     &              (-B0(k2,MA02,MZ2) + B0(k2,MZ2,MZ2) + 
     &                2*(B1(k2,MA02,MZ2) - B1(k2,MZ2,MZ2) + 
     &                   k2*DB1(k2,MA02,MZ2) + 
     &                   CW2*
     &                    (-B0(k2,MHp2,MW2) + B0(k2,MW2,MW2) + 
     &                      2*
     &                       (B1(k2,MHp2,MW2) + 
     &                       k2*DB1(k2,MHp2,MW2)) - 
     &                      2*(B1(k2,MW2,MW2) + k2*DB1(k2,MW2,MW2))
     &                      )) - 2*k2*DB1(k2,MZ2,MZ2)))))/
     &      (CW2**2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     &    CA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     &    SA*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     &    SA*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     &    CA*UChaC(Cha3,1)*VChaC(Cha4,2)

        dse_h0HH = dse_h0HH + 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (B1(k2,MCha2(Cha3),MCha2(Cha4)) + 
     &             k2*DB1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          DB0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup3*dup4*MCha(Cha3)*MCha(Cha4) + 
     &             dup2*dup4*MCha2(Cha3) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup7 = -(SW*ZNeu(Neu4,1)) + CW*ZNeu(Neu4,2)

	dup8 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)

	dup9 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup10 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup11 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup12 = -(SW*ZNeuC(Neu4,1)) + CW*ZNeuC(Neu4,2)

	dup13 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)

	dup14 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

        dse_h0HH = dse_h0HH + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (-(k2*DB1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &             ((-(dup5*dup8) + 
     &                  dup7*(SA*ZNeu(Neu3,3) + CA*ZNeu(Neu3,4)))*
     &                (dup10*dup14 + CA*dup11*ZNeuC(Neu3,3) + 
     &                  dup12*SA*ZNeuC(Neu3,4)) - 
     &               (dup5*dup9 + CA*dup6*ZNeu(Neu3,3) + 
     &                  dup7*SA*ZNeu(Neu3,4))*
     &                (dup10*dup13 + 
     &                  dup11*(SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4))
     &                  ))) + 
     &          B1(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           ((dup5*dup8 - 
     &                dup7*(SA*ZNeu(Neu3,3) + CA*ZNeu(Neu3,4)))*
     &              (dup10*dup14 + CA*dup11*ZNeuC(Neu3,3) + 
     &                dup12*SA*ZNeuC(Neu3,4)) + 
     &             (dup5*dup9 + CA*dup6*ZNeu(Neu3,3) + 
     &                dup7*SA*ZNeu(Neu3,4))*
     &              (dup10*dup13 + 
     &                dup11*(SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4))))
     &            + DB0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (-(MNeu2(Neu3)*
     &                ((-(dup5*dup8) + dup7*SA*ZNeu(Neu3,3) + 
     &                     CA*dup7*ZNeu(Neu3,4))*
     &                   (dup10*dup14 + CA*dup11*ZNeuC(Neu3,3) + 
     &                     dup12*SA*ZNeuC(Neu3,4)) - 
     &                  (dup5*dup9 + CA*dup6*ZNeu(Neu3,3) + 
     &                     dup7*SA*ZNeu(Neu3,4))*
     &                   (dup10*dup13 + 
     &                     dup11*
     &                      (SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4))))
     &                ) - MNeu(Neu3)*MNeu(Neu4)*
     &              ((-(dup5*dup8) + dup7*SA*ZNeu(Neu3,3) + 
     &                   CA*dup7*ZNeu(Neu3,4))*
     &                 (dup5*dup9 + CA*dup6*ZNeu(Neu3,3) + 
     &                   dup7*SA*ZNeu(Neu3,4)) - 
     &                (dup10*dup14 + CA*dup11*ZNeuC(Neu3,3) + 
     &                   dup12*SA*ZNeuC(Neu3,4))*
     &                 (dup10*dup13 + 
     &                   dup11*
     &                    (SA*ZNeuC(Neu3,3) + CA*ZNeuC(Neu3,4))))))
     &        )/(CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_h0HH =", dse_h0HH ENDL
#endif

	end


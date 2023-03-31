* dse_ZZ.F
* generated 30-May-2012 12:16
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"
#include "types.h"


	ComplexType function dse_ZZ(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Neu3, Neu4
	integer Sfe3, Sfe4

	ComplexType dup1, dup2, dup3, dup4, dup5, dup6, tmp1, tmp2

	dse_ZZ = 0

	LOOP(Gen3, gM,3,1)

        dse_ZZ = dse_ZZ + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        (2*(9 + SW2*(-24 + 32*SW2))*
     &           DB00(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)) + 
     &          (-9 + (24 - 32*SW2)*SW2)*
     &           (B1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)) + 
     &             k2*DB1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3))) - 
     &          9*DB0(k2,Mf2(tM,Gen3),Mf2(tM,Gen3))*Mf2(tM,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp1 = DB00(k2,MASf2(All3,tM),MASf2(All4,tM))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_ZZ = dse_ZZ - 
     &    1/(12.D0*Pi)*(Alfa1L*tmp1*
     &        (-3*UASf(All4,Ind2,tM)*UASfC(All3,Ind2,tM) + 
     &          SW2*(4*UASf(All4,Ind2,tM)*UASfC(All3,Ind2,tM) + 
     &             4*UASf(All4,3 + Ind2,tM)*UASfC(All3,3 + Ind2,tM)
     &             ))*(-3*UASf(All3,Ind1,tM)*UASfC(All4,Ind1,tM) + 
     &          SW2*(4*UASf(All3,Ind1,tM)*UASfC(All4,Ind1,tM) + 
     &             4*UASf(All3,3 + Ind1,tM)*UASfC(All4,3 + Ind1,tM)
     &             )))/(CW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_ZZ =", dse_ZZ ENDL
#endif

	if( mssmpart .eq. 1 ) return

	LOOP(Gen3, gM,3,1)

        dse_ZZ = dse_ZZ + 
     &    1/(24.D0*Pi)*(Alfa1L*
     &        (2*(9 + SW2*(-12 + 8*SW2))*
     &           DB00(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          (-9 + (12 - 8*SW2)*SW2)*
     &           (B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &             k2*DB1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))) - 
     &          9*DB0(k2,Mf2(bX,Gen3),Mf2(bX,Gen3))*Mf2(bX,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = DB00(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dse_ZZ = dse_ZZ - 
     &    1/(12.D0*Pi)*(Alfa1L*tmp2*
     &        (-3*UASf(All4,Ind2,bX)*UASfC(All3,Ind2,bX) + 
     &          SW2*(2*UASf(All4,Ind2,bX)*UASfC(All3,Ind2,bX) + 
     &             2*UASf(All4,3 + Ind2,bX)*UASfC(All3,3 + Ind2,bX)
     &             ))*(-3*UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) + 
     &          SW2*(2*UASf(All3,Ind1,bX)*UASfC(All4,Ind1,bX) + 
     &             2*UASf(All3,3 + Ind1,bX)*UASfC(All4,3 + Ind1,bX)
     &             )))/(CW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_ZZ =", dse_ZZ ENDL
#endif

	if( mssmpart .eq. 2 ) return

        dse_ZZ = dse_ZZ - 
     &    3/(8.D0*Pi)*(Alfa1L*
     &        (B1(k2,0.D0,0.D0) - 2*DB00(k2,0.D0,0.D0) + k2*DB1(k2,0.D0,0.D0)))
     &       /(CW2*SW2)

	LOOP(Gen3, 1,3,1)

        dse_ZZ = dse_ZZ - 
     &    1/(8.D0*Pi)*(Alfa1L*
     &        (2*DB00(k2,MSf2(1,1,Gen3),MSf2(1,1,Gen3)) + 
     &          (1 + SW2*(-4 + 8*SW2))*
     &           (B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) - 
     &             2*DB00(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &             k2*DB1(k2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     &          DB0(k2,Mf2(2,Gen3),Mf2(2,Gen3))*Mf2(2,Gen3)))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        dse_ZZ = dse_ZZ - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        DB00(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (-(USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     &          SW2*(2*USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     &             2*USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))*
     &        (-(USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     &          SW2*(2*USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &             2*USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3))))/
     &      (CW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_ZZ =", dse_ZZ ENDL
#endif

	if( mssmpart .eq. 3 ) return

        dse_ZZ = dse_ZZ - 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        (-(MW2*(SBA2*DB0(k2,Mh02,MZ2) + 
     &               CBA2*DB0(k2,MHH2,MZ2))) + 
     &          CW2*((CW2*(5*CW2*k2 + 2*MW2) - 2*MW2*SW2)*
     &              DB0(k2,MW2,MW2) + 
     &             SBA2*(DB00(k2,Mh02,MZ2) + DB00(k2,MHH2,MA02)) + 
     &             CBA2*(DB00(k2,Mh02,MA02) + DB00(k2,MHH2,MZ2)) + 
     &             (CW2 - SW2)**2*DB00(k2,MHp2,MHp2) + 
     &             (CW2*(9*CW2 - 2*SW2) + SW2**2)*
     &              DB00(k2,MW2,MW2) + 
     &             CW2**2*(5*B0(k2,MW2,MW2) + 
     &                2*(B1(k2,MW2,MW2) + k2*DB1(k2,MW2,MW2))))))/
     &      (CW2**2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2)) - 
     &    UCha(Cha4,1)*UChaC(Cha3,1)

        dup2 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     &    UCha(Cha3,1)*UChaC(Cha4,1)

        dup3 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     &    VCha(Cha4,1)*VChaC(Cha3,1)

        dup4 = SW2*Delta(Cha3,Cha4) - 
     &    1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2)) - 
     &    VCha(Cha3,1)*VChaC(Cha4,1)

        dse_ZZ = dse_ZZ + 
     &    1/(2.D0*Pi)*(Alfa1L*
     &        (-((dup1*dup2 + dup3*dup4)*
     &             (B1(k2,MCha2(Cha3),MCha2(Cha4)) - 
     &               2*DB00(k2,MCha2(Cha3),MCha2(Cha4)) + 
     &               k2*DB1(k2,MCha2(Cha3),MCha2(Cha4)))) + 
     &          DB0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup2*dup4*MCha(Cha3)*MCha(Cha4) - 
     &             dup3*dup4*MCha2(Cha3) + 
     &             dup1*(dup3*MCha(Cha3)*MCha(Cha4) - 
     &                dup2*MCha2(Cha3)))))/(CW2*SW2)

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - ZNeu(Neu4,4)*ZNeuC(Neu3,4)

	dup6 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - ZNeu(Neu3,4)*ZNeuC(Neu4,4)

        dse_ZZ = dse_ZZ + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (-(DB0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &             (dup6**2*MNeu(Neu3)*MNeu(Neu4) + 
     &               dup5*(dup5*MNeu(Neu3)*MNeu(Neu4) + 
     &                  2*dup6*MNeu2(Neu3)))) + 
     &          dup5*(4*dup6*DB00(k2,MNeu2(Neu3),MNeu2(Neu4)) + 
     &             2*(B1(k2,MNeu2(Neu3),MNeu2(Neu4)) + 
     &                k2*DB1(k2,MNeu2(Neu3),MNeu2(Neu4)))*
     &              (-(ZNeu(Neu3,3)*ZNeuC(Neu4,3)) + 
     &                ZNeu(Neu3,4)*ZNeuC(Neu4,4)))))/(CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_ZZ =", dse_ZZ ENDL
#endif

	end


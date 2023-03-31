
	C8LNeu = 0

	LOOP(All4, 1,6,1)

	tmp1 = A0(MASf2(All4,bTR))

	LOOP(Neu4, 1,4,1)

	dup1 = CW*SW*ZNeu(Neu4,1) - 3*CW2*ZNeu(Neu4,2)

	dup2 = SW2*ZNeu(Neu4,1) - 3*CW*SW*ZNeu(Neu4,2)

        C8LNeu = C8LNeu - 
     &   1/18.D0*(MW2*UASfC(All4,2,bTR)*
     &       (1/2.D0*((-2*tmp1 + 
     &              MB2*(1 + 
     &                 2*C0(MB2,0.D0,0.D0,MASf2(All4,bTR),MNeu2(Neu4),
     &                   MASf2(All4,bTR))*MASf2(All4,bTR)) + 
     &              2*(A0(MNeu2(Neu4)) + 
     &                 B0(MB2,MASf2(All4,bTR),MNeu2(Neu4))*
     &                  (MASf2(All4,bTR) - MNeu2(Neu4))))*
     &            (CB*MW*UASf(All4,3,bTR)*
     &               (-3*ZNeu(Neu4,2)*
     &                  (CW*SW*ZNeuC(Neu4,1) - 3*CW2*ZNeuC(Neu4,2))
     &                   + ZNeu(Neu4,1)*
     &                  (SW2*ZNeuC(Neu4,1) - 3*CW*SW*ZNeuC(Neu4,2))
     &                 ) + 
     &              3*dup1*Mf(bTR,3)*UASf(All4,6,bTR)*ZNeuC(Neu4,3)
     &              )) + MB*MNeu(Neu4)*
     &          (B0(MB2,MASf2(All4,bTR),MNeu2(Neu4))*
     &             (-2*CB*dup2*MW*UASf(All4,6,bTR)*ZNeu(Neu4,1) - 
     &               3*dup1*Mf(bTR,3)*UASf(All4,3,bTR)*ZNeu(Neu4,3)
     &               ) + ((tmp1 - A0(MNeu2(Neu4)))*
     &               (2*CB*dup2*MW*UASf(All4,6,bTR)*ZNeu(Neu4,1) + 
     &                 3*dup1*Mf(bTR,3)*UASf(All4,3,bTR)*
     &                  ZNeu(Neu4,3)))/
     &             (MASf2(All4,bTR) - MNeu2(Neu4)))))/
     &     (CB*CW2*MB2**2*MW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Neu4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C8LNeu =", C8LNeu ENDL
#endif

	C8RNeu = 0

	LOOP(All4, 1,6,1)

	tmp2 = A0(MASf2(All4,bTR))

	LOOP(Neu4, 1,4,1)

        dup3 = 2*CB*MW*SW*UASf(All4,6,bTR)*ZNeu(Neu4,1) + 
     &   3*CW*Mf(bTR,3)*UASf(All4,3,bTR)*ZNeu(Neu4,3)

        C8RNeu = C8RNeu - 
     &   1/9.D0*(MW2*SW2*UASfC(All4,5,bTR)*ZNeuC(Neu4,1)*
     &       (1/2.D0*(dup3*(-2*tmp2 + 2*A0(MNeu2(Neu4)) + 
     &              MB2*(1 + 
     &                 2*C0(MB2,0.D0,0.D0,MASf2(All4,bTR),MNeu2(Neu4),
     &                   MASf2(All4,bTR))*MASf2(All4,bTR)) + 
     &              B0(MB2,MASf2(All4,bTR),MNeu2(Neu4))*
     &               (2*MASf2(All4,bTR) - 2*MNeu2(Neu4)))) + 
     &         MB*MNeu(Neu4)*
     &          (-B0(MB2,MASf2(All4,bTR),MNeu2(Neu4)) + 
     &            (tmp2 - A0(MNeu2(Neu4)))/
     &             (MASf2(All4,bTR) - MNeu2(Neu4)))*
     &          (CB*MW*UASf(All4,3,bTR)*
     &             (SW*ZNeuC(Neu4,1) - 3*CW*ZNeuC(Neu4,2)) + 
     &            3*CW*Mf(bTR,3)*UASf(All4,6,bTR)*ZNeuC(Neu4,3))))/
     &     (CB*CW2*MB2**2*MW*SW*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Neu4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C8RNeu =", C8RNeu ENDL
#endif


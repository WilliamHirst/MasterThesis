
	C7LCha = 0

	LOOP(All4, 1,6,1)

	tmp5 = A0(MASf2(All4,3))

	LOOP(Cha4, 1,2,1)

	tmp1 = A0(MCha2(Cha4))

	tmp2 = B0(MB2,MASf2(All4,3),MCha2(Cha4))

	tmp3 = C0(MB2,0.D0,0.D0,MASf2(All4,3),MCha2(Cha4),MASf2(All4,3))

	tmp4 = C0(MB2,0.D0,0.D0,MASf2(All4,3),MCha2(Cha4),MCha2(Cha4))

	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        dup1 = sqrt2*(MW*SB*UASfC(All4,Ind1,3)*VCha(Cha4,1)) - 
     &   Mf(3,Ind1)*UASfC(All4,3 + Ind1,3)*VCha(Cha4,2)

        C7LCha = C7LCha - 
     &   1/6.D0*(MB*SB2*CKM(Ind2,3)*CKMC(Ind1,2)*MCha(Cha4)*
     &        Mf(bTR,3)*UASf(All4,Ind2,3)*UCha(Cha4,2)*
     &        ((dup1*(-tmp1 + tmp5))/
     &           (MASf2(All4,3) - MCha2(Cha4)) + 
     &          tmp2*(-(sqrt2*
     &                (MW*SB*UASfC(All4,Ind1,3)*VCha(Cha4,1))) + 
     &             Mf(3,Ind1)*UASfC(All4,3 + Ind1,3)*VCha(Cha4,2)))
     &         + MB2*(1/2.D0*
     &           (CB*CKM(Ind1,3)*CKMC(Ind2,2)*
     &              (-2*tmp1 + MB2*(5 + 4*tmp3*MASf2(All4,3)) + 
     &                2*(tmp5 + 
     &                   tmp2*(-MASf2(All4,3) + MCha2(Cha4))))*
     &              (SB2*UASf(All4,Ind1,3)*
     &                 (-(sqrt2*
     &                      (MW*Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*
     &                       VCha(Cha4,2))) + 
     &                   2*MW2*SB*UASfC(All4,Ind2,3)*VCha(Cha4,1))*
     &                 VChaC(Cha4,1) + 
     &                Mf(3,Ind1)*UASf(All4,3 + Ind1,3)*
     &                 (-(sqrt2*
     &                      (MW*SB2*UASfC(All4,Ind2,3)*
     &                       VCha(Cha4,1))) + 
     &                   SB*Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*
     &                    VCha(Cha4,2))*VChaC(Cha4,2)))/MB2 - 
     &          3*tmp4*(dup1*MB*SB2*CKM(Ind2,3)*CKMC(Ind1,2)*
     &              MCha(Cha4)*Mf(bTR,3)*UASf(All4,Ind2,3)*
     &              UCha(Cha4,2) + 
     &             CB*CKM(Ind1,3)*CKMC(Ind2,2)*MCha2(Cha4)*
     &              (SB2*UASf(All4,Ind1,3)*
     &                 (sqrt2*
     &                    (MW*Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*
     &                      VCha(Cha4,2)) - 
     &                   2*MW2*SB*UASfC(All4,Ind2,3)*VCha(Cha4,1))*
     &                 VChaC(Cha4,1) + 
     &                Mf(3,Ind1)*UASf(All4,3 + Ind1,3)*
     &                 (sqrt2*
     &                    (MW*SB2*UASfC(All4,Ind2,3)*VCha(Cha4,1))-
     &                     SB*Mf(3,Ind2)*UASfC(All4,3 + Ind2,3)*
     &                    VCha(Cha4,2))*VChaC(Cha4,2)))))/
     &     (CB*MB2**2*SB*SB2*CKM(3,3)*CKMC(3,2))

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)

	ENDLOOP(Cha4)

	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DCONST "C7LCha =", C7LCha ENDL
#endif


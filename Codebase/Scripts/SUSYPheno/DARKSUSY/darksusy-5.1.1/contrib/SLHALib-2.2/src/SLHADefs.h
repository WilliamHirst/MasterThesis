#if 0
	SLHADefs.h
		declarations for SLHALib data
		generated 21 Dec 2007 17:40
#endif

#ifndef __SLHADefs_h__
#define __SLHADefs_h__

#define invalid -999

#define OffsetModSel 0
#define LengthModSel 12
#define BlockModSel(i) SlhaData(i)
#define ModSel_Model Slhadata(1)
#define ModSel_Content Slhadata(2)
#define ModSel_RPV Slhadata(3)
#define ModSel_CPV Slhadata(4)
#define ModSel_FV Slhadata(5)
#define ModSel_GridPts Slhadata(6)
#define ModSel_Qmax Slhadata(7)
#define ModSel_PDG(n) Slhadata(7+n)

#define OffsetSMInputs 12
#define LengthSMInputs 16
#define BlockSMInputs(i) SlhaData(12+i)
#define SMInputs_invAlfaMZ Slhadata(13)
#define SMInputs_GF Slhadata(14)
#define SMInputs_AlfasMZ Slhadata(15)
#define SMInputs_MZ Slhadata(16)
#define SMInputs_Mf(t,g) Slhadata(12+t+4*(g))
#define SMInputs_MfFlat(i) Slhadata(16+i)
#define   SMInputs_Mnue SMInputs_Mf(1,1)
#define   SMInputs_Me SMInputs_Mf(2,1)
#define   SMInputs_Mu SMInputs_Mf(3,1)
#define   SMInputs_Md SMInputs_Mf(4,1)
#define   SMInputs_Mnumu SMInputs_Mf(1,2)
#define   SMInputs_Mmu SMInputs_Mf(2,2)
#define   SMInputs_Mc SMInputs_Mf(3,2)
#define   SMInputs_Ms SMInputs_Mf(4,2)
#define   SMInputs_Mnutau SMInputs_Mf(1,3)
#define   SMInputs_Mtau SMInputs_Mf(2,3)
#define   SMInputs_Mt SMInputs_Mf(3,3)
#define   SMInputs_Mb SMInputs_Mf(4,3)

#define OffsetMinPar 28
#define LengthMinPar 6
#define BlockMinPar(i) SlhaData(28+i)
#define MinPar_M0 Slhadata(29)
#define   MinPar_Lambda MinPar_M0
#define MinPar_M12 Slhadata(30)
#define   MinPar_Mmess MinPar_M12
#define   MinPar_M32 MinPar_M12
#define MinPar_TB Slhadata(31)
#define MinPar_signMUE Slhadata(32)
#define MinPar_A Slhadata(33)
#define   MinPar_N5 MinPar_A
#define MinPar_cgrav Slhadata(34)

#define OffsetExtPar 34
#define LengthExtPar 42
#define BlockExtPar(i) SlhaData(34+i)
#define ExtPar_Q SlhaData(35)
#define ExtPar_M1 Slhadata(36)
#define ExtPar_M2 Slhadata(37)
#define ExtPar_M3 Slhadata(38)
#define ExtPar_Af(t) Slhadata(37+t)
#define   ExtPar_Atau ExtPar_Af(2)
#define   ExtPar_At ExtPar_Af(3)
#define   ExtPar_Ab ExtPar_Af(4)
#define ExtPar_MHu2 Slhadata(42)
#define ExtPar_MHd2 Slhadata(43)
#define ExtPar_MUE Slhadata(44)
#define ExtPar_MA02 Slhadata(45)
#define ExtPar_TB Slhadata(46)
#define ExtPar_MA0 Slhadata(47)
#define ExtPar_MHp Slhadata(48)
#define ExtPar_MSS(g,q) Slhadata(45+g+3*(q))
#define   ExtPar_MSL(g) ExtPar_MSS(g,1)
#define   ExtPar_MSE(g) ExtPar_MSS(g,2)
#define   ExtPar_MSQ(g) ExtPar_MSS(g,3)
#define   ExtPar_MSU(g) ExtPar_MSS(g,4)
#define   ExtPar_MSD(g) ExtPar_MSS(g,5)
#define ExtPar_N5(g) Slhadata(63+g)
#define ExtPar_lambda Slhadata(67)
#define ExtPar_kappa Slhadata(68)
#define ExtPar_Alambda Slhadata(69)
#define ExtPar_Akappa Slhadata(70)
#define ExtPar_lambdaS Slhadata(71)
#define ExtPar_xiF Slhadata(72)
#define ExtPar_xiS Slhadata(73)
#define ExtPar_MUEprime Slhadata(74)
#define ExtPar_mS2prime Slhadata(75)
#define ExtPar_mS2 Slhadata(76)

#define OffsetQExtPar 76
#define LengthQExtPar 26
#define BlockQExtPar(i) SlhaData(76+i)
#define QExtPar_QM1 Slhadata(77)
#define QExtPar_QM2 Slhadata(78)
#define QExtPar_QM3 Slhadata(79)
#define QExtPar_QAf(t) Slhadata(78+t)
#define   QExtPar_QAtau QExtPar_QAf(2)
#define   QExtPar_QAt QExtPar_QAf(3)
#define   QExtPar_QAb QExtPar_QAf(4)
#define QExtPar_QMHu2 Slhadata(83)
#define QExtPar_QMHd2 Slhadata(84)
#define QExtPar_QMUE Slhadata(85)
#define QExtPar_QMA02 Slhadata(86)
#define QExtPar_QTB Slhadata(87)
#define QExtPar_QMSS(g,q) Slhadata(84+g+3*(q))
#define   QExtPar_QMSL(g) QExtPar_QMSS(g,1)
#define   QExtPar_QMSE(g) QExtPar_QMSS(g,2)
#define   QExtPar_QMSQ(g) QExtPar_QMSS(g,3)
#define   QExtPar_QMSU(g) QExtPar_QMSS(g,4)
#define   QExtPar_QMSD(g) QExtPar_QMSS(g,5)

#define OffsetNMSSMRun 102
#define LengthNMSSMRun 11
#define BlockNMSSMRun(i) SlhaData(102+i)
#define NMSSMRun_Q SlhaData(103)
#define NMSSMRun_lambda Slhadata(104)
#define NMSSMRun_kappa Slhadata(105)
#define NMSSMRun_Alambda Slhadata(106)
#define NMSSMRun_Akappa Slhadata(107)
#define NMSSMRun_lambdaS Slhadata(108)
#define NMSSMRun_xiF Slhadata(109)
#define NMSSMRun_xiS Slhadata(110)
#define NMSSMRun_MUEprime Slhadata(111)
#define NMSSMRun_mS2prime Slhadata(112)
#define NMSSMRun_mS2 Slhadata(113)

#define OffsetMass 113
#define LengthMass 53
#define BlockMass(i) SlhaData(113+i)
#define Mass_Mf(t,g) Slhadata(109+t+4*(g))
#define Mass_MfFlat(i) Slhadata(113+i)
#define Mass_MSf(s,t,g) Slhadata(115+s+8*(g)+2*(t))
#define Mass_MSfFlat(i) Slhadata(125+i)
#define Mass_MZ Slhadata(150)
#define Mass_MW Slhadata(151)
#define Mass_Mh0 Slhadata(152)
#define Mass_MHH Slhadata(153)
#define Mass_MA0 Slhadata(154)
#define Mass_MHp Slhadata(155)
#define   Mass_MH1 Mass_Mh0
#define   Mass_MH2 Mass_MHH
#define Mass_MH3 Slhadata(156)
#define   Mass_MA1 Mass_MA0
#define Mass_MA2 Slhadata(157)
#define Mass_MNeu(n) Slhadata(157+n)
#define Mass_MCha(c) Slhadata(162+c)
#define Mass_MGl Slhadata(165)
#define Mass_MGrav Slhadata(166)

#define OffsetDMass 166
#define LengthDMass 4
#define BlockDMass(i) SlhaData(166+i)
#define DMass_DeltaMh0 Slhadata(167)
#define DMass_DeltaMHH Slhadata(168)
#define DMass_DeltaMA0 Slhadata(169)
#define DMass_DeltaMHp Slhadata(170)

#define OffsetNMix 170
#define LengthNMix 16
#define BlockNMix(i) SlhaData(170+i)
#define NMix_ZNeu(n1,n2) Slhadata(166+n1+4*(n2))
#define NMix_ZNeuFlat(i) Slhadata(170+i)

#define OffsetUMix 186
#define LengthUMix 4
#define BlockUMix(i) SlhaData(186+i)
#define UMix_UCha(c1,c2) Slhadata(184+c1+2*(c2))
#define UMix_UChaFlat(i) Slhadata(186+i)

#define OffsetVMix 190
#define LengthVMix 4
#define BlockVMix(i) SlhaData(190+i)
#define VMix_VCha(c1,c2) Slhadata(188+c1+2*(c2))
#define VMix_VChaFlat(i) Slhadata(190+i)

#define OffsetSfMix 194
#define LengthSfMix 12
#define BlockSfMix(i) SlhaData(194+i)
#define SfMix_USf(s1,s2,t) Slhadata(184+s1+2*(s2)+4*(t))
#define SfMix_USfFlat(i,t) Slhadata(186+i+4*(t))

#define OffsetStauMix 194
#define LengthStauMix 4
#define BlockStauMix(i) SlhaData(194+i)
#define   StauMix_USf(s1,s2) SfMix_USf(s1,s2,2)
#define   StauMix_USfFlat(i) SfMix_USfFlat(i,2)

#define OffsetStopMix 198
#define LengthStopMix 4
#define BlockStopMix(i) SlhaData(198+i)
#define   StopMix_USf(s1,s2) SfMix_USf(s1,s2,3)
#define   StopMix_USfFlat(i) SfMix_USfFlat(i,3)

#define OffsetSbotMix 202
#define LengthSbotMix 4
#define BlockSbotMix(i) SlhaData(202+i)
#define   SbotMix_USf(s1,s2) SfMix_USf(s1,s2,4)
#define   SbotMix_USfFlat(i) SfMix_USfFlat(i,4)

#define OffsetAlpha 206
#define LengthAlpha 1
#define BlockAlpha(i) SlhaData(206+i)
#define Alpha_Alpha Slhadata(207)

#define OffsetDAlpha 207
#define LengthDAlpha 1
#define BlockDAlpha(i) SlhaData(207+i)
#define DAlpha_DeltaAlpha Slhadata(208)

#define OffsetHMix 208
#define LengthHMix 5
#define BlockHMix(i) SlhaData(208+i)
#define HMix_Q SlhaData(209)
#define HMix_MUE Slhadata(210)
#define HMix_TB Slhadata(211)
#define HMix_VEV Slhadata(212)
#define HMix_MA02 Slhadata(213)

#define OffsetGauge 213
#define LengthGauge 4
#define BlockGauge(i) SlhaData(213+i)
#define Gauge_Q SlhaData(214)
#define Gauge_g1 Slhadata(215)
#define Gauge_g2 Slhadata(216)
#define Gauge_g3 Slhadata(217)

#define OffsetMSoft 217
#define LengthMSoft 21
#define BlockMSoft(i) SlhaData(217+i)
#define MSoft_Q SlhaData(218)
#define MSoft_M1 Slhadata(219)
#define MSoft_M2 Slhadata(220)
#define MSoft_M3 Slhadata(221)
#define MSoft_MHu2 Slhadata(222)
#define MSoft_MHd2 Slhadata(223)
#define MSoft_MSS(g,q) Slhadata(220+g+3*(q))
#define   MSoft_MSL(g) MSoft_MSS(g,1)
#define   MSoft_MSE(g) MSoft_MSS(g,2)
#define   MSoft_MSQ(g) MSoft_MSS(g,3)
#define   MSoft_MSU(g) MSoft_MSS(g,4)
#define   MSoft_MSD(g) MSoft_MSS(g,5)

#define OffsetAf 238
#define LengthAf 30
#define BlockAf(i) SlhaData(238+i)
#define Af_Q(t) Slhadata(219+10*(t))
#define Af_Af(g1,g2,t) Slhadata(216+g1+3*(g2)+10*(t))
#define Af_AfFlat(i,t) Slhadata(219+i+10*(t))

#define OffsetAe 238
#define LengthAe 11
#define BlockAe(i) SlhaData(238+i)
#define   Ae_Q Af_Q(2)
#define   Ae_Af(g1,g2) Af_Af(g1,g2,2)
#define   Ae_AfFlat(i) Af_AfFlat(i,2)
#define   Ae_Atau Ae_Af(3,3)

#define OffsetAu 249
#define LengthAu 11
#define BlockAu(i) SlhaData(249+i)
#define   Au_Q Af_Q(3)
#define   Au_Af(g1,g2) Af_Af(g1,g2,3)
#define   Au_AfFlat(i) Af_AfFlat(i,3)
#define   Au_At Au_Af(3,3)

#define OffsetAd 260
#define LengthAd 11
#define BlockAd(i) SlhaData(260+i)
#define   Ad_Q Af_Q(4)
#define   Ad_Af(g1,g2) Af_Af(g1,g2,4)
#define   Ad_AfFlat(i) Af_AfFlat(i,4)
#define   Ad_Ab Ad_Af(3,3)

#define OffsetYf 271
#define LengthYf 30
#define BlockYf(i) SlhaData(271+i)
#define Yf_Q(t) Slhadata(252+10*(t))
#define Yf_Yf(g1,g2,t) Slhadata(249+g1+3*(g2)+10*(t))
#define Yf_YfFlat(i,t) Slhadata(252+i+10*(t))

#define OffsetYe 271
#define LengthYe 11
#define BlockYe(i) SlhaData(271+i)
#define   Ye_Q Yf_Q(2)
#define   Ye_Yf(g1,g2) Yf_Yf(g1,g2,2)
#define   Ye_YfFlat(i) Yf_YfFlat(i,2)
#define   Ye_Ytau Ye_Yf(3,3)

#define OffsetYu 282
#define LengthYu 11
#define BlockYu(i) SlhaData(282+i)
#define   Yu_Q Yf_Q(3)
#define   Yu_Yf(g1,g2) Yf_Yf(g1,g2,3)
#define   Yu_YfFlat(i) Yf_YfFlat(i,3)
#define   Yu_Yt Yu_Yf(3,3)

#define OffsetYd 293
#define LengthYd 11
#define BlockYd(i) SlhaData(293+i)
#define   Yd_Q Yf_Q(4)
#define   Yd_Yf(g1,g2) Yf_Yf(g1,g2,4)
#define   Yd_YfFlat(i) Yf_YfFlat(i,4)
#define   Yd_Yb Yd_Yf(3,3)

#define OffsetRVLamLLEIn 304
#define LengthRVLamLLEIn 27
#define BlockRVLamLLEIn(i) SlhaData(304+i)
#define RVLamLLEIn_lamLLE(i,j,k) Slhadata(292+i+3*(j)+9*(k))
#define RVLamLLEIn_lamLLEFlat(i) Slhadata(304+i)

#define OffsetRVLamLQDIn 331
#define LengthRVLamLQDIn 27
#define BlockRVLamLQDIn(i) SlhaData(331+i)
#define RVLamLQDIn_lamLQD(i,j,k) Slhadata(319+i+3*(j)+9*(k))
#define RVLamLQDIn_lamLQDFlat(i) Slhadata(331+i)

#define OffsetRVLamUDDIn 358
#define LengthRVLamUDDIn 27
#define BlockRVLamUDDIn(i) SlhaData(358+i)
#define RVLamUDDIn_lamUDD(i,j,k) Slhadata(346+i+3*(j)+9*(k))
#define RVLamUDDIn_lamUDDFlat(i) Slhadata(358+i)

#define OffsetRVLamLLE 385
#define LengthRVLamLLE 28
#define BlockRVLamLLE(i) SlhaData(385+i)
#define RVLamLLE_Q SlhaData(386)
#define RVLamLLE_lamLLE(i,j,k) Slhadata(374+i+3*(j)+9*(k))
#define RVLamLLE_lamLLEFlat(i) Slhadata(386+i)

#define OffsetRVLamLQD 413
#define LengthRVLamLQD 28
#define BlockRVLamLQD(i) SlhaData(413+i)
#define RVLamLQD_Q SlhaData(414)
#define RVLamLQD_lamLQD(i,j,k) Slhadata(402+i+3*(j)+9*(k))
#define RVLamLQD_lamLQDFlat(i) Slhadata(414+i)

#define OffsetRVLamUDD 441
#define LengthRVLamUDD 28
#define BlockRVLamUDD(i) SlhaData(441+i)
#define RVLamUDD_Q SlhaData(442)
#define RVLamUDD_lamUDD(i,j,k) Slhadata(430+i+3*(j)+9*(k))
#define RVLamUDD_lamUDDFlat(i) Slhadata(442+i)

#define OffsetRVTLLEIn 469
#define LengthRVTLLEIn 27
#define BlockRVTLLEIn(i) SlhaData(469+i)
#define RVTLLEIn_TLLE(i,j,k) Slhadata(457+i+3*(j)+9*(k))
#define RVTLLEIn_TLLEFlat(i) Slhadata(469+i)

#define OffsetRVTLQDIn 496
#define LengthRVTLQDIn 27
#define BlockRVTLQDIn(i) SlhaData(496+i)
#define RVTLQDIn_TLQD(i,j,k) Slhadata(484+i+3*(j)+9*(k))
#define RVTLQDIn_TLQDFlat(i) Slhadata(496+i)

#define OffsetRVTUDDIn 523
#define LengthRVTUDDIn 27
#define BlockRVTUDDIn(i) SlhaData(523+i)
#define RVTUDDIn_TUDD(i,j,k) Slhadata(511+i+3*(j)+9*(k))
#define RVTUDDIn_TUDDFlat(i) Slhadata(523+i)

#define OffsetRVTLLE 550
#define LengthRVTLLE 28
#define BlockRVTLLE(i) SlhaData(550+i)
#define RVTLLE_Q SlhaData(551)
#define RVTLLE_TLLE(i,j,k) Slhadata(539+i+3*(j)+9*(k))
#define RVTLLE_TLLEFlat(i) Slhadata(551+i)

#define OffsetRVTLQD 578
#define LengthRVTLQD 28
#define BlockRVTLQD(i) SlhaData(578+i)
#define RVTLQD_Q SlhaData(579)
#define RVTLQD_TLQD(i,j,k) Slhadata(567+i+3*(j)+9*(k))
#define RVTLQD_TLQDFlat(i) Slhadata(579+i)

#define OffsetRVTUDD 606
#define LengthRVTUDD 28
#define BlockRVTUDD(i) SlhaData(606+i)
#define RVTUDD_Q SlhaData(607)
#define RVTUDD_TUDD(i,j,k) Slhadata(595+i+3*(j)+9*(k))
#define RVTUDD_TUDDFlat(i) Slhadata(607+i)

#define OffsetRVKappaIn 634
#define LengthRVKappaIn 3
#define BlockRVKappaIn(i) SlhaData(634+i)
#define RVKappaIn_kappa(i) Slhadata(634+i)

#define OffsetRVKappa 637
#define LengthRVKappa 4
#define BlockRVKappa(i) SlhaData(637+i)
#define RVKappa_Q SlhaData(638)
#define RVKappa_kappa(i) Slhadata(638+i)

#define OffsetRVDIn 641
#define LengthRVDIn 3
#define BlockRVDIn(i) SlhaData(641+i)
#define RVDIn_D(i) Slhadata(641+i)

#define OffsetRVD 644
#define LengthRVD 4
#define BlockRVD(i) SlhaData(644+i)
#define RVD_Q SlhaData(645)
#define RVD_D(i) Slhadata(645+i)

#define OffsetRVSnVEVIn 648
#define LengthRVSnVEVIn 3
#define BlockRVSnVEVIn(i) SlhaData(648+i)
#define RVSnVEVIn_VEV(i) Slhadata(648+i)

#define OffsetRVSnVEV 651
#define LengthRVSnVEV 4
#define BlockRVSnVEV(i) SlhaData(651+i)
#define RVSnVEV_Q SlhaData(652)
#define RVSnVEV_VEV(i) Slhadata(652+i)

#define OffsetRVM2LH1In 655
#define LengthRVM2LH1In 3
#define BlockRVM2LH1In(i) SlhaData(655+i)
#define RVM2LH1In_M2LH1(i) Slhadata(655+i)

#define OffsetRVM2LH1 658
#define LengthRVM2LH1 4
#define BlockRVM2LH1(i) SlhaData(658+i)
#define RVM2LH1_Q SlhaData(659)
#define RVM2LH1_M2LH1(i) Slhadata(659+i)

#define OffsetRVNMix 662
#define LengthRVNMix 49
#define BlockRVNMix(i) SlhaData(662+i)
#define RVNMix_ZNeu(n1,n2) Slhadata(655+n1+7*(n2))
#define RVNMix_ZNeuFlat(i) Slhadata(662+i)

#define OffsetRVUMix 711
#define LengthRVUMix 25
#define BlockRVUMix(i) SlhaData(711+i)
#define RVUMix_UCha(c1,c2) Slhadata(706+c1+5*(c2))
#define RVUMix_UChaFlat(i) Slhadata(711+i)

#define OffsetRVVMix 736
#define LengthRVVMix 25
#define BlockRVVMix(i) SlhaData(736+i)
#define RVVMix_VCha(c1,c2) Slhadata(731+c1+5*(c2))
#define RVVMix_VChaFlat(i) Slhadata(736+i)

#define OffsetRVHMix 761
#define LengthRVHMix 25
#define BlockRVHMix(i) SlhaData(761+i)
#define RVHMix_UH(h1,h2) Slhadata(756+h1+5*(h2))
#define RVHMix_UHFlat(i) Slhadata(761+i)

#define OffsetRVAMix 786
#define LengthRVAMix 25
#define BlockRVAMix(i) SlhaData(786+i)
#define RVAMix_UA(h1,h2) Slhadata(781+h1+5*(h2))
#define RVAMix_UAFlat(i) Slhadata(786+i)

#define OffsetRVLMix 811
#define LengthRVLMix 64
#define BlockRVLMix(i) SlhaData(811+i)
#define RVLMix_CLep(l1,l2) Slhadata(803+l1+8*(l2))
#define RVLMix_CLepFlat(i) Slhadata(811+i)

#define OffsetVCKMIn 875
#define LengthVCKMIn 4
#define BlockVCKMIn(i) SlhaData(875+i)
#define VCKMIn_lambda Slhadata(876)
#define VCKMIn_A Slhadata(877)
#define VCKMIn_rho Slhadata(878)
#define VCKMIn_eta Slhadata(879)

#define OffsetVCKM 879
#define LengthVCKM 10
#define BlockVCKM(i) SlhaData(879+i)
#define VCKM_Q SlhaData(880)
#define VCKM_VCKM(g1,g2) Slhadata(877+g1+3*(g2))
#define VCKM_VCKMFlat(i) Slhadata(880+i)

#define OffsetUPMNSIn 889
#define LengthUPMNSIn 6
#define BlockUPMNSIn(i) SlhaData(889+i)
#define UPMNSIn_theta12 Slhadata(890)
#define UPMNSIn_theta23 Slhadata(891)
#define UPMNSIn_theta13 Slhadata(892)
#define UPMNSIn_delta13 Slhadata(893)
#define UPMNSIn_alpha1 Slhadata(894)
#define UPMNSIn_alpha2 Slhadata(895)

#define OffsetUPMNS 895
#define LengthUPMNS 10
#define BlockUPMNS(i) SlhaData(895+i)
#define UPMNS_Q SlhaData(896)
#define UPMNS_UPMNS(g1,g2) Slhadata(893+g1+3*(g2))
#define UPMNS_UPMNSFlat(i) Slhadata(896+i)

#define OffsetMSS2In 905
#define LengthMSS2In 45
#define BlockMSS2In(i) SlhaData(905+i)
#define MSS2In_MSS2(g1,g2,q) Slhadata(893+g1+3*(g2)+9*(q))
#define MSS2In_MSS2Flat(i,q) Slhadata(896+i+9*(q))

#define OffsetMSL2In 905
#define LengthMSL2In 9
#define BlockMSL2In(i) SlhaData(905+i)
#define   MSL2In_MSL2(g1,g2) MSS2In_MSS2(g1,g2,1)
#define   MSL2In_MSL2Flat(i) MSS2In_MSS2Flat(i,1)

#define OffsetMSE2In 914
#define LengthMSE2In 9
#define BlockMSE2In(i) SlhaData(914+i)
#define   MSE2In_MSE2(g1,g2) MSS2In_MSS2(g1,g2,2)
#define   MSE2In_MSE2Flat(i) MSS2In_MSS2Flat(i,2)

#define OffsetMSQ2In 923
#define LengthMSQ2In 9
#define BlockMSQ2In(i) SlhaData(923+i)
#define   MSQ2In_MSQ2(g1,g2) MSS2In_MSS2(g1,g2,3)
#define   MSQ2In_MSQ2Flat(i) MSS2In_MSS2Flat(i,3)

#define OffsetMSU2In 932
#define LengthMSU2In 9
#define BlockMSU2In(i) SlhaData(932+i)
#define   MSU2In_MSU2(g1,g2) MSS2In_MSS2(g1,g2,4)
#define   MSU2In_MSU2Flat(i) MSS2In_MSS2Flat(i,4)

#define OffsetMSD2In 941
#define LengthMSD2In 9
#define BlockMSD2In(i) SlhaData(941+i)
#define   MSD2In_MSD2(g1,g2) MSS2In_MSS2(g1,g2,5)
#define   MSD2In_MSD2Flat(i) MSS2In_MSS2Flat(i,5)

#define OffsetMSS2 950
#define LengthMSS2 50
#define BlockMSS2(i) SlhaData(950+i)
#define MSS2_Q(q) Slhadata(941+10*(q))
#define MSS2_MSS2(g1,g2,q) Slhadata(938+g1+3*(g2)+10*(q))
#define MSS2_MSS2Flat(i,q) Slhadata(941+i+10*(q))

#define OffsetMSL2 950
#define LengthMSL2 10
#define BlockMSL2(i) SlhaData(950+i)
#define   MSL2_Q MSS2_Q(1)
#define   MSL2_MSL2(g1,g2) MSS2_MSS2(g1,g2,1)
#define   MSL2_MSL2Flat(i) MSS2_MSS2Flat(i,1)

#define OffsetMSE2 960
#define LengthMSE2 10
#define BlockMSE2(i) SlhaData(960+i)
#define   MSE2_Q MSS2_Q(2)
#define   MSE2_MSE2(g1,g2) MSS2_MSS2(g1,g2,2)
#define   MSE2_MSE2Flat(i) MSS2_MSS2Flat(i,2)

#define OffsetMSQ2 970
#define LengthMSQ2 10
#define BlockMSQ2(i) SlhaData(970+i)
#define   MSQ2_Q MSS2_Q(3)
#define   MSQ2_MSQ2(g1,g2) MSS2_MSS2(g1,g2,3)
#define   MSQ2_MSQ2Flat(i) MSS2_MSS2Flat(i,3)

#define OffsetMSU2 980
#define LengthMSU2 10
#define BlockMSU2(i) SlhaData(980+i)
#define   MSU2_Q MSS2_Q(4)
#define   MSU2_MSU2(g1,g2) MSS2_MSS2(g1,g2,4)
#define   MSU2_MSU2Flat(i) MSS2_MSS2Flat(i,4)

#define OffsetMSD2 990
#define LengthMSD2 10
#define BlockMSD2(i) SlhaData(990+i)
#define   MSD2_Q MSS2_Q(5)
#define   MSD2_MSD2(g1,g2) MSS2_MSS2(g1,g2,5)
#define   MSD2_MSD2Flat(i) MSS2_MSS2Flat(i,5)

#define OffsetTfIn 1000
#define LengthTfIn 27
#define BlockTfIn(i) SlhaData(1000+i)
#define TfIn_Tf(g1,g2,t) Slhadata(979+g1+3*(g2)+9*(t))
#define TfIn_TfFlat(i,t) Slhadata(982+i+9*(t))

#define OffsetTeIn 1000
#define LengthTeIn 9
#define BlockTeIn(i) SlhaData(1000+i)
#define   TeIn_Tf(g1,g2) TfIn_Tf(g1,g2,2)
#define   TeIn_TfFlat(i) TfIn_TfFlat(i,2)

#define OffsetTuIn 1009
#define LengthTuIn 9
#define BlockTuIn(i) SlhaData(1009+i)
#define   TuIn_Tf(g1,g2) TfIn_Tf(g1,g2,3)
#define   TuIn_TfFlat(i) TfIn_TfFlat(i,3)

#define OffsetTdIn 1018
#define LengthTdIn 9
#define BlockTdIn(i) SlhaData(1018+i)
#define   TdIn_Tf(g1,g2) TfIn_Tf(g1,g2,4)
#define   TdIn_TfFlat(i) TfIn_TfFlat(i,4)

#define OffsetTf 1027
#define LengthTf 30
#define BlockTf(i) SlhaData(1027+i)
#define Tf_Q(t) Slhadata(1008+10*(t))
#define Tf_Tf(g1,g2,t) Slhadata(1005+g1+3*(g2)+10*(t))
#define Tf_TfFlat(i,t) Slhadata(1008+i+10*(t))

#define OffsetTe 1027
#define LengthTe 10
#define BlockTe(i) SlhaData(1027+i)
#define   Te_Q Tf_Q(2)
#define   Te_Tf(g1,g2) Tf_Tf(g1,g2,2)
#define   Te_TfFlat(i) Tf_TfFlat(i,2)

#define OffsetTu 1037
#define LengthTu 10
#define BlockTu(i) SlhaData(1037+i)
#define   Tu_Q Tf_Q(3)
#define   Tu_Tf(g1,g2) Tf_Tf(g1,g2,3)
#define   Tu_TfFlat(i) Tf_TfFlat(i,3)

#define OffsetTd 1047
#define LengthTd 10
#define BlockTd(i) SlhaData(1047+i)
#define   Td_Q Tf_Q(4)
#define   Td_Tf(g1,g2) Tf_Tf(g1,g2,4)
#define   Td_TfFlat(i) Tf_TfFlat(i,4)

#define OffsetASfMix 1057
#define LengthASfMix 144
#define BlockASfMix(i) SlhaData(1057+i)
#define ASfMix_UASf(s1,s2,t) Slhadata(1015+s1+6*(s2)+36*(t))
#define ASfMix_UASfFlat(i,t) Slhadata(1021+i+36*(t))

#define OffsetSnuMix 1057
#define LengthSnuMix 36
#define BlockSnuMix(i) SlhaData(1057+i)
#define   SnuMix_UASf(s1,s2) ASfMix_UASf(s1,s2,1)
#define   SnuMix_UASfFlat(i) ASfMix_UASfFlat(i,1)

#define OffsetSelMix 1093
#define LengthSelMix 36
#define BlockSelMix(i) SlhaData(1093+i)
#define   SelMix_UASf(s1,s2) ASfMix_UASf(s1,s2,2)
#define   SelMix_UASfFlat(i) ASfMix_UASfFlat(i,2)

#define OffsetUSqMix 1129
#define LengthUSqMix 36
#define BlockUSqMix(i) SlhaData(1129+i)
#define   USqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,3)
#define   USqMix_UASfFlat(i) ASfMix_UASfFlat(i,3)

#define OffsetDSqMix 1165
#define LengthDSqMix 36
#define BlockDSqMix(i) SlhaData(1165+i)
#define   DSqMix_UASf(s1,s2) ASfMix_UASf(s1,s2,4)
#define   DSqMix_UASfFlat(i) ASfMix_UASfFlat(i,4)

#define OffsetSnsMix 1201
#define LengthSnsMix 9
#define BlockSnsMix(i) SlhaData(1201+i)
#define SnsMix_US(g1,g2) Slhadata(1198+g1+3*(g2))
#define SnsMix_USFlat(i) Slhadata(1201+i)

#define OffsetSnaMix 1210
#define LengthSnaMix 9
#define BlockSnaMix(i) SlhaData(1210+i)
#define SnaMix_UA(g1,g2) Slhadata(1207+g1+3*(g2))
#define SnaMix_UAFlat(i) Slhadata(1210+i)

#define OffsetCVHMix 1219
#define LengthCVHMix 16
#define BlockCVHMix(i) SlhaData(1219+i)
#define CVHMix_UH(h1,h2) Slhadata(1215+h1+4*(h2))
#define CVHMix_UHFlat(i) Slhadata(1219+i)

#define OffsetNMNMix 1235
#define LengthNMNMix 25
#define BlockNMNMix(i) SlhaData(1235+i)
#define NMNMix_ZNeu(n1,n2) Slhadata(1230+n1+5*(n2))
#define NMNMix_ZNeuFlat(i) Slhadata(1235+i)

#define OffsetNMHMix 1260
#define LengthNMHMix 9
#define BlockNMHMix(i) SlhaData(1260+i)
#define NMHMix_UH(h1,h2) Slhadata(1257+h1+3*(h2))
#define NMHMix_UHFlat(i) Slhadata(1260+i)

#define OffsetNMAMix 1269
#define LengthNMAMix 9
#define BlockNMAMix(i) SlhaData(1269+i)
#define NMAMix_UA(h1,h2) Slhadata(1266+h1+3*(h2))
#define NMAMix_UAFlat(i) Slhadata(1269+i)

#define OffsetPrecObs 1278
#define LengthPrecObs 11
#define BlockPrecObs(i) SlhaData(1278+i)
#define PrecObs_DeltaRho Slhadata(1279)
#define PrecObs_MWMSSM Slhadata(1280)
#define PrecObs_MWSM Slhadata(1281)
#define PrecObs_SW2effMSSM Slhadata(1282)
#define PrecObs_SW2effSM Slhadata(1283)
#define PrecObs_gminus2mu Slhadata(1284)
#define PrecObs_EDMeTh Slhadata(1285)
#define PrecObs_EDMn Slhadata(1286)
#define PrecObs_EDMHg Slhadata(1287)
#define PrecObs_bsgammaMSSM Slhadata(1288)
#define PrecObs_bsgammaSM Slhadata(1289)

#define OffsetSPInfo 1289
#define LengthSPInfo 92
#define BlockSPInfo(i) SlhaData(1289+i)
#define SPInfo_Severity SlhaData(1290)
#define SPInfo_NLines SlhaData(1291)
#define SPInfo_Code(n) SlhaData(1291+n)
#define SPInfo_Text(i,n) SlhaData(1301+i+5*(n))
#define SPInfo_TextFlat(i) SlhaData(1306+i)
#define   SPInfo_Len 80

#define OffsetDCInfo 1381
#define LengthDCInfo 92
#define BlockDCInfo(i) SlhaData(1381+i)
#define DCInfo_Severity SlhaData(1382)
#define DCInfo_NLines SlhaData(1383)
#define DCInfo_Code(n) SlhaData(1383+n)
#define DCInfo_Text(i,n) SlhaData(1393+i+5*(n))
#define DCInfo_TextFlat(i) SlhaData(1398+i)
#define   DCInfo_Len 80

#define OffsetDecays 1473
#define LengthDecays 4096
#define BlockDecays(i) SlhaData(1473+i)
#define Decays_Data(n) Slhadata(1473+n)

#define nslhadata 5569

#endif

!
      SUBROUTINE GOUT601(G,QEND,PRCH,OUTFILE)
!
!Purpose: Writes all the results in a readable format
!         PRCH is a switch to alter the output.
!             =1 -> output with no tilde-couplings or lambdas,
!                   most appropriate for high scale output
!             =2 -> output with tilde-couplings and lambdas,
!                   no regular Yukawas,
!                   most appropriate for low scale output
!             =3 -> everything,
!                   only appropriate at m_H
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,GUTOUT,WKOUT
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      DOUBLE PRECISION QEND
      DOUBLE COMPLEX G(601)
      CHARACTER*20 OUTFILE
      INTEGER I,SW,PRCH
!
!Now write out the results
!
      OPEN(15,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(15,*)
      WRITE(15,*)'Q IS',QEND,' COMP IS',COMP
      WRITE(15,*)
!
      WRITE(15,*)'COUPLINGS:'
      WRITE(15,10)G(1),G(2),G(3)
      WRITE(15,*)
!
      IF(PRCH.NE.2)THEN
        WRITE(15,*)'f_U:'
        DO I=1,3
          WRITE(15,10)G(4+3*(I-1)),G(5+3*(I-1)),G(6+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'f_D:'
        DO I=1,3
          WRITE(15,10)G(13+3*(I-1)),G(14+3*(I-1)),G(15+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'f_E:'
        DO I=1,3
          WRITE(15,10)G(22+3*(I-1)),G(23+3*(I-1)),G(24+3*(I-1))
        END DO
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'GAUGINO MASSES M:'
      WRITE(15,10)G(31),G(32),G(33)
      WRITE(15,*)
!
      IF(COMP.EQ.1)THEN
        WRITE(15,*)'GAUGINO MASSES M'':'
        WRITE(15,10)G(599),G(600),G(601)
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'a_U:'
      DO I=1,3
        WRITE(15,10)G(34+3*(I-1)),G(35+3*(I-1)),G(36+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'a_D:'
      DO I=1,3
        WRITE(15,10)G(43+3*(I-1)),G(44+3*(I-1)),G(45+3*(I-1)) 
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'a_E:'
      DO I=1,3
        WRITE(15,10)G(52+3*(I-1)),G(53+3*(I-1)),G(54+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'SQUARED HIGGS MASSES:'
      WRITE(15,5)G(61)-G(108)**2,G(62)-G(108)**2
      WRITE(15,*)
!
      WRITE(15,*)'M_Q^2:'
      DO I=1,3
        WRITE(15,10)G(63+3*(I-1)),G(64+3*(I-1)),G(65+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_L^2:'
      DO I=1,3
        WRITE(15,10)G(72+3*(I-1)),G(73+3*(I-1)),G(74+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_U^2:'
      DO I=1,3
        WRITE(15,10)G(81+3*(I-1)),G(82+3*(I-1)),G(83+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_D^2:'
      DO I=1,3
        WRITE(15,10)G(90+3*(I-1)),G(91+3*(I-1)),G(92+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_E^2:'
      DO I=1,3
        WRITE(15,10)G(99+3*(I-1)),G(100+3*(I-1)),G(101+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MU AND B:'
      WRITE(15,5)G(108),G(109)
      WRITE(15,*)
      WRITE(15,*)'V_U AND V_D:'
      WRITE(15,5)G(110),G(111)
      WRITE(15,*)
!
      IF(PRCH.NE.1)THEN
        WRITE(15,*)'lambda_U:'
        DO I=1,3
          WRITE(15,10)G(112+3*(I-1)),G(113+3*(I-1)),G(114+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'lambda_D:'
        DO I=1,3
          WRITE(15,10)G(121+3*(I-1)),G(122+3*(I-1)),G(123+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'lambda_E:'
        DO I=1,3
          WRITE(15,10)G(130+3*(I-1)),G(131+3*(I-1)),G(132+3*(I-1))
        END DO
        WRITE(15,*)
  
        WRITE(15,*)'GTP_Q:'
        DO I=1,3
          WRITE(15,10)G(139+3*(I-1)),G(140+3*(I-1)),G(141+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_L:'
        DO I=1,3
          WRITE(15,10)G(148+3*(I-1)),G(149+3*(I-1)),G(150+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_U:'
        DO I=1,3
          WRITE(15,10)G(157+3*(I-1)),G(158+3*(I-1)),G(159+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_D:'
        DO I=1,3
          WRITE(15,10)G(166+3*(I-1)),G(167+3*(I-1)),G(168+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_E:'
        DO I=1,3
          WRITE(15,10)G(175+3*(I-1)),G(176+3*(I-1)),G(177+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTPH_U AND GTPH_D:'
        WRITE(15,5)G(184),G(185)
        WRITE(15,*)
!
        WRITE(15,*)'GT_Q:'
        DO I=1,3
          WRITE(15,10)G(186+3*(I-1)),G(187+3*(I-1)),G(188+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GT_L:'
        DO I=1,3
          WRITE(15,10)G(195+3*(I-1)),G(196+3*(I-1)),G(197+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTH_U AND GTH_D:'
        WRITE(15,5)G(204),G(205)
        WRITE(15,*)
!
        WRITE(15,*)'GTS_Q:'
        DO I=1,3
          WRITE(15,10)G(206+3*(I-1)),G(207+3*(I-1)),G(208+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTS_U:'
        DO I=1,3
          WRITE(15,10)G(215+3*(I-1)),G(216+3*(I-1)),G(217+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTS_D:'
        DO I=1,3
          WRITE(15,10)G(224+3*(I-1)),G(225+3*(I-1)),G(226+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTQ_U:'
        DO I=1,3
          WRITE(15,10)G(233+3*(I-1)),G(234+3*(I-1)),G(235+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTQ_D:'
        DO I=1,3
          WRITE(15,10)G(242+3*(I-1)),G(243+3*(I-1)),G(244+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTL_E:'
        DO I=1,3
          WRITE(15,10)G(251+3*(I-1)),G(252+3*(I-1)),G(253+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTU_U:'
        DO I=1,3
          WRITE(15,10)G(260+3*(I-1)),G(261+3*(I-1)),G(262+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTD_D:'
        DO I=1,3
          WRITE(15,10)G(269+3*(I-1)),G(270+3*(I-1)),G(271+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTE_E:'
        DO I=1,3
          WRITE(15,10)G(278+3*(I-1)),G(279+3*(I-1)),G(280+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'sGTPH_U AND cGTPH_D:'
        WRITE(15,5)G(287),G(288)
        WRITE(15,*)
!
        WRITE(15,*)'sGTH_U AND cGTH_D:'
        WRITE(15,5)G(289),G(290)
        WRITE(15,*)
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'MSSM COUPLINGS:'
      WRITE(15,10)G(291),G(292),G(293)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_U:'
      DO I=1,3
        WRITE(15,10)G(294+3*(I-1)),G(295+3*(I-1)),G(296+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_D:'
      DO I=1,3
        WRITE(15,10)G(303+3*(I-1)),G(304+3*(I-1)),G(305+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_E:'
      DO I=1,3
        WRITE(15,10)G(312+3*(I-1)),G(313+3*(I-1)),G(314+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM GAUGINO MASSES:'
      WRITE(15,10)G(321),G(322),G(323)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_U:'
      DO I=1,3
        WRITE(15,10)G(324+3*(I-1)),G(325+3*(I-1)),G(326+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_D:'
      DO I=1,3
        WRITE(15,10)G(333+3*(I-1)),G(334+3*(I-1)),G(335+3*(I-1)) 
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_E:'
      DO I=1,3
        WRITE(15,10)G(342+3*(I-1)),G(343+3*(I-1)),G(344+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM SQUARED HIGGS MASSES:'
      WRITE(15,5)G(351),G(352)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_Q^2:'
      DO I=1,3
        WRITE(15,10)G(353+3*(I-1)),G(354+3*(I-1)),G(355+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_L^2:'
      DO I=1,3
        WRITE(15,10)G(362+3*(I-1)),G(363+3*(I-1)),G(364+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_U^2:'
      DO I=1,3
        WRITE(15,10)G(371+3*(I-1)),G(372+3*(I-1)),G(373+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_D^2:'
      DO I=1,3
        WRITE(15,10)G(380+3*(I-1)),G(381+3*(I-1)),G(382+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_E^2:'
      DO I=1,3
        WRITE(15,10)G(389+3*(I-1)),G(390+3*(I-1)),G(391+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM MU AND B:'
      WRITE(15,5)G(398),G(399)
      WRITE(15,*)
!
      IF(PRCH.NE.1)THEN
        WRITE(15,*)'TRI_U:'
        DO I=1,3
          WRITE(15,10)G(400+3*(I-1)),G(401+3*(I-1)),G(402+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'TRI_D:'
        DO I=1,3
          WRITE(15,10)G(409+3*(I-1)),G(410+3*(I-1)),G(411+3*(I-1)) 
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'TRI_E:'
        DO I=1,3
          WRITE(15,10)G(418+3*(I-1)),G(419+3*(I-1)),G(420+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'M_HUD:'
        WRITE(15,7)G(427)
        WRITE(15,*)
!
        WRITE(15,*)'SM VEV AND LAMBDA_SM:'
        WRITE(15,5)G(428),G(429)
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_U:'
        DO I=1,3
          WRITE(15,10)G(430+3*(I-1)),G(431+3*(I-1)),G(432+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_D:'
        DO I=1,3
          WRITE(15,10)G(439+3*(I-1)),G(440+3*(I-1)),G(441+3*(I-1)) 
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_E:'
        DO I=1,3
          WRITE(15,10)G(448+3*(I-1)),G(449+3*(I-1)),G(450+3*(I-1))
        END DO
        WRITE(15,*)
      END IF
!
      CLOSE(15)
!
 5    FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')',2X,'(',D11.4,',',D11.4,')'
     $       ) !Used to be 11.4
 7    FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')')
 10   FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')',2X,'(',D11.4,',',D11.4,')'
     $      ,2X,'(',D11.4,',',D11.4,')')
!
      RETURN
      END

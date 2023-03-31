*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

**************************************************************
*** include file                                           ***
*** this file contains common blocks for dsas files        ***
*** AUTHOR: Piero Ullio, piero@physto.se                   ***
*** Date: 99-06-22                                         ***
**************************************************************

c family indices (set in dsinit)
      integer ivfam(0:50),ivtype(0:50)
      common /asfamily/ ivfam,ivtype
c internal family indices
      integer iifam(4),itype(4)
      common /asfamilyint/ iifam,itype
c mass parameters
      real*8 mass1,mass2,mass3,mass4
      common /asparmass/ mass1,mass2,mass3,mass4
c symmetry factor of the final state
      integer s34
c and type of fermions in the final state
      integer  q3,q4
      common /aspartype/ q3,q4,s34
c kinematical input values
      real*8 p12,costheta
      common /askin/ p12,costheta
c kinematical values defined in first step
      real*8 Svar,Tvar,Uvar,k34,ep1,ep2,ek3,ek4,scal(4,4)
      common /askinder/ Svar,Tvar,Uvar,k34,ep1,ep2,ek3,ek4,
     & scal
c terms to compute a fermion line
      complex*16 ASxpl(4),ASxpr(4),ASyl,ASyr
      common /asampli/ ASxpl,ASxpr,ASyl,ASyr
c terms to compute a fermion line
      complex*16 ASxplc(6,4),ASxprc(6,4),ASylc(6),ASyrc(6),
     & colfactor(6,6)
      common /asamplic/ ASxplc,ASxprc,ASylc,ASyrc,colfactor
c some roundoff terms
      real*8 thstep,fertoll
      common /astoll/ thstep,fertoll 

c dsasdwdcossfsf final/intermediate state arrays/variables
      real*8  cfactini,cfactfin,gg1,gg2
      integer kp1s,kp2s,kp3in(30),kp4in(30),chcol,nsfertc,ksfertc(6),
     &  nsferuc,ksferuc(6),nsfertn,ksfertn(6),nsferun,ksferun(6),
     &  kf2,kf2o,ick1,ick2
      logical chon(30),gluonin,gammain,neutcurr,nosneutrinov
      common/aswsfsfcom/cfactini,cfactfin,gg1,gg2,kp1s,kp2s,kp3in,kp4in,
     &  chcol,nsfertc,ksfertc,nsferuc,ksferuc,nsfertn,ksfertn,nsferun,
     &  ksferun,kf2,kf2o,ick1,ick2,chon,gluonin,gammain,neutcurr,
     &  nosneutrinov

c dsasdwdcossfchi final/intermediate state arrays/variables
      real*8  gg1c,gg2c
      integer kp1c,kp2c,ciaux,kcfers,ncfers,kcfersv(3),ncferd,kcferd(3),
     & ncsfert,kcsfertn(2),ncsfertc,kcsfertc(6)
      logical cgammain,cgluonin
      common/aswsfchicom/gg1c,gg2c,kp1c,kp2c,ciaux,kcfers,ncfers,
     & kcfersv,ncferd,kcferd,ncsfert,kcsfertn,ncsfertc,kcsfertc,
     & cgammain,cgluonin

c logical value to set whether you want to have a warning statement in 
c case of negative partial results
      logical aszeroprint
      common/aswprintcom/aszeroprint

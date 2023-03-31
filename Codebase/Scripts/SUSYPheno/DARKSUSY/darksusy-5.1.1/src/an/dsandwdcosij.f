      function dsandwdcosij(p,costheta,kp1,kp2)
      implicit none
      include 'dsmssm.h'
      real*8 dsandwdcosij,p,costheta,ans
      real*8 dsandwdcosnn,dsandwdcoscn,dsandwdcoscc,
     &     dsasdwdcossfsf,dsasdwdcossfchi
      logical dsisneutralino,dsischargino,dsissfermion
      integer kp1,kp2
      ans=0.d0
      if (dsisneutralino(kp1).and.dsisneutralino(kp2)) then
         ans=dsandwdcosnn(p,costheta,kp1,kp2)
      else if (dsisneutralino(kp1).and.dsischargino(kp2)) then
         ans=dsandwdcoscn(p,costheta,kp2,kp1)
      else if (dsischargino(kp1).and.dsisneutralino(kp2)) then
         ans=dsandwdcoscn(p,costheta,kp1,kp2)
      else if (dsischargino(kp1).and.dsischargino(kp2)) then
         ans=dsandwdcoscc(p,costheta,kp1,kp2)
      else if (dsissfermion(kp1).and.dsissfermion(kp2)) then
         ans=dsasdwdcossfsf(p,costheta,kp1,kp2)
      else if ((dsissfermion(kp1).and.dsisneutralino(kp2)).or.
     &         (dsissfermion(kp1).and.dsischargino(kp2))) then
         ans=dsasdwdcossfchi(p,costheta,kp1,kp2)
      else if ((dsisneutralino(kp1).and.dsissfermion(kp2)).or.
     &         (dsischargino(kp1).and.dsissfermion(kp2))) then
         ans=dsasdwdcossfchi(p,costheta,kp2,kp1)
      else
         write(*,*) 'DS: ERROR in dsandwdcosij: unimplemented initial state ',
     &  kp1,'(',pname(kp1),'), ',kp2,'(',pname(kp2),')'
         stop
      endif
      dsandwdcosij=ans
      return
      end

      function dsisneutralino(kp)
      implicit none
      include 'dsmssm.h'
      logical dsisneutralino,ans
      integer kp,i
      ans=.false.
      if (kp.ge.kn(1).and.kp.le.kn(4)) ans=.true.
      dsisneutralino=ans
      return
      end

      function dsischargino(kp)
      implicit none
      include 'dsmssm.h'
      logical dsischargino,ans
      integer kp,i
      dsischargino=kp.eq.kcha1.or.kp.eq.kcha2
      return
      end

      function dsissfermion(kp)
      implicit none
      include 'dsmssm.h'
      logical dsissfermion
      integer kp,i
      dsissfermion=kp.eq.ksnue.or.kp.eq.kse1.or.kp.eq.kse2.or.kp.eq.ksnumu
     &    .or.kp.eq.ksmu1.or.kp.eq.ksmu2.or.kp.eq.ksnutau.or.kp.eq.kstau1
     &    .or.kp.eq.kstau2.or.kp.eq.ksu1.or.kp.eq.ksu2.or.kp.eq.ksd1.or.
     &    kp.eq.ksd2.or.kp.eq.ksc1.or.kp.eq.ksc2.or.kp.eq.kss1.or.
     &    kp.eq.kss2.or.kp.eq.kst1.or.kp.eq.kst2.or.kp.eq.ksb1.or.kp.eq.ksb2
      return
      end


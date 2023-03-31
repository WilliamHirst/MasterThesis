      function dsgm2muon()
      implicit none
c supersymmetric contribution to (g-2)_muon
c  output:
c    gm2amp : susy contribution to g-2 amplitude  = (g-2)/2
c according to T Moroi hep-ph/9512396 v3
c (T. Moroi, PRD 1996; (E) 1997)

c author: paolo gondolo 2001-02-08
c reference: e a baltz and p gondolo, hep-ph/0102147
      include 'dsmssm.h'
      include 'dsge.h'
      real*8 dsgm2muon,daneu,dacha,tmp,x,dsabsq
      integer i,j
      daneu=0.d0
      do i=1,2
        do j=1,4
          x=(mass(kn(j))/mass(ksmu(i)))**2
          tmp=-(mass(kmu)/mass(ksmu(i)))**2/(6.d0*(1.d0-x)**4) *
     &       (dsabsq(gl(kmu,kn(j),ksmu(i)))+
     &        dsabsq(gr(kmu,kn(j),ksmu(i)))) *
     &       (1.d0-6.d0*x+3.d0*x**2+2.d0*x**3-6.d0*x**2*log(x))
     &       -mass(kmu)*mass(kn(j))/mass(ksmu(i))**2/(1.d0-x)**3 *
     &       dble(gl(kmu,kn(j),ksmu(i))*conjg(gr(kmu,kn(j),ksmu(i)))) *
     &       (1.d0-x**2+2.d0*x*log(x))
          daneu=daneu+tmp
        enddo
      enddo
      dacha=0.d0
      do i=1,2
        x=(mass(kcha(i))/mass(ksnumu))**2
        tmp=(mass(kmu)/mass(ksnumu))**2/(3.d0*(1.d0-x)**4) *
     &     (dsabsq(gl(ksnumu,kmu,kcha(i)))+
     &      dsabsq(gr(ksnumu,kmu,kcha(i)))) *
     &     (1.d0+1.5d0*x-3.d0*x**2+0.5d0*x**3+3.d0*x*log(x))
     &     -mass(kmu)*mass(kcha(i))/mass(ksnumu)**2/(1.d0-x)**3 *
     &     dble(gl(ksnumu,kmu,kcha(i))*conjg(gr(ksnumu,kmu,kcha(i)))) *
     &     (3.d0-4.d0*x+x**2+2.d0*log(x))
        dacha=dacha+tmp
      enddo
      dsgm2muon=(daneu+dacha)/(16.d0*pi**2)
      return
      end

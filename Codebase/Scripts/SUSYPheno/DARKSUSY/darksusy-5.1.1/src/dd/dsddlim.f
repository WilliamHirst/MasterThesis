      function dsddlim(mx,iexp)
c
c     limits on scattering cross section on nucleons from
c     direct dark matter searches
c
c     input: mx - wimp mass
c            iexp - experiment
c                   1 = future cawo_4 cresst
c                   2 = future genius
c                   3 = dama 1997, plb389, 757
c     output: dsddlim - upper limit or sensitivity limit on wimp-nucleon
c                     spin-independent cross section in pb
c                     (returns 10^99 if no limit)
c
c     author: paolo gondolo 1999
c
      implicit none
      include 'dsmpconst.h'
      real*8 dsddlim,mx
      integer iexp
      integer opt ! opt=1 signal < sqrt(background), opt=2 signal < background
      integer i,n,nmax ! max number of nuclear species
      parameter (nmax=5)
      real*8 expos,eresol,eth,rho,uo,urms,mmol,fn,en,qn,rn,sn,
     &     mnucx,vn,mnuc,mpx,xo,xn,z,dsffn,eta,cn,cc,drdedsigma,b,k,erf
c      external erf
      real*8 an(nmax),stech(nmax),quench(nmax),alphan(nmax)

      if (iexp.eq.1) then ! future cresst with calcium tungstate cawo_4
         opt=1
         expos=300.d0           ! exposure in kg yr
         b=0.1d0*1.d-4 ! background in counts/kg/kev/day after rejection
         eth=15.d0 ! energy threshold in kev
         eresol=0.23d0 ! energy resolution in kev
         n=3 ! number of nuclear species
         an(1)=40.078d0 ! atomic weight ------- ca
         stech(1)=1.d0 ! stoichiometric coefficient
         quench(1)=1.d0 ! quenching factor
         alphan(1)=1.d0 ! exponent in evis = quenching erecoil^alphan
         an(2)=183.84d0 ! atomic weight ------- w
         stech(2)=1.d0 ! stoichiometric coefficient
         quench(2)=1.d0 ! quenching factor
         alphan(2)=1.d0 ! exponent in evis = quenching erecoil^alphan
         an(3)=16.d0 ! atomic weight ------- o
         stech(3)=4.d0 ! stoichiometric coefficient
         quench(3)=1.d0 ! quenching factor
         alphan(3)=1.d0 ! exponent in evis = quenching erecoil^alphan
      else if (iexp.eq.2) then ! future genius
         opt=2
         expos=0.d0 ! irrelevant for opt=2
         b=0.01d0/365.25d0 ! background in counts/kg/kev/day
         eth=11.d0 ! energy threshold in kev
         eresol=0.d0 ! irrelevant for opt=2
         n=1 ! number of nuclear species
         an(1)=76.d0 ! atomic weight ------- ge76
         stech(1)=1.d0 ! stoichiometric coefficient
         quench(1)=0.14d0 ! quenching factor
         alphan(1)=1.19d0 ! exponent in evis = quenching erecoil^alphan
      else if (iexp.eq.3) then ! dama 1997
         opt=2
         expos=0.d0 ! irrelevant for opt=2
         b=0.286d0 ! background in counts/kg/kev/day
         eth=2.d0 ! energy threshold in kev
         eresol=0.d0 ! irrelevant for opt=2
         n=2 ! number of nuclear species
         an(1)=23.d0 ! atomic weight ------- na
         stech(1)=1.d0 ! stoichiometric coefficient
         quench(1)=0.30d0 ! quenching factor
         alphan(1)=1.0d0 ! exponent in evis = quenching erecoil^alphan
         an(2)=127.d0 ! atomic weight ------- i
         stech(2)=1.d0 ! stoichiometric coefficient
         quench(2)=0.09d0 ! quenching factor
         alphan(2)=1.0d0 ! exponent in evis = quenching erecoil^alphan
      else
         write (*,*) 'dsddlim: invalid code for experiment'
         stop
      endif

      rho = 0.3d0 ! local halo density in gev/cm^3
      uo = 245.d0 ! earth speed through dark halo
      urms = 270.d0 ! local halo 1d velocity dispersion in km/s

      mmol = 0.d0
      do i=1,n
         mmol=mmol+stech(i)*an(i)
      enddo
      cc = 0.d0
      do i=1,n
         fn=stech(i)*an(i)/mmol
         en=(eth/quench(i))**(1.d0/alphan(i))
         mnuc=an(i)*0.9315d0
         qn=sqrt(2.d0*mnuc*en)/197.327d0
         rn=an(i)**(1.d0/3.d0)
         sn=1.0d0
         mnucx=mnuc*mx/(mnuc+mx)
         vn=sqrt(mnuc*en*1.d-6/2.d0/mnucx**2)*3.d5
         mpx=m_p*mx/(m_p+mx)
         xo=uo/sqrt(2.d0)/urms
         xn=vn/sqrt(2.d0)/urms
         z=qn*rn
         dsffn = (3.d0*(sin(z)/z**2-cos(z)/z)/z)**2*exp(-(qn*sn)**2)
         eta = (erf(xn+xo)-erf(xn-xo))/xo
         cn = (en/eth/alphan(i))*fn*an(i)**2*dsffn*eta
         cc = cc+cn
      enddo
      drdedsigma = 4.36d5*cc*rho/(mx*mpx**2*4.d0*sqrt(2.d0)*urms)
                     ! counts/kg/kev/day/pb
      if (opt.eq.1) then
         k = sqrt(b/(expos*365.25d0*eresol)) !  counts/kg/kev/day
      else
         k = b
      endif
      dsddlim = 1.d99
      if (drdedsigma.le.0.d0) return
      dsddlim = k/drdedsigma
      return
      end


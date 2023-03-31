      real*8 function dsrdthav(x,wrate)
c_______________________________________________________________________
c  the thermal average of the effective annihilation cross section.
c  input:
c    x - mass/temperature (real)
c    wrate - invariant annihilation rate (real)
c  output:
c    dsrdthav - thermal averged cross section
c  common:
c    'dsrdcom.h' - included common blocks
c  uses qrkck or dgadap.
c  called by dsrdrhs
c  author: joakim edsjo (edsjo@physto.se) 98-05-01
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wavs,dsrdfuncs
      real*8 x,wrate,dsrdfunc,dsai,dsbi,
     &  epsin,epsout,gpindp
      integer i,iop
      external dsrdfunc,wrate,dsrdfuncs,gpindp
      real*8 wav,xmin,wtmp,utop

c-----------------------------------------------------------------------
c..added by je 97-01-17 for gadap integration
      real*8 rdx
      common /gadint2/ rdx

c-------------- thermally averaged cross section times relative velocity

c     integrate as far as umax in dsrdtab allows  ! je 97-01-07
      utop=umax*0.9999d0    ! must be <= umax, or rather the maximal u
                            ! corresponding to pmax set in dsrdtab

c... In dsrdtab, the maximum tabulated p is given, calculate the maximal
c... u that this pmax corresponds to
c      write(*,*) 'dsrdthav: pmax=',pmax
c      utop=sqrt((sqrt(mco(1)**2+pmax**2)-2*mco(1))*xmin/mco(1))
c      utop=0.9999d0*utop
c      write(*,*) '  utop=',utop,'  umax=',umax

c...define xmin. xmin >= xstart in dsrdens required.
      xmin=max(xinit,1.0001d0*mco(1)/tgev(1))  ! je corr 97-03-14
      if (xmin.ne.xinit) rdwar=ibset(rdwar,3)

      nlo=1
      nhi=2

      if (thavint.eq.1) then ! gadap integration (thavint=1)
        rdx=x
        wav=0.0d0
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &        mco(1)**2)-1.0d0))
          if (i.gt.1.and.dsai.gt.dsbi) then
             wavs=0.5d0*(dsrdfuncs(dsai)+dsrdfuncs(dsbi))*(dsai-dsbi)
             wavs=wavs/mco(1)**2
             wav=wav+wavs/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)
          call dgadap(dsai,dsbi,dsrdfuncs,0.0001d0,wavs)
c...divide out mco(1)^2 introduced in dsrdfunc
          wavs=wavs/mco(1)**2
          wav=wav+wavs/1.0d15
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      elseif (thavint.eq.2) then ! runge-kutta integration (thavint=2)

        wav=0.0d0
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &        mco(1)**2)-1.0d0))
          if (i.gt.1.and.dsai.gt.dsbi) then
             wtmp=0.5d0*(dsrdfuncs(dsai)+dsrdfuncs(dsbi))*(dsai-dsbi)
             wtmp=wtmp/mco(1)**2
             wav=wav+wtmp/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)
          call dsrdqrkck(dsrdfunc,x,wrate,dsai,dsbi,wtmp)
c...divide out mco(1)^2 introduced in dsrdfunc
          wtmp=wtmp/mco(1)**2
          wav=wav+wtmp
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      else  ! gpindp (thavint=3)

        rdx=x
        wav=0.0d0
        epsin=0.001d0
        epsout=0.001d0
        iop=1
        do i=1,nlim
          dsai=sqrt(2.0d0*x*(sqrt(1.0d0+(plow(i))**2/
     &        mco(1)**2)-1.0d0))
          if (i.gt.1.and.dsai.gt.dsbi) then
             wavs=0.5d0*(dsrdfuncs(dsai)+dsrdfuncs(dsbi))*(dsai-dsbi)
             wavs=wavs/mco(1)**2
             wav=wav+wavs/1.0d15
          endif
          dsbi=min(sqrt(2.0d0*x*(sqrt(1.0d0+(phigh(i))**2
     &      /mco(1)**2)-1.0d0)),utop)
          wavs = gpindp(dsai,dsbi,epsin,epsout,dsrdfuncs,iop)
c...divide out mco(1)^2 introduced in dsrdfunc
          wavs=wavs/mco(1)**2

          wav=wav+wavs/1.0d15
          if (iop.eq.0) then
            rderr=ibset(rderr,8)
            write(*,*) 'error in dsrdrhs: gpindp integration failed'
            write(*,*) '  for model ',rdtag
            goto 130
          endif
          if (dsbi.ge.0.9999d0*utop) goto 130
        enddo

      endif

 130  dsrdthav=wav

      return
      end









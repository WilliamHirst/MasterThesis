       subroutine dsfeynhiggsfast(hwar,HM,mh,hc,halpha,drho,
     &  ATop,ABot,inmt,inmb,my,M2,mqtl,mqtr,mqbl,mqbr,
     &  mgluino,mh3,tanb)
c
c -----------------------------------------------------------------       
c Implementation of FeynHiggsFast in DarkSUSY by J. Edsjo 2000-04-25
c
c All names of subroutines and functions that clashed with the
c full FeynHiggs names have _fast appended to them. In most cases,
c they are probably the same routines though, but to be on the
c safe side, the names were changed.
c Output: mh is lighter scalar Higgs mass, HM is heavier Higgs mass 
c -----------------------------------------------------------------
c
c     FeynHiggsFast
c     =============
c      
c       Calculation of the masses of the neutral CP-even
c       Higgs bosons in the MSSM
c       
c       Author: Sven Heinemeyer
c
c       Based on hep-ph/9903404
c       by S. Heinemeyer, W. Hollik, G. Weiglein
c      
c       In case of problems or questions,
c       contact Sven Heinemeyer :-)
c       email: Sven.Heinemeyer@desy.de
c       
c       FeynHiggs homepage:
c       http://www-itp.physik.uni-karlsruhe.de/feynhiggs/
c
c --------------------------------------------------------------
c
c Warnings implemented by J. Edsjo, 2000-04-25
c
c                   Bit  Value
c  Bits of hwar:    0 -  1:   Potential numerical problems at 1-loop
c                   1 -  2:   Potential numerical problems at 2-loop
c                   2 -  4:   Error with not used H2 mass expression 1-loop
c                   3 -  8:   Error with not used H2 mass expression 2-loop
c                   4 -  16:  Error with not used H2 mass expression 2-loop
c                   5 -  32:  1-loop Higgs sector not OK
c                   6 -  64:  2-loop Higgs sector not OK
c                   7 -  128: Stop or sbottom masses not OK
c---------------------------------------------------------------

      implicit real*8(a-z)

      include 'dsidtag.h'

c---DarkSusy declarations - begin
      real*8 halpha,hc,drho,
     &  ATop,ABot,inmt,inmb,my,M2,mqtl,mqtr,mqbl,mqbr,mgluino,
     &  mh3,tanb
      real*8 alpha1,alpha2,alpha3,v
      integer hwar,k
      logical first
      data first/.true./
      save first

c---DarkSusy declarations - end

      double precision mlhlle1, mlhlle2, 
     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $       alphadiag1, alphadiag2, mhptree, mhp
      integer selec2,selec3
      double precision delrholimit, delrholimitnew, delrhores,
     $                 delrho1loop, delrho2loopgluon,
     $                 delrho2loopgluino
      double precision msusytrnew
      double precision mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                 el, gf
      double precision cf, eps, pi
      complex*16 i
      double precision  MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma

      common /smpara/ mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                el, gf
      common /para/ cf, eps, pi, i
      common /susypara_fast/ MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma
      common /selec/ selec2, selec3
      common/higgsmass/mhh, mlh
      common /lle_fast/ mlhlle1, mlhlle2,
     $             mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $             alphadiag1, alphadiag2, mhptree, mhp

c---DarkSusy addition begin
      integer jehwar,prlev
      character*12 jehtag
      common /jehiggs/ jehwar,prlev,jehtag
c---DarkSusy addition end

      if (first) then 
         write(*,*) ('-',k=1,70)
         write(*,*) 
     &     'DS WARNING: You have requested to use FeynHiggsFast',
     &     ' (higloop=6)'
         write(*,*)
     &     'for the Higgs mass calculations. This code is obsolete and'
         write(*,*)
     &     'we strongly recommend that you use FeynHiggs (higloop=5)'
         write(*,*)
     &     'instead. This warning will only be printed once.'
         write(*,*) ('-',k=1,70)
         first=.false.
      endif

c---Setup
      jehwar=0
      jehtag=idtag

      prlev=0    ! 0=don't print things, 1=do print things

      alphadiag2=0.0d0

c ----------------------------------------------------------------
c setting the SM parameters 
c ----------------------------------------------------------------


      pi = 3.1415926535897d0

c pat scott 08-11-03
c added updated input data

      mz  = 91.1876d0
      mw =  80.39d0     ! The exp. values of mz and mw are the ones to
                        ! use for the corrections to be correct
      gf = 1.16637d-5
      as = .1095d0
      alphasmz = 0.1185d0
      cf = 4d0/3d0
      eps = 1d-10
      i = (0d0, 1d0)

      mt = inmt
      mb = inmb
      mb=2.97d0  ! Use running mb instead (JE 2000-04-25)
c      mmt = inmt   !lbe DarkSusy change
c      mbb = inm    !lbe DarkSusy change

      delrholimit = 1.3d-3


c -----------------------------------------------------------------
c
c  switches
c  ========
c
c  selec2 = 1: full 1-loop + 2-loop QCD with top pole mass
c           2: same as 1, but in addition with
c                         Yukawa term added for light Higgs
c           3: full 1-loop + 2-loop QCD with running top mass
c           4: same as 3, but in addition with
c                          Yukawa term added for light Higgs
c              The choice for the top mass (pole/running) will
c              also be applied on the top mass appearing in
c              the stop mass matrix

      selec2=4

c
c  selec3 = 1: input is Msusy, MtLR, ...
c              output is MSt1, MSt2, MSb1, MSb2, stt, stb, ...
c           2: input is MSt2, delmst, stt
c              output is MSt1, MSb1, MSb2, stb, ...
c

      selec3=1

c -----------------------------------------------------------------

c...Transfer MSSM parameters
      tb = tanb
      msusy = 1.d10   ! pg Darksusy change (not necessary as input)
      msusytl = mqtl  ! pg Darksusy change
      msusytr = mqtr  ! pg Darksusy change
      msusybl = mqbl  ! pg Darksusy change
      msusybr = mqbr  ! pg Darksusy change
      mue = my
      ma=mh3
      mtlr=ATop-my/tanb      ! M_tilde_LR
      mblr=ABot-my*tanb      ! M_tilde_LR

      if (prlev.ne.0) then
        write(*,*)
        write(*,*) "Your parameters:"
        if (selec3.eq.1) then
           write(*,*) "Msusy(top-left), Msusy(top-right), MtLR"
        else
           write(*,*) "MSt1, MSt2, stt"
        endif
        write(*,*) "tb, MT, MA, Mue"
        if (selec3.eq.1) then
           write(*,*) real(msusytl), real(msusytr), real(mtlr)
        else
           write(*,*) real(mst2 - delmst), real(mst2), real(stt)
        endif
        write(*,*) real(tb), real(mt), real(ma), real(mue)
        write(*,*)
      endif

      call feynhiggsfastsub()

      if (prlev.ne.0) then
        write(*,*) "-------------------------------------------------"
        write(*,*) "-------------------------------------------------"
        write(*,*) "The results: light Higgs    heavy Higgs     alpha"
        write(*,*) "-------------------------------------------------"
        if ((real(mst1).ge.65d0).and.(real(msb2).ge.65d0)) then
          write(*,*) "mh-tree : ", real(mlh),real(mhh), real(alpha)
          write(*,*) "-------------------------------------------------"
          write(*,*) "mh-1loop: ", real(mlhlle1) ,
     $           "                              (formula for MA >> MZ)"
          write(*,*) "mh-1loop: ", real(mlhllediag1),real(mhhllediag1),
     $                            real(alphadiag1),
     $              "(matrix diagonalisation)"
          write(*,*) "-------------------------------------------------"
          write(*,*) "mh-2loop: ", real(mlhlle2) ,
     $           "                              (formula for MA >> MZ)"
          write(*,*) "mh-2loop: ", real(mlhllediag2),real(mhhllediag2),
     $                            real(alphadiag2),
     $              "(matrix diagonalisation)"
          write(*,*) "-------------------------------------------------"
          write(*,*) "charged Higgs, tree: ", real(mhptree)
          write(*,*) "             1-loop: ", real(mhp)
        else
          write(*,*) "unphysical sfermion masses:"
          write(*,*) "MSt1, MSt2: ", real(mst1), real(mst2)
          write(*,*) "MSb1, MSb2: ", real(msb1), real(msb2)
        endif
      endif

      if (real(mst1).lt.65d0.or.real(msb2).lt.65d0) then
          jehwar=ibset(jehwar,7)
          hwar=jehwar
          goto 150   ! exit
      endif

      if (prlev.ne.0) then
        write(*,*) "-------------------------------------------------"
        write(*,*) "-------------------------------------------------"
      endif

      call delrho_fast(mst1, mst2, stt, msb1, msb2, stb, gf, as, 
     $            cf, gs, el, mz, mw,
     $            delrho1loop, delrho2loopgluon) ! ,delrho2loopgluino)
      delrhores = delrho1loop + delrho2loopgluon
      if (prlev.ne.0) then 
        if (dabs(delrhores).gt.delrholimit) then
           write(*,*) 'WARNING: Delta rho > experimental limit'
        endif
        write(*,*) 'Delta rho 1-loop         : ', delrho1loop
        write(*,*) 'Delta rho 2-loop (gluon) : ', delrho2loopgluon
        write(*,*) 'Delta rho total          : ', delrhores
        write(*,*) "-------------------------------------------------"
      endif

c...Give back results
      mh = real(mlhllediag2)    ! H_2 mass
      HM = real(mhhllediag2)    ! H_1 mass
      hc = real(mhp)            ! H+- mass
      halpha = real(alphadiag2) ! alpha
      drho=delrhores

      hwar=jehwar
      if (and(hwar,64).eq.64) then  ! 2-loop Higgs sector not OK
        mh=0.0d0
        HM=0.0d0
        hc=0.0d0
        halpha=0.0d0
        drho=0.0d0
      endif

      return

c...Take care of errors - we only get here from error goto's
 150  continue
      mh=0.0d0
      HM=0.0d0
      hc=0.0d0
      halpha=0.0d0
      drho=0.0d0
      return

      end


c------------------------------------------------------------------------
c      SUBROUTINES FOR FeynHiggsFast
c------------------------------------------------------------------------

      subroutine feynhiggsfastsub()

      implicit real*8(a-z)
      double precision mlhlle1, mlhlle2, 
     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $       alphadiag1, alphadiag2, mhptree, mhp
      integer selec2,selec3
      double precision mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                 el, gf
      double precision cf, eps, pi
      complex*16 i
      double precision  MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma

      common /smpara/ mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                el, gf
      common /para/ cf, eps, pi, i
      common /susypara_fast/ MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma
      common/higgsmass/mhh, mlh
      common /selec/selec2, selec3
      common /lle_fast/ mlhlle1, mlhlle2,
     $             mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $             alphadiag1, alphadiag2, mhptree, mhp

c---DarkSusy addition begin
      integer jehwar,prlev
      character*12 jehtag
      common /jehiggs/ jehwar,prlev,jehtag
c---DarkSusy addition end


c      if((dabs(mue)/1000d0 + tb/25d0).gt.2d0)then
c         write(*,*) 'WARNING - WARNING - WARNING - WARNING'
c         write(*,*) 'Mue-tan(beta) combination enforces m_b = 0'
c         write(*,*) 'WARNING - WARNING - WARNING - WARNING'
c         mb = 0.001
c      endif


      mz2 = mz**2
      mw2 = mw**2
      sw2 = 1.d0 - (mw/mz)**2
      el = dsqrt(8d0*gf*MW**2*(1-MW**2/MZ**2)/dsqrt(2d0))
      gs = dsqrt(4*pi*as)
      gmue = gf

      mb2 = mb**2

c      write(*,*) 'mtrun:', alphasmz, mt, mmz, pi
      as = alphasmz/(1d0 + (11d0 - 10d0/3d0)/(4d0*pi) * alphasmz
     $                       * dlog(mt**2/mz**2))
      mtrun = mt/(1d0 + 4d0/(3d0*pi) * as)

      if (selec2.ge.3) then
         mt = mtrun
         if (prlev.ne.0) write(*,*) "using running top mass:", mtrun
      endif

      if (selec3.eq.2) then
         call def3_fast(tb,mst2,delmst,stt,mw,mz,mt,
     $             msusytl,msusytr,msusybl,msusybr,mtlr)
      endif

      au = mtlr + mue/tb
c      ad = au
c --> note: At = Ab is an arbitrary choice
c      mblr = ad - mue*tb

c...JE addition - ad different from au
      ad = mblr + mue*tb

      beta = datan(tb)

      mlh2 = 0.5d0*(ma**2+mz**2 - dsqrt((ma**2+mz**2)**2 - 4.d0*
     &              mz**2*ma**2*dcos(2.d0*beta)**2))
      mlh = dsqrt(mlh2)
      Mhh2 = 0.5d0*(ma**2+mz**2 + dsqrt((ma**2+mz**2)**2 -
     &              4.d0*mz**2*ma**2*dcos(2.d0*beta)**2))
      mhh = dsqrt(mhh2)

      alpha = dacos(-dcos(2.d0*beta)*(ma**2-mz**2)/(mhh2-mlh2))/2.d0
      if (alpha.gt.0.d0) then
       alpha = -alpha
      endif

      mhptree = dsqrt(ma**2 + mw**2)

      call def2_fast
c$$$      write(*,*) '-------------------------------------------'
c$$$      write(*,*) 'Sfermions:'
c$$$      write(*,*) 'MSt1, MSt2, stt'
c$$$      write(*,*) mst1, mst2, stt
c$$$      write(*,*) 'MSb1, MSb2, stb'
c$$$      write(*,*) msb1, msb2, stb
c$$$      write(*,*) '-------------------------------------------'


c...JE addition to avoid errors
      if (mst1.le.0.0d0.or.mst2.le.0.0d0) then
        jehwar=ibset(jehwar,7)
        return
      endif

      if (msb1.le.0.0d0.or.msb2.le.0.0d0) then
        jehwar=ibset(jehwar,7)
        return
      endif
c -------------------------------------------------------------

      mtold = mt
c      write(*,*) 'vor lle'
c$$$      write(*,*) 'LLE parameters:'
c$$$      write(*,*) real(msusytl), real(msusytr), real(mtlr), real(mtold),
c$$$     $           real(mt), real(mb), real(ma), real(tb), real(gf),
c$$$     $           real(as), real(cf), real(el), real(mz), real(mw)
      call mhiggslle_fast(msusytl, msusytr, mtlr, mblr, mtold, mt, mb,
     $               mst1, mst2, stt,
     $               ma, tb, alpha, gf, as, cf, el, mz, mw, mue,
     $               mlhlle1, mlhlle2, 
     $               mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $               alphadiag1, alphadiag2, mhp)
c      write(*,*) 'nach lle'

c -------------------------------------------------------------



 999  continue


      end


c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine def2_fast

      double precision mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                 el, gf
      double precision cf, eps, pi
      complex*16 i
      double precision  MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma

      common /smpara/ mw, mz, sw2, mt, mtrun, mb, as, alphasmz, gs, 
     $                el, gf
      common /para/ cf, eps, pi, i
      common /susypara_fast/ MSt1, MSt2, stt, ctt, delmst,
     $                  MSb1, MSb2, stb, ctb,
     $                  msusy, msusytl, msusytr, msusybl, msusybr,
     $                  mtlr, mblr, au, ad, 
     $                  mmm, mgl, mue, 
     $                  beta, alpha, tb, ma

      double precision ml2t, mr2t, ml2b, mr2b, sgnt, sgnb, tt, bb
     $               , m12, m22, b, c2b

      b = datan(tb)
      c2b = dcos(2d0 * b)

c      write(*,*) 'def2_fast:'
c      write(*,*) msusytl, msusytr, msusybl, msusybr, mz, sw2, mt,mb,c2b
      ml2t = msusytl**2 + mz**2*c2b*( .5d0 - 2d0/3d0*sw2) + mt**2
      mr2t = msusytr**2 + mz**2*c2b*(        2d0/3d0*sw2) + mt**2
      ml2b = msusybl**2 + mz**2*c2b*(-.5d0 + 1d0/3d0*sw2) + mb**2
      mr2b = msusybr**2 + mz**2*c2b*(      - 1d0/3d0*sw2) + mb**2
c      write(*,*) 'diagonal entries:'
c      write(*,*) ml2t, mr2t
      sgnt = (ml2t - mr2t)/dabs(ml2t - mr2t)
      sgnb = (ml2b - mr2b)/dabs(ml2b - mr2b)

      tt = datan(-sgnt*2d0*mt*mtlr/
     $           (-sgnt*(ml2t - mr2t) 
     $            -dsqrt((ml2t - mr2t)**2 + 4d0 * mt**2 * mtlr**2)))

      stt = dsin(tt)
      ctt = dcos(tt)

      bb = datan(-sgnb*2d0*mb*mblr/
     $           (-sgnb*(ml2b - mr2b) 
     $            -dsqrt((ml2b - mr2b)**2 + 4d0 * mb**2 * mblr**2)))

      stb = dsin(bb)
      ctb = dcos(bb)

      m12 = .5d0*(ml2t + mr2t 
     $            + sgnt * dsqrt((ml2t - mr2t)**2 + 4d0*mt**2*mtlr**2))
      m22 = .5d0*(ml2t + mr2t 
     $            - sgnt * dsqrt((ml2t - mr2t)**2 + 4d0*mt**2*mtlr**2))

      
      if ((m12.lt.0d0) .or. (m22.lt.0d0))  then
         MSt1 = 0d0
         MSt2 = 0d0
         goto 100
      endif
      MSt1 = dsqrt(m12)
      MSt2 = dsqrt(m22)

 100  continue

      m12 = .5d0*(ml2b + mr2b 
     $            + sgnb * dsqrt((ml2b - mr2b)**2 + 4d0*mb**2*mblr**2))
      m22 = .5d0*(ml2b + mr2b 
     $            - sgnb * dsqrt((ml2b - mr2b)**2 + 4d0*mb**2*mblr**2))

      
      if ((m12.lt.0d0) .or. (m22.lt.0d0))  then
         MSb1 = 0d0
         MSb2 = 0d0
         goto 200
      endif
      MSb1 = dsqrt(m12)
      MSb2 = dsqrt(m22)

 200  continue

c$$$      write(*,*) "Msusy    :",msusy
c$$$      write(*,*) "Mlr      :",mlr
c$$$      write(*,*) "MT       :", MT
c$$$      write(*,*) "theta_top:",tt
c$$$      write(*,*) "-Pi/4    :",-Pi/4d0
c$$$      write(*,*) real(msusytl), real(msusytr), 
c$$$     $           real(msusybl), real(msusybr)
c$$$      write(*,*) real(mtlr), real(mblr)
c$$$      write(*,*) 'Squark masses:'
c$$$      write(*,*) real(MSt1), real(MSt2), real(stt)
c$$$      write(*,*) real(MSb1), real(MSb2), real(stb)


      end

c-------------------------------------------------------------------

      subroutine def3_fast(tb,mst2,delmst,stt,mw, mz, mt,
     $                xmsusytl,xmsusytr,xmsusybl,xmsusybr,xmtlr)

      double precision mst2,delmst,stt,ctt
      double precision xmsusytl,xmsusytr,xmsusybl,xmsusybr,xmtlr,
     $                 xmsusytl2, xmsusytr2
      double precision mw,mz,mt,tb,b,c2b,sw2,dt1,dt2
      double precision rot(1:2,1:2), rottrans(1:2,1:2), 
     $                 mdiag(1:2,1:2), mnondiag(1:2,1:2), z(1:2,1:2)

c$$$      write(*,*) 'Parameters in def3_fast:'
c$$$      write(*,*) tb, mst2, stt, delmst, mw, mz, mt

      b = datan(tb)
      c2b = dcos(2d0*b)
      sw2 = 1d0 - mw**2/mz**2
      dt1 = mz**2 * c2b * (.5d0 - 2d0/3d0 * sw2)
      dt2 = mz**2 * c2b * (       2d0/3d0 * sw2)

      mdiag(1,1) = (mst2 - delmst)**2
      mdiag(1,2) = 0d0
      mdiag(2,1) = 0d0
      mdiag(2,2) = mst2**2
      ctt = dsqrt(1d0 - stt**2)
      rot(1,1) = ctt
      rot(1,2) = stt
      rot(2,1) = -stt
      rot(2,2) = ctt
      rottrans(1,1) = ctt
      rottrans(1,2) = -stt
      rottrans(2,1) = stt
      rottrans(2,2) = ctt

      call matmult_fast(rottrans, mdiag, z)
      call matmult_fast(z, rot, mnondiag)

      xmsusytl2 = mnondiag(1,1) - mt**2 - dt1
      xmsusytr2 = mnondiag(2,2) - mt**2 - dt2
      if ((xmsusytl2.ge.0d0).and.(xmsusytr2.ge.0d0)) then
         xmsusytl = dsqrt(mnondiag(1,1) - mt**2 - dt1)
         xmsusytr = dsqrt(mnondiag(2,2) - mt**2 - dt2)
         xmsusybl = xmsusytl
         xmsusybr = xmsusybl
         xmtlr = mnondiag(1,2)/mt
      else
         xmsusytl = 0d0
         xmsusytr = 0d0
         xmsusybl = 0d0
         xmsusybr = 0d0
         xmtlr = 0d0
      endif

c      write(*,*) 'parameters from def3_fast:'
c      write(*,*) xmsusytl, xmsusytr, xmsusybl, xmtlr

      end


c-------------------------------------------------------------------

      subroutine matmult_fast(a, b, c)

      double precision a(1:2,1:2), b(1:2,1:2), c(1:2,1:2)

      c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
      c(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
      c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
      c(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

      end

c------------------------------------------------------------------------

      subroutine diagonalization_fast(a, b, c, d)

      double precision a(1:2,1:2), b(1:2,1:2), c(1:2,1:2), d(1:2,1:2),
     $                 e(1:2,1:2)

      call matmult_fast(a, b, e)
      call matmult_fast(e, c, d)

      end

c------------------------------------------------------------------------

c------------------------------------------------------------------------


c     subroutine for the fast (and approximate) calculation of the 
c     neutral CP-even Higgs-boson masses in the MSSM
c     Two approaches are used:
c     1.) express mlh^2 through a closed formula. Here problems might 
c         arise for small MA
c     2.) the Higgs mass matrix ist determined numerically up to 2-loop,
c         the diagonalisation yields the two Higgs masses

      subroutine mhiggslle_fast(msusytl, msusytr, mtlr, mblr, mt,
     $           mtrun, mb, 
     $           mst1, mst2, stt,
     $           ma, tb, alpha, gf, as, cf, el, mz, mw, mue,
     $           mlhlle1, mlhlle2, 
     $           mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $           alphadiag1, alphadiag2, mhp)

c     meaning of variables:
c     msusyt[l,r]: SSB terms in the stop mass matrix
c     ms[1,2]: The Susy mass scale: ms1^2 = Msusy^2 + MT^2 (for 1-loop contr.)
c                                   ms2^2 = Msusy^2 + MTrun^2 (2-loop contr.)
c     ms[1,2]2: == ms[1,2]^2
c     mtlr: the off-diagonal entry in the stop mass matrix
c     mt, mb: fermion pole masses
c     mtrun: running top mass, can be used for 2-loop contribution
c     ma, tb: the usual Higgs sector varialbes M_A and tan(beta)
c     gf: Fermi's constant
c     as: strong coupling constant at the scale mt
c     cf = 4/3
c     el: electromagnetic coupling, derived from gf
c     mz, mw: gauge boson masses
c     mlhlle[1,2]: the mass of the light Higgs boson in approach 1
c     m[l,h]hllediag[1,2]: the Higgs masses in approach 2


      double precision msusy, msusytl, msusytr, mtlr, mblr, mt, mtrun, 
     $                 mb, ma, theta,
     $                 gf, as, pi, cf, el, mw, mz, pref, al,
     $                 sw, sw2, sw4, cw2, et, eb, mue,
     $                 ma2, mz2, mw2, mt2, mt4, msfkt_fast,
     $                 msusy2, ms1, ms2, ms12, ms22,
     $                 alpha, tb, b, sb, cb, sb2, cb2, c2b
      double precision MSt1, MSt2, stt
      double precision nc, ng, pt, pb, pf, pg, pgp, p1h, p2h, p2hp
      double precision mlhlle1, mlhlle2, mlhlle1sq, mlhlle2sq, 
     $       mlhllediag1, mlhllediag2, mhhllediag1, mhhllediag2,
     $       mlhllediag1sq, mlhllediag2sq, mhhllediag1sq, mhhllediag2sq,
     $       alphadiag1, alphadiag2, mhp
      double precision m11, m12, m221loop, m222loop
      double precision m11botmix, m12botmix, m22botmix, 
     $                 m11top, m12top, m22top, m22log, 
     $                 m11rest, m12rest, m22rest
      double precision xttilde, vvv, ttt, delmlhsq
      double precision rot(1:2,1:2),rott(1:2,1:2),mhnondiag(1:2,1:2),
     $                 mhdiag(1:2,1:2)
      integer selec2, selec3
      common /selec/selec2, selec3

c---DarkSusy addition begin
      integer jehwar,prlev
      character*12 jehtag
      common /jehiggs/ jehwar,prlev,jehtag
c---DarkSusy addition end


c      write(*,*) 'LLE parameters:'
c      write(*,*) mt, mtrun
c      write(*,*) mb
c      write(*,*) real(mtlr), real(mblr), real(mue)

      pi = 3.14159265358979d0
      al = el**2/(4d0 * pi)
      cw2 = mw**2/mz**2
      sw2 = 1d0 - mw**2/mz**2
      sw = dsqrt(sw2)
      sw4 = sw2**2
      et = 2d0/3d0
      eb = -1d0/3d0
      b = datan(tb)
      sb = dsin(b)
      sb2 = sb**2
      cb = dcos(b)
      cb2 = cb**2
      c2b = dcos(2d0 * b)
      ma2 = ma**2
      mz2 = mz**2
      mw2 = mw**2
      mt2 = mt**2
      mt4 = mt**4
      if (msusytl.eq.msusytr) then
         ms1 = dsqrt(msusytl**2 + mt**2)
         ms2 = dsqrt(msusytl**2 + mtrun**2)
      else
         ms1 = msfkt_fast(msusytl, msusytr, mt)
         ms2 = msfkt_fast(msusytl, msusytr, mtrun)
      endif
      ms12 = ms1**2
      ms22 = ms2**2
      msusy2 = ms12 - mt**2
      msusy = dsqrt(msusy2)

      nc = 3d0
      ng = 3d0
      pt = nc * (1d0 - 4d0 * et * sw2 + 8d0 * et**2 * sw4)
      pb = nc * (1d0 + 4d0 * eb * sw2 + 8d0 * eb**2 * sw4)
      pf = nc * (ng - 1d0) * (2d0 - 4d0 * sw2 
     $           + 8d0 * (et**2 + eb**2) * sw4) 
     $           + ng * (2d0 - 4d0 * sw2 + 8d0 * sw4)
      pg = -44d0 + 106d0 * sw2 - 62d0 * sw4
      pgp = 10d0 + 34d0 * sw2 - 26d0 * sw4
      p1h = -9d0 * c2b**4 + (1d0 - 2d0 * sw2 + 2d0 * sw4) * c2b**2
      p2h = -10d0 + 2d0 * sw2 - 2d0 * sw4
      p2hp = 8d0 - 22d0 * sw2 + 10d0 * sw4

      theta = 0d0
      if (ma2.ge.mz2) theta = 1d0


c --> 1. approach: the 'direct' formula for mh2

      mlhlle1sq =  
     $    (MA**2 + MZ**2 - Sqrt(MA**4 + MZ**4 - 
     -        2*MA**2*MZ**2*(1 - 8*SB**2 + 8*SB**4)))/2. + 
     -   (al*MT**4*((-3*MT**4*MtLR**8)/(55.99999999999999*Ms1**12) + 
     n    (1d0/4d0 * MZ**2/MT**2 - 11d0/80d0 * MZ**4/MT**4)
     -       + (MtLR**2*(1.5 - (3*MT**2)/(4.*Ms1**2) + 
     $                       0d0 - 
     -             MZ**2/(2.*MT**2)))/Ms1**2 + 
     -     (MtLR**6*((-3*MT**2)/(40.*Ms1**2) + (3*MT**4)/(10.*Ms1**4) - 
     -                                    0d0  ))/Ms1**6 + 
     -    (MtLR**4*(-0.125 + MT**2/(2.*Ms1**2) - (3*MT**4)/(8.*Ms1**4)-
     -                                0d0 + 
     $                             0d0 ))/Ms1**4
     -        )*(1 + ((-MA**2 + MZ**2)*(1 - 2*SB**2))/
     -         Sqrt(MA**4 + MZ**4 - 2*MA**2*MZ**2*
     $     (1 - 8*SB**2 + 8*SB**4))))/
     -    (2.*MW**2*Pi*SB**2*sw2)  

      mlhlle1sq = mlhlle1sq -
     -   (al*MT**4*(1.5 - (3*MZ**2*SB**2)/(4.*MT**2) + 
     -        ((32*MW**4 - 40*MW**2*MZ**2 + 17*MZ**4)*
     $     SB**2)/(72.*MT**4) + 
     -        (4*(MA**2 + MZ**2)*SB**2*(1 - SB**2)*
     -           ((3*MZ**2)/(8.*MT**2) - 
     -             ((32*MW**4 - 40*MW**2*MZ**2 + 17*MZ**4)*SB**2)/
     $     (72.*MT**4))
     -           )/Sqrt(MA**4 + MZ**4 - 2*MA**2*MZ**2*
     $     (1 - 8*SB**2 + 8*SB**4))
     -          + ((-MA**2 + MZ**2)*(1 - 2*SB**2)*
     -           (1.5 - (3*MZ**2*SB**2)/(4.*MT**2) - 
     -             ((32*MW**4 - 40*MW**2*MZ**2 + 17*MZ**4)*SB**2*
     -                (1 - 2*SB**2))/(72.*MT**4)))/
     -     Sqrt(MA**4 + MZ**4 - 2*MA**2*MZ**2*
     $     (1 - 8*SB**2 + 8*SB**4)))*
     -      Log(MT**2/Ms1**2))/(2.*MW**2*Pi*SB**2*sw2)

      mlhlle1sq = mlhlle1sq +
     -   (al*MZ**4*(Log(Msusy**2/MZ**2)*
     -         ((36*MB**4)/MZ**4 - (18*MB**2*Cos(2*b))/MZ**2 + 
     -           (3*(1 - (4*sw2)/3. + (8*sw4)/9.) + 
     -              6*(2 - 4*sw2 + (40*sw4)/9.) + 
     -              3*(2 - 4*sw2 + 8*sw4))*Cos(2*b)**2 - 
     -           2*(18 + 12*sw2 - 16*sw4)*Cos(b)**2*Sin(b)**2 + 
     -           (-54 + 108*sw2 - 64*sw4)*
     $     (Cos(b)**4 + Sin(b)**4)) - 
     -        theta * Log(MA**2/MZ**2)*(-((1 - 2*sw2 + 2*sw4)*
     $     Cos(2*b)**2) + 
     -           9*Cos(2*b)**4 - 
     -           2*(8 - 22*sw2 + 10*sw4)*Cos(b)**2*Sin(b)**2 + 
     -           (-10 + 2*sw2 - 2*sw4)*
     $     (Cos(b)**4 + Sin(b)**4))))/
     -    (23.99999999999999*MW**2*Pi*sw2)

      mlhlle1sq = mlhlle1sq +
     -   (3*al*((2*MB**4*MbLR**2*(1 - MbLR**2/
     $     (11.99999999999999*Msusy**2)))/
     -         Msusy**2 - (MB**2*MZ**2*Cos(2*b)*
     -           (MbLR**2 + (-(Mue**2*Tan(b)**2) + 
     $     (MbLR + Mue*Tan(b))**2)/3.)
     -           )/(2.*Msusy**2)))/(4.*MW**2*Pi*sw2)
      if (mlhlle1sq.gt.0d0) then
         mlhlle1 = dsqrt(mlhlle1sq)
      else
         jehwar=ibset(jehwar,2)
         mlhlle1 = 119.9999d0
      endif


c      write(*,*) 'lle 2loop'
      if (selec2.ge.3) then
c         write(*,*) '1-loop with mtrun => changed 2-loop formula'
      mlhlle2sq = (
     $     mlhlle1sq
     -   + (al*as*MT**4*(1 + ((-MA**2 + MZ**2)*(1 - 2*SB**2))/
     -         Sqrt(MA**4 + MZ**4 - 2*MA**2*MZ**2*
     $     (1 - 8*SB**2 + 8*SB**4)))*
     -      ((6*MtLR)/Ms2 - (12*MtLR**4)/(17.*Ms2**4) + 
     -        (+2d0 - (3*MtLR**2)/Ms2**2)*Log(Ms2**2/MT**2) - 
     -        3*Log(Ms2**2/MT**2)**2
     $      - 4d0 + 8d0 * MtLR**2/Ms2**2))/
     $     (2.*MW**2*Pi**2*SB**2*sw2) 
     $     )
      else
c         write(*,*) '1-loop with mt pole => original 2-loop formula'
      mlhlle2sq = (
     $     mlhlle1sq
     -   + (al*as*MT**4*(1 + ((-MA**2 + MZ**2)*(1 - 2*SB**2))/
     -         Sqrt(MA**4 + MZ**4 - 2*MA**2*MZ**2*
     $     (1 - 8*SB**2 + 8*SB**4)))*
     -      ((6*MtLR)/Ms1 - (3*MtLR**4)/(4.*Ms1**4) + 
     -        (-6 - (3*MtLR**2)/Ms1**2)*Log(Ms1**2/MT**2) - 
     -        3*Log(Ms1**2/MT**2)**2))/
     $     (2.*MW**2*Pi**2*SB**2*sw2) 
     $     )
      endif

      if (mlhlle2sq.gt.0d0) then
         mlhlle2 = dsqrt(mlhlle2sq)
      else
         jehwar=ibset(jehwar,3)
         mlhlle2 = 119.9999d0
      endif


c ----------------------------------------------------------------------

c --> 2. approach: the Higgs mass matrix


c      write(*,*) '2. approach'
      pref = el**2 * mz2 / (96d0 * Pi**2 * cw2 * sw2)
      mwmzfak = 32d0 * mw**4 - 40d0 * mw2 * mz2 + 17d0 * mz**4

c      write(*,*) '2. approach'
c --> the top sector at 1-loop
      m11top = -el**2/(288d0 * mw2 * pi**2 * sw2) 
     $         * mwmzfak * cb2 * dlog(mt2/ms12)
      m12top = el**2 * cb/(288d0 * sb * mw2 * pi**2 * sw2) 
     $         * (- 27 * mt2 * mz2 + mwmzfak * sb2) * dlog(mt2/ms12)
      m22top = el**2/(288d0 * sb2 * mw2 * pi**2 * sw2) *
     $  ((-108d0 * mt**4 + 54d0 * mt2 * mz2 * sb2 - mwmzfak * sb**4) *
     $     dlog(mt2/ms12)
     n     + (2d0 * MZ**2/MT**2 - 11d0/10d0 * MZ**4/MT**4) * 9d0 * MT**4
     $     +9d0 * (MtLR/Ms1)**2 * (MT/Ms1)**2 *
     $ (-6d0 * mt**4                   + 4d0 * ms12 * (3d0 * mt2 - mz2))
     $     +9d0/10d0 * (MtLR/Ms1)**4 * (MT/Ms1)**4 *
     $     (-10d0 * ms1**4 + 4 * ms1**2 * (10d0 * mt2            ) +
     $     15d0 * mt2 * (-2d0 * mt2      ))
     $     +9d0/35d0 * (MtLR/Ms1)**6 * (MT/Ms1)**6 *
     $     (-21d0 * ms1**4 + 12d0 * ms1**2 * (7d0 * mt2            ) +
     $     35d0 * mt2 * (-2d0 * mt2      ))
     $     -288d0/469762048d0 * (MtLR/Ms1)**8 * (MT/Ms1)**8 *
     $     (6291456d0 * ms1**4 + ms1**2 * (-25165824d0 * mt2
     $                      ) + 11010048d0 * (2d0 * mt**4            )))
      

c      write(*,*) '2. approach'
c --> all other sectors except for the sbottom mixing
      m11rest = pref * cb2 *
     $          ((  12d0 * nc * mb**4/(mz**4 * cb**4) 
     $           - 6d0 * nc * mb**2/(mz**2 * cb**2) 
     $           + pb + pf + pg + p2h ) * dlog(msusy2/mz2)
     $          + theta * (p1h - p2h) * dlog(ma2/mz2) )
      m12rest = -pref * sb * cb *
     $          ((-3d0 * nc * mb**2/(mz2 * cb2) 
     $            + pb + pf + pgp + p2hp) * dlog(msusy2/mz2) 
     $           - theta * (p1h + p2hp) * dlog(ma2/mz2) )
      m22rest = pref * sb2 * 
     $          ((pb + pf + pg + p2h) * dlog(msusy2/mz2) 
     $           + theta * (p1h - p2h) * dlog(ma2/mz2) )
       
        
c      write(*,*) '2. approach'
c --> the bottom mixing at 1-loop 
      m11botmix = 
     &  (AL*MW**(-2)*NC/SW2*(-(MB**2*MSUSY**(-2)*MZ**2*(MTLR+
     &  MUE/tan(B))*(MBLR+((MTLR+MUE/tan(B)))/3D0))+4D0*CB**(-2)*
     &  MB**4*MBLR*MSUSY**(-2)*(MTLR+MUE/tan(b))*(1D0-(MBLR*MSUSY**
     &  (-2)*(MTLR+MUE/tan(b)))/12D0)))/(8D0*PI)
      m12botmix =
     &  -((AL*MW**(-2)*NC/SW2*(4D0*CB**(-2)*MB**4*MBLR*MSUSY**
     &  (-2)*MUE*(1D0-(MBLR*MSUSY**(-2)*(MTLR+MUE/tan(b)))/6D0)-MB**
     &  2*MSUSY**(-2)*MZ**2*(MBLR*(MTLR+2D0*MUE/tan(b))+((MUE**2+
     &  (MTLR+MUE/tan(b))**2))/3D0)*TAN(B)))/(16D0*PI))
      m22botmix = 
     &  (AL*MW**(-2)*NC/SW2*(-((1D0*CB**(-2)*MB**4*MBLR**2*
     &  MSUSY**(-4)*MUE**2)/3D0)-MB**2*MSUSY**(-2)*MUE*MZ**2*TAN(B)*
     &  (MBLR+(MUE*TAN(B))/3D0)))/(8D0*PI)


c --> the top sector at 2-loop
      if (selec2.ge.3) then
c         write(*,*) '1-loop with mtrun => changed 2-loop formula'
         m22log = (al*as*MTrun**4/dsin(b)**2 * 
     -      ((6*MtLR)/Ms2 - (12*MtLR**4)/(17.*Ms2**4) + 
     -        (+2d0 - (3*MtLR**2)/Ms2**2)*dLog(Ms2**2/MT**2) - 
     -        3*dLog(Ms2**2/MT**2)**2
     $      - 4d0 + 8d0 * MtLR**2/Ms2**2))/
     $        (MW**2*Pi**2*SW2)
      else
c         write(*,*) '1-loop with mt pole => original 2-loop formula'
         m22log = (al*as*MT**4/dsin(b)**2 * 
     $        ((6*MtLR)/Ms2 - (3*MtLR**4)/(4.*Ms2**4) + 
     $        (-6 - (3*MtLR**2)/Ms2**2)*dLog(Ms2**2/MT**2) 
     $        - 3*dLog(Ms2**2/MT**2)**2
     $        ))/(MW**2*Pi**2*SW2)
      endif

c$$$      write(*,*) 'Parameters in LLESub:'
c$$$      write(*,*) 'Msusy, MT, MtLR,'
c$$$      write(*,*) 'M22top1loop, M22rest, M22botmix, M22top2loop'
c$$$      write(*,*) real(msusy), real(mt), real(mtlr)
c$$$      write(*,*) real(m22top), real(m22rest), real(m22botmix), 
c$$$     $           real(m22log)
c$$$      write(*,*) cf, el, al, sw, pi, sb, mtrun, mw, llexpansionp2mtrun

c     including the leading Yukawa term for the light Higgs mass
c     this term is taken from Carena, Espinosa, Quiros, Wagner
c     Nucl. Phys. B461 (1996) 407

      if (prlev.ne.0) write(*,*) 'adding Yukawa term'
      xttilde = (((MSt2**2 - MSt1**2)/(4d0*mt**2) * 
     $            (2d0 * stt * dsqrt(1d0 - stt**2))**2 )**2 *
     $           (2d0 - (MSt2**2 + MSt1**2)/(MSt2**2 - MSt1**2) *
     $                  dlog(MSt2**2/MSt1**2)) +
     $           (MSt2**2 - MSt1**2)/(2d0 * mt**2) *
     $            (2d0 * stt * dsqrt(1d0 - stt**2))**2 *
     $                  dlog(MSt2**2/MSt1**2) )
      ttt = .5d0 * dlog(MSt2**2 * MSt1**2/mt**4)
      vvv = 174.1d0
      delmlhsq = ( 3d0/(4d0 * pi**2) * mt**4/vvv**2 *
     $            ( 1d0/(16d0 * pi**2) * 3d0/2d0 * mt**2/vvv**2 *
     $             ( xttilde * ttt + ttt**2)))


c----------------------------------------------------------------


      m11 = m11top + m11rest + m11botmix
      m12 = m12top + m12rest + m12botmix
      m221loop = m22top + m22rest + m22botmix
      m222loop = m221loop + m22log
      if ((selec2.eq.2).or.(selec2.eq.4)) then
         m222loop = m222loop + delmlhsq/sb2
      endif


c----------------------------------------------------------------


      mlhllediag1sq = 
     $     (ma2 + mz2 + m221loop + m11 -
     $      dsqrt(-8d0 * m12 * (ma2 + mz2) * sb * cb +
     $            (ma2 + mz2 + m221loop + m11)**2
     $            - 4d0 * (-m12**2 + mz2 * m221loop * cb2 
     $                     + m221loop * m11 + mz2 * sb2 * m11 
     $                     + ma2 * (mz2 - 4d0 * mz2 * sb2 
     $                     + m221loop * sb2 + 4d0 * mz2 * sb**4
     $                     + m11 * cb2))))/2d0
      mlhllediag2sq = 
     $     (ma2 + mz2 + m222loop + m11 -
     $      dsqrt(-8d0 * m12 * (ma2 + mz2) * sb * cb +
     $            (ma2 + mz2 + m222loop + m11)**2
     $            - 4d0 * (-m12**2 + mz2 * m222loop * cb2 
     $                     + m222loop * m11 + mz2 * sb2 * m11 
     $                     + ma2 * (mz2 - 4d0 * mz2 * sb2 
     $                     + m222loop * sb2 + 4d0 * mz2 * sb**4
     $                     + m11 * cb2))))/2d0
      mhhllediag1sq = 
     $     (ma2 + mz2 + m221loop + m11 +
     $      dsqrt(-8d0 * m12 * (ma2 + mz2) * sb * cb +
     $            (ma2 + mz2 + m221loop + m11)**2
     $            - 4d0 * (-m12**2 + mz2 * m221loop * cb2 
     $                     + m221loop * m11 + mz2 * sb2 * m11 
     $                     + ma2 * (mz2 - 4d0 * mz2 * sb2 
     $                     + m221loop * sb2 + 4d0 * mz2 * sb**4
     $                     + m11 * cb2))))/2d0
      mhhllediag2sq = 
     $     (ma2 + mz2 + m222loop + m11 +
     $      dsqrt(-8d0 * m12 * (ma2 + mz2) * sb * cb +
     $            (ma2 + mz2 + m222loop + m11)**2
     $            - 4d0 * (-m12**2 + mz2 * m222loop * cb2 
     $                     + m222loop * m11 + mz2 * sb2 * m11 
     $                     + ma2 * (mz2 - 4d0 * mz2 * sb2 
     $                     + m222loop * sb2 + 4d0 * mz2 * sb**4
     $                     + m11 * cb2))))/2d0


      if (mlhllediag1sq.gt.0d0) then 
         mlhllediag1 = dsqrt(mlhllediag1sq)
         alphadiag1 = datan((-(ma2 + mz2) * sb * cb + m12)/
     $                      (mz2 * cb2 + ma2 * sb2 + m11 
     $                       - mlhllediag1sq))
      else
         mlhllediag1 = 119.9999d0
         alphadiag1 = 0d0
         jehwar=ibset(jehwar,5)
      endif

      if (mlhllediag2sq.gt.0d0) then
         mlhllediag2 = dsqrt(mlhllediag2sq)
         alphadiag2 = datan((-(ma2 + mz2) * sb * cb + m12)/
     $                      (mz2 * cb2 + ma2 * sb2 + m11 
     $                       - mlhllediag2sq))
      else
         mlhllediag2 = 119.9999d0
         alphadiag2 = 0d0
         jehwar=ibset(jehwar,6)
      endif

      mhhllediag1 = dsqrt(mhhllediag1sq)
      mhhllediag2 = dsqrt(mhhllediag2sq)


c----------------------------------------------------------------

      if ((selec2.eq.2).or.(selec2.eq.4)) then
c     including the leading Yukawa term for the light Higgs mass
c     this term is taken from Carena, Espinosa, Quiros, Wagner
c     Nucl. Phys. B461 (1996) 407

      if (mlhlle2**2 + delmlhsq.gt.0.0d0) then
        mlhlle2 = dsqrt(mlhlle2**2 + delmlhsq)
      else
        jehwar=ibset(jehwar,4)
        mlhlle2=119.9999d0
      endif

c      write(*,*) real(mlhlle2), real(mlhllediag2)

      endif

c----------------------------------------------------------------

c --> rotation check
      rot(1,1) = dcos(alpha)
      rot(1,2) = -dsin(alpha)
      rot(2,1) = -rot(1,2)
      rot(2,2) = rot(1,1)
      rott(1,1) = rot(1,1)
      rott(1,2) = rot(2,1)
      rott(2,1) = rot(1,2)
      rott(2,2) = rot(2,2)
      mhnondiag(1,1) = ma2 * sb2 + mz2 * cb2
      mhnondiag(1,2) = - (ma2 + mz2) * sb * cb
      mhnondiag(2,1) = mhnondiag(1,2)
      mhnondiag(2,2) = ma2 * cb2 + mz2 * sb2
      call diagonalization_fast(rott, mhnondiag, rot, mhdiag)
      if (dabs(mhdiag(1,2)).le.1d-6) mhdiag(1,2) = 0d0
      if (dabs(mhdiag(2,1)).le.1d-6) mhdiag(2,1) = 0d0
c      write(*,*) 'tree level diagonalization:'
c      write(*,*) real(dsqrt(mhdiag(1,1))), real(mhdiag(1,2))
c      write(*,*) real(mhdiag(2,1)), real(dsqrt(mhdiag(2,2)))

      rot(1,1) = dcos(alphadiag1)
      rot(1,2) = -dsin(alphadiag1)
      rot(2,1) = -rot(1,2)
      rot(2,2) = rot(1,1)
      rott(1,1) = rot(1,1)
      rott(1,2) = rot(2,1)
      rott(2,1) = rot(1,2)
      rott(2,2) = rot(2,2)
      mhnondiag(1,1) = mhnondiag(1,1) + m11
      mhnondiag(1,2) = mhnondiag(1,2) + m12
      mhnondiag(2,1) = mhnondiag(1,2)
      mhnondiag(2,2) = mhnondiag(2,2) + m221loop
      call diagonalization_fast(rott, mhnondiag, rot, mhdiag)
      if (dabs(mhdiag(1,2)).le.1d-6) mhdiag(1,2) = 0d0
      if (dabs(mhdiag(2,1)).le.1d-6) mhdiag(2,1) = 0d0
c      write(*,*) '1-loop diagonalization:'
c      write(*,*) real(dsqrt(mhdiag(1,1))), real(mhdiag(1,2))
c      write(*,*) real(mhdiag(2,1)), real(dsqrt(mhdiag(2,2)))
      if ((dabs(mlhllediag1sq - mhdiag(2,2)).gt.1d-4).or.
     $    (dabs(mhhllediag1sq - mhdiag(1,1)).gt.1d-4).or.
     $    (dabs(mhdiag(1,2)).gt.1d-6).or.
     $    (dabs(mhdiag(2,1)).gt.1d-6)) then
         write(*,*) 'Model: ',jehtag
         write(*,*) "WARNING --- WARNING --- WARNING --- WARNING"
         write(*,*) "Error in 1-loop rotation !!!!!!!!!!!!!!!!!!"
         write(*,*) "WARNING --- WARNING --- WARNING --- WARNING"
         write(*,*) '1-loop diagonalization:'
         write(*,*) 'Rotation matrix:'
         write(*,*) real(rot(1,1)),real(rot(1,2))
         write(*,*) real(rot(2,1)),real(rot(2,2))
         write(*,*) 'Non-diagonal matrix:'
         write(*,*) real(mhnondiag(1,1)),real(mhnondiag(1,2))
         write(*,*) real(mhnondiag(2,1)),real(mhnondiag(2,2))
         write(*,*) 'Diagonalized matrix:'
         write(*,*) real(mhdiag(1,1)), real(mhdiag(1,2))
         write(*,*) real(mhdiag(2,1)), real(mhdiag(2,2))
         jehwar=ibset(jehwar,5)
      endif

      rot(1,1) = dcos(alphadiag2)
      rot(1,2) = -dsin(alphadiag2)
      rot(2,1) = -rot(1,2)
      rot(2,2) = rot(1,1)
      rott(1,1) = rot(1,1)
      rott(1,2) = rot(2,1)
      rott(2,1) = rot(1,2)
      rott(2,2) = rot(2,2)
      mhnondiag(1,1) = mhnondiag(1,1)
      mhnondiag(1,2) = mhnondiag(1,2)
      mhnondiag(2,1) = mhnondiag(1,2)
      mhnondiag(2,2) = mhnondiag(2,2) - m221loop + m222loop
      call diagonalization_fast(rott, mhnondiag, rot, mhdiag)
      if (dabs(mhdiag(1,2)).le.1d-6) mhdiag(1,2) = 0d0
      if (dabs(mhdiag(2,1)).le.1d-6) mhdiag(2,1) = 0d0
c      write(*,*) '2-loop diagonalization:'
c      write(*,*) real(dsqrt(mhdiag(1,1))), real(mhdiag(1,2))
c      write(*,*) real(mhdiag(2,1)), real(dsqrt(mhdiag(2,2)))
      if ((dabs(mlhllediag2sq - mhdiag(2,2)).gt.1d-4).or.
     $    (dabs(mhhllediag2sq - mhdiag(1,1)).gt.1d-4).or.
     $    (dabs(mhdiag(1,2)).gt.1d-6).or.
     $    (dabs(mhdiag(2,1)).gt.1d-6)) then
         write(*,*) 'Model: ',jehtag
         write(*,*) "WARNING --- WARNING --- WARNING --- WARNING"
         write(*,*) "Error in 2-loop rotation !!!!!!!!!!!!!!!!!!"
         write(*,*) "WARNING --- WARNING --- WARNING --- WARNING"
         write(*,*) '2-loop diagonalization:'
         write(*,*) 'Rotation matrix:'
         write(*,*) real(rot(1,1)),real(rot(1,2))
         write(*,*) real(rot(2,1)),real(rot(2,2))
         write(*,*) 'Non-diagonal matrix:'
         write(*,*) real(mhnondiag(1,1)),real(mhnondiag(1,2))
         write(*,*) real(mhnondiag(2,1)),real(mhnondiag(2,2))
         write(*,*) 'Diagonalized matrix:'
         write(*,*) real(mhdiag(1,1)), real(mhdiag(1,2))
         write(*,*) real(mhdiag(2,1)), real(mhdiag(2,2))
         jehwar=ibset(jehwar,6)
      endif

c --> end of rotation check

c-------------------------------------------------------------------
c --> charged Higgs at one loop
      mhp = dsqrt(ma2 + mw2 + 3d0/(4d0 * pi) * al/(sw2 * mw2) *
     $            (  2d0 * mt2 * mb**2/(sb2 * cb2) 
     $             - mw2 * (mt2/sb2 + mb**2/cb2) + 2d0/3d0 * mw2**2)
     $            * dlog(msusy/mt)
     $            + mw2/(6d0 * pi) * 15d0 * al/cw2 * dlog(msusy/mw))


      end


c-------------------------------------------------------------------

      double precision function msfkt_fast(msusytl, msusytr, mtrun)

      double precision msusytl, msusytr, mtrun

      msfkt_fast = dsqrt( dsqrt(  msusytl**2 * msusytr**2 
     $         + mtrun**2 * (msusytl**2 + msusytr**2) + mtrun**4))

      end


c-------------------------------------------------------------------




      


c------------------------------------------------------------------------
c------------------------------------------------------------------------


      subroutine delrho_fast(mst1, mst2, stt, msb1, msb2, stb,
     $           gf, as, cf, gs, el, mz, mw, 
     $           res1loop, res2loopgluon) ! , res2loopgluino)

      double precision mst1, mst2, stt, ctt,  
     $                 msb1, msb2, stb, ctb, 
     $                 as, pi, gf, cf, gs, el, mw, mz, pref
      double precision oneloop, twoloop, f0_fast, f1_fast,
     $                 res1loop, res2loopgluon

      pi = 3.14159265358979d0
      ctt = dsqrt(1d0 - stt**2)
      ctb = dsqrt(1d0 - stb**2)
      pref = (3d0 * cf * gs**2 * el**2 * mz**2)/
     $       (2**6 * mw**2 * (mw**2 - mz**2) * pi**4)

c$$$      write(*,*) 'paramters in delrho_fast:'
c$$$      write(*,*) 'mst1, mst2, stt, msb1, msb2, stb, gf, as'
c$$$      write(*,*) real(mst1), real(mst2), real(stt), real(msb1), 
c$$$     $           real(msb2), real(stb), real(gf), real(as),
c$$$     $           pref

      oneloop = 3d0 * gf/(8d0 * dsqrt(2d0) * pi**2) *
     $          ( - stt**2 * ctt**2 * f0_fast(mst1**2, mst2**2)
     $            - stb**2 * ctb**2 * f0_fast(msb1**2, msb2**2)
     $            + ctt**2 * ctb**2 * f0_fast(mst1**2, msb1**2)
     $            + ctt**2 * stb**2 * f0_fast(mst1**2, msb2**2)
     $            + stt**2 * ctb**2 * f0_fast(mst2**2, msb1**2)
     $            + stt**2 * stb**2 * f0_fast(mst2**2, msb2**2) )

      twoloop = gf * as/(4d0 * dsqrt(2d0) * pi**3) *
     $          ( - stt**2 * ctt**2 * f1_fast(mst1**2, mst2**2)
     $            - stb**2 * ctb**2 * f1_fast(msb1**2, msb2**2)
     $            + ctt**2 * ctb**2 * f1_fast(mst1**2, msb1**2)
     $            + ctt**2 * stb**2 * f1_fast(mst1**2, msb2**2)
     $            + stt**2 * ctb**2 * f1_fast(mst2**2, msb1**2)
     $            + stt**2 * stb**2 * f1_fast(mst2**2, msb2**2) )

      res1loop = oneloop
      res2loopgluon = twoloop

      end


c-------------------------------------------------------------------

      double precision function f0_fast(x,y)

      double precision x, y

      if (x.ne.y) then
c         write(*,*) 'f0_fast:', x, y
      f0_fast = x + y - (2d0 * x * y)/(x - y) * dlog(x/y)
      else
      f0_fast = 0d0
      endif

c$$$      f0_fast = x + y - (2d0 * x * y)/(x - y) * dlog(x/y)

      end


c-------------------------------------------------------------------

      double precision function f1_fast(x2,y2)

      double precision x2, y2
      complex*16 ff, x, y, cspen_fast
      
      x = (dsqrt(x2) - (0d0,1d0) * 10d-10)**2
      y = (dsqrt(y2) - (0d0,1d0) * 10d-10)**2

      if (x.ne.y) then
c         write(*,*) 'f1_fast:', x, y
      ff = x + y 
     $   - (2d0 * x * y)/(x - y) * cdlog(x/y) * (2d0 + x/y * cdlog(x/y))
     $   + (x + y)*x**2/(x - y)**2 * (cdlog(x/y))**2 
     $   - 2d0 * (x - y) * cspen_fast(1d0 - x/y)

      f1_fast = dreal(ff)
      else
      f1_fast = 0d0
      endif

c$$$      ff = x + y 
c$$$     $   - (2d0 * x * y)/(x - y) * cdlog(x/y) * (2d0 + x/y * cdlog(x/y))
c$$$     $   + (x + y)*x**2/(x - y)**2 * (cdlog(x/y))**2 
c$$$     $   - 2d0 * (x - y) * cspen_fast(1d0 - x/y)
c$$$
c$$$      f1_fast = dreal(ff)

      end


c------------------------------------------------------------------------
c------------------------------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN_FAST(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPEN_FAST,W,SUM,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)/
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
c      write(*,*) 'z:',z
      Z =Z*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c      write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
        CSPEN_FAST=-CDLOG(1D0-Z)
c        write(*,*) 'cspen_fast:', cspen_fast
        RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
        CSPEN_FAST=1.64493406684822643D0
c        write(*,*) 'cspen_fast:', cspen_fast
        RETURN
      END IF
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CDLOG(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 2
c      write(*,*) 'u:',u
c      write(*,*) 'sum:',sum
      DO 1 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2
      SUM=SUM+U*B(K)
 1    CONTINUE
 2    CSPEN_FAST=SUM
c        write(*,*) 'cspen_fast:', cspen_fast
      RETURN
10    W=-CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 12

      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12
      SUM=SUM+U*B(K)
11    CONTINUE
12    CSPEN_FAST=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
c        write(*,*) 'cspen_fast:', cspen_fast
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 22
      DO 21 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22
      SUM=SUM+U*B(K)
21    CONTINUE
22    CSPEN_FAST=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen_fast:', cspen_fast
      RETURN
30    W=CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 32
      DO 31 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32
      SUM=SUM+U*B(K)
31    CONTINUE
32    CSPEN_FAST=SUM+3.28986813369645287D0
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
50    CONTINUE
c        write(*,*) 'cspen_fast:', cspen_fast
      END


c------------------------------------------------------------------------
c------------------------------------------------------------------------


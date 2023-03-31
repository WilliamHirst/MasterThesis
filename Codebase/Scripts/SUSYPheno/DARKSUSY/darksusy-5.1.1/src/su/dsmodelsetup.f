      subroutine dsmodelsetup(unphys,hwarning)
c_______________________________________________________________________
c  set up global variables for the supersymmetric model routines.
c  If rate routines are going to be called afterwards, dsprep should
c  be called after this routine, like it is done in dssusy.
c  You should only call this routine directly yourself if you know
c  what you are doing. If you are the least unsure, call dssusy instead.
c  uses sconst, sfesct, chasct, higsct,
c    neusct, vertx, hwidths.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'

      integer i,unphys,hwarning
      real*8 aux,mscale,dsralph3,dsrmq

c-----------------------------------------------------------------------

      if (tanbe.eq.0.0d0) then
         write (*,*)
     &'dsmodelsetup: susy parameters are not properly initialized'
         write (*,'(a,a)')
     &'        pls set tanbe,ma,mu,m2,m1,m3,mass2u,mass2q,mass2d,',
     &        'mass2l,mass2e,asoftu,asoftd,asofte'
         stop
      endif

      unphys=0
      mass(0)=1.d10

c------------------------------------------------------ global constants
      call dssuconst

c------------------------------------------ particle spectrum and mixing
      call dsspectrum(unphys,hwarning)
      if (unphys.ne.0.or.hwarning.ne.0) return
      
c------------- reset running things

      call dssuconst_yukawa_running

c-------------------------------------------------- some useful vertices

      call dsvertx

c--------------------------------------------------- and particle widths

c...Add Higgs widths and QCD correction to them and to vertices

      call dshigwid

c...Add sparticle widths

      call dsspwid

      end



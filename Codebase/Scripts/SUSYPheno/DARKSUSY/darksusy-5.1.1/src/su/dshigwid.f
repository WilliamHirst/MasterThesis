      subroutine dshigwid
c_______________________________________________________________________
c   Wrapper routine to choose which Higgs width routines to use
c   Author: Joakim Edsjo, edsjo@fysik.su.se
c   Date: 2009-03-13
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsidtag.h'

c...Calculate with standard formulae from literature
c...Do this if higwid=1 or if the model is an mSUGRA model
      if (higwid.eq.1.or.modeltype.eq.3) then
         call dshigwid1

c...Calculate with FeynHiggs. Note: only applicable if higloop=5.
c...Note that the actual widhts are extracted in dsfeynhiggs, so
c...we don't need to do anything here.
      elseif (higwid.eq.5) then
         if (higloop.ne.5) then
            write(*,*) 'DS ERROR in dshigwid:'
            write(*,*) 'higwid=5 requires higloop=5, but it is set to ',
     &        higloop
            write(*,*) 'Stopping...'
            stop
         endif

c...Invalid higwid option         
      else
         write(*,*) 'DS ERROR in dshigwid:'
         write(*,*) 'Invalid higwid=',higwid
         write(*,*) 'Stopping...'
         stop         
      endif

      return
      end

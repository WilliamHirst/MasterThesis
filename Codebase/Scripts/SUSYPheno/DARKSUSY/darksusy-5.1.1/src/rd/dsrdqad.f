      subroutine dsrdqad(wrate,mgev,oh2,ierr)
c_______________________________________________________________________
c  present density in units of the critical density times the
c    hubble constant squared. quick and dirty method
c  input:
c    wrate - invariant annihilation rate (real, external)
c    mgev - relic and coannihilating mass in gev
c  output:
c    oh2 - relic density parameter times h**2 (real)
c    ierr - error code (integer)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdtab, dsrdeqn.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996
c  modified: joakim edsjo (edsjo@physto.se) 97-05-12
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      integer ierr
      real*8 wrate,mgev,oh2,xfreeze,yend,w0,p1,w1,a,b,c,tfreeze
      external wrate
      integer k
c-----------------------------------------------------------------------
      k=0
      if (rdinit.ne.1234) then
        call dswrite(0,0,
     &    'dsrdqad: warning: no previous call to dsrdinit')
        call dsrdinit
      endif
      rderr=0
      oh2=0.0d0
c---------------------------------------------------------- get a+b*v^2
      w0=wrate(0.d0)/4.0d0      ! je change 97-05-12
      p1=0.25d0*mgev
      w1=wrate(p1)/4.0d0        ! je change 97-05-12
      a=w0/mgev**2
      b=max((w1-w0)/p1**2,0.0d0)
c------------------------------------------------------ find freeze-out
      c=0.5d0
      xfreeze=20.0d0
c------------------------------------------------ locate tfreeze in dof
      tfreeze=mgev/xfreeze
      khi=nf
      klo=1
  100 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (tgev(k).lt.tfreeze) then
          khi=k
        else
          klo=k
        endif
        goto 100
      endif
c-------------------------------------------------------- recurrence
      xfreeze=log(0.0764*1.22e19*(a+6.d0*b/xfreeze)*c*(c+2.d0)*mgev/
     &     fg(k)/sqrt(xfreeze))
c------------------------------------------------ locate tfreeze in dof
      tfreeze=mgev/xfreeze
      khi=nf
      klo=1
  110 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (tgev(k).lt.tfreeze) then
          khi=k
        else
          klo=k
        endif
        goto 110
      endif
      xfreeze=log(0.0764*1.22d19*(a+6.d0*b/xfreeze)*c*(c+2.d0)*mgev/
     &     fg(k)/sqrt(xfreeze))
c------------------------------------------------ locate tfreeze in dof
      tfreeze=mgev/xfreeze
      khi=nf
      klo=1
  120 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (tgev(k).lt.tfreeze) then
          khi=k
        else
          klo=k
        endif
        goto 120
      endif
c--------------------------------------------- number-to-entropy ratio
      yend=1.d0/(0.264d0*fg(k)*1.22d19*mgev*
     &     (a+3.d0*b/xfreeze)/xfreeze)
c----------------------------------------- normalize to critical density
      oh2=0.72240d8*fh(nf)*mgev*yend
      if (.not.(oh2.ge.0.0d0)) then
        write(*,*) 'error in dsrdqad.f: oh2 = nan'
        write(*,*) '  for model ',rdtag
        write(*,*) 'interesting values follow:'
        write(*,*) '  w0 = ',w0
        write(*,*) '  w1 = ',w1
        write(*,*) '  a = ',a
        write(*,*) '  b = ',b
        write(*,*) '  xfreeze = ',xfreeze
        write(*,*) '  yend = ',yend
        write(*,*)
      endif
 999  ierr=rderr
      end

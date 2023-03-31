      subroutine dsrdens(wrate,npart,kpart,mgev,dof,nrs,rm,rw,
     &  nt,tm,oh2,tf,ierr,iwar)
c_______________________________________________________________________
c  present density in units of the critical density times the
c    hubble constant squared.
c  input:
c    wrate - invariant annihilation rate (real, external)
c    npart - number of particles coannihilating
c    kpart - particle codes (integer)
c    mgev  - relic and coannihilating mass in gev
c    dof   - internal degrees of freedom of the particles
c    nrs   - number of resonances to take special care of
c    rm    - mass of resonances in gev
c    rw    - width of resonances in gev
c    nt    - number of thresholds to take special care of
c            do not include coannihilation thresholds (that's automatic)
c    tm    - sqrt(s) of the thresholds in gev
c  output:
c    oh2   - relic density parameter times h**2 (real*8)
c    tf    - freeze-out temperature in gev (real*8)
c    ierr  - error code (integer)
c      bit 0 (1) = array capacity exceeded. increase nrmax in dsrdcom.h
c          1 (2) = a zero vector is given to dsrdnormlz.
c          2 (4) = step size underflow in dsrdeqn
c          3 (8) = stepsize smaller than minimum hmin in dsrdeqn
c          4 (16) = too many steps in dsrdeqn
c          5 (32) = step size underflow in dsrdqrkck
c          6 (64) = step size smaller than miminum in dsrdqrkck
c          7 (128) = too many steps in dsrdqrkck
c          8 (256) = gpindp integration failed in dsrdthav
c          9 (512) = threshold array too small. increase tharsi in dsrdcom.h
c    iwar  - warning code (integer)
c      bit 0 (1) = a difference of >5waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          1 (2) = a difference of >10waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          2 (4) = a difference of >15waccd in the ratio of w_spline
c                  and w_linear is obtained due to delta_p<dpmin.
c          3 (8) = wimp too heavy, d.o.f. table needs to be
c                  extended to higher temperatures. now the solution
c                  is started at a higher x than xinit (=2).
c          4 (16) = spline interpolated value too high (overflow) during
c                   check of interpolation accuracty (dsrdwintpch)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdtab, dsrdeqn.
c  authors: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996 and
c           joakim edsjo (edsjo@physto.se) 30-april-98
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      integer npart,nrs,nt,nfcn,ierr,iwar
      real*8 wrate,dsrdwintp,mgev(npart),dof(npart),oh2,
     &   xstart,xend,yend,rm(nrs),rw(nrs),tm(nt),xf,tf
      integer kpart(npart)
      external wrate,dsrdwintp
      integer k
      real*8 tstart
c      external dsrdcom ! set up common block variables
c-----------------------------------------------------------------------

      call dsrdcom ! set up common block variables

      rderr=0
      rdwar=0

      if (rdinit.ne.1234) then
        call dswrite(0,0,
     &    'dsrdens: warning: no previous call to dsrdinit')
        call dsrdinit
      endif

c      write(*,*) 'dsrdens: npart=',npart
      
      call dsrdstart(npart,kpart,mgev,dof,nrs,rm,rw,nt,tm)

      if (rderr.ne.0) goto 999

c-------------------------------------------------------------- initial x

      xstart=max(xinit,1.0001d0*mco(1)/tgev(1))  ! je corr 97-03-14
      if (xstart.ne.xinit) rdwar=ibset(rdwar,3)

c--------------------------------------------------- locate tstart in dof

      tstart=mgev(1)/xstart
      khi=nf
      klo=1
  100 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (tgev(k).lt.tstart) then
          khi=k
        else
          klo=k
        endif
        goto 100
      endif

c------------------------------------------- tabulate the invariant rate

      call dsrdtab(wrate,xstart)
      if (rderr.ne.0) goto 999

c------------------------------------------ determine integration limits

      call dsrdthlim

c---------------------------------------------------------- call dsrdeqn

      call dsrdeqn(dsrdwintp,xstart,xend,yend,xf,nfcn)
      tf=mco(1)/xf

c----------------------------------------- normalize to critical density

 999  if (rderr.eq.0) then
        oh2=0.70365d8*fh(nf)*mco(1)*yend  ! je 970404, t_0=2.726 k
      else
        oh2=0.0d0     ! something went wrong, changed to 0.0 020816, je
        tf=-1.0d0
      endif
      ierr=rderr
      iwar=rdwar
c      write (*,*) 'PG-DEBUG dsrdens: oh2,yend,mco(1),fh(nf)',
c     &     oh2,yend,mco(1),fh(nf)

      end

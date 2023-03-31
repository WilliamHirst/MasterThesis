************************************************************************
                      subroutine dsrdwfunc(x,wrate)
c_______________________________________________________________________
c  write out dsrdfunc for the given x = mass/temperature
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdfunc
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-01-20
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 x,du,w,wrate,dsrdfunc,utop,xmin,b,u
      external wrate,dsrdfunc
      integer uout,i,nmax
      character*25 filefun
      character*12 idtagold
      data idtagold /'000000000000'/
      save idtagold
      parameter (uout=14)
c------------------------------------------------------------------------
c...check if called before for the same model in which case dsrdfunc will
c...produce only zeros (since rderr.ne.0)
      if (idtagold.eq.rdtag) then
        return
      else
        idtagold=rdtag
      endif

      filefun='fun'//rdtag//'.dat'
      open (unit=uout,file=filefun,status='unknown',form='formatted')
c      write (uout,110) x
 110  format('#...u...   .dsrdfunc.   x = ',1x,g14.8)

      utop=umax*0.9999d0    ! must be <= utop in dsrdtab.f
      xmin=2.0d0              ! must be >= xstart in dsrdens.f
      b=min(sqrt(sqrt(4.0d0*x**2+
     &  utop**2*(4.0d0*x**2/xmin+x**2/xmin**2*utop**2))
     &   -2.0d0*x),15.0d0)

      nmax=5000
      du=b/dble(nmax)
      do i=1,nmax
        u=dble(i)*du
        w=dsrdfunc(u,x,wrate)
        write(uout,'(2(x,e14.8))') u,w*1.0d10
      enddo
      close(uout)

 1001 format(i1)
 1002 format(i2)
      end









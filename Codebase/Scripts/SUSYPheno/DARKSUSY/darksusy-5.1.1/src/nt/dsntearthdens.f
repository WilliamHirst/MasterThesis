c      program test
c      implicit none
c      integer i
c      real*8 radius,dsntearthdens,dsntearthmassint,dsntearthmass,
c     &  dsntearthpotint,dsntearthpot,tmp,dsntearthvesc,
c     &  dsntearthdenscomp
c      real*8 depth(42)
c      data depth/0.0d0,3.0d0,15.0d0,24.0d0,80.0d0,
c     &  219.99d0,220.0d0,399.99d0,400.0d0,500.0d0,
c     &  600.0d0,669.99d0,670.0d0,770.0d0,1000.0d0,
c     &  1250.0d0,1500.0d0,1750.0d0,2000.0d0,2250.0d0,
c     &  2500.0d0,2750.0d0,2899.99d0,2900.0d0,3000.0d0,
c     &  3250.0d0,3500.0d0,3750.0d0,4000.0d0,4250.0d0,
c     &  4500.0d0,4750.0d0,5000.0d0,5149.99d0,5150.0d0,
c     &  5250.0d0,5500.0d0,5750.0d0,6000.0d0,6250.0d0,
c     &  6371.0d0,6379.0d0/   ! in km
c
cc      do i=42,1,-1
cc        radius=max((6378.140-depth(i))*1.0d3,0.0d0)
cc        write(*,*) radius,dsntearthpotint(radius)
cc      enddo
c
c      do i=0,1100
c        radius=dble(i)/dble(1000.0d0)*6378.14d0*1.0d3
c        write(*,'(10(x,e12.6))') radius,dsntearthdens(radius),
c     &    dsntearthmassint(radius)/1.0d24,
c     &    dsntearthmass(radius)/1.0d24,
c     &    dsntearthpotint(radius),dsntearthpot(radius),
c     &    dsntearthvesc(radius)
c
c        write(*,'(10(x,e12.6))') radius,
c     &    dsntearthdenscomp(radius,16),
c     &    dsntearthdenscomp(radius,24),
c     &    dsntearthdenscomp(radius,28),
c     &    dsntearthdenscomp(radius,56)
c      enddo
c
c      end



***********************************************************************
*** dsntearthdens gives the density in the earth as a function of radius
*** the radius should be given in m and the density is returned in
*** g/cm^3
*** author: joakim edsjo
*** date: march 19, 1999
***********************************************************************

      real*8 function dsntearthdens(r)
      implicit none
      include 'dsearth.h'

      real*8 r,dd,dpl
      integer i,j

      data eadepth/0.0d0,3.0d0,15.0d0,24.0d0,80.0d0,
     &  219.99d0,220.0d0,399.99d0,400.0d0,500.0d0,
     &  600.0d0,669.99d0,670.0d0,770.0d0,1000.0d0,
     &  1250.0d0,1500.0d0,1750.0d0,2000.0d0,2250.0d0,
     &  2500.0d0,2750.0d0,2899.99d0,2900.0d0,3000.0d0,
     &  3250.0d0,3500.0d0,3750.0d0,4000.0d0,4250.0d0,
     &  4500.0d0,4750.0d0,5000.0d0,5149.99d0,5150.0d0,
     &  5250.0d0,5500.0d0,5750.0d0,6000.0d0,6250.0d0,
     &  6371.0d0,6379.0d0/   ! in km
      data eadens/2.60d0,2.60d0,2.90d0,3.38d0,3.37d0,
     &  3.36d0,3.44d0,3.54d0,3.72d0,3.85d0,
     &  3.98d0,3.99d0,4.38d0,4.44d0,4.58d0,
     &  4.72d0,4.86d0,4.99d0,5.12d0,5.25d0,
     &  5.37d0,5.50d0,5.57d0,9.90d0,10.08d0,
     &  10.44d0,10.77d0,11.06d0,11.32d0,11.55d0,
     &  11.75d0,11.93d0,12.09d0,12.17d0,12.76d0,
     &  12.82d0,12.92d0,13.00d0,13.06d0,13.09d0,
     &  13.09d0,13.09d0/


      dd=(r_earth-r)*1.0d-3
      if (dd.gt.r_earth*1.d-3.or.dd.lt.0.0d0) then
        dsntearthdens=0.0d0
        return
      endif

      call dshunt(eadepth,42,dd,j)
      if (j.le.42.and.j.ge.1) goto 20

      dsntearthdens=0.0d0
      return

 20   dpl=(dd-eadepth(j))/(eadepth(j+1)-eadepth(j))

      dsntearthdens=eadens(j)*(1.0d0-dpl)+eadens(j+1)*dpl

      return

      end

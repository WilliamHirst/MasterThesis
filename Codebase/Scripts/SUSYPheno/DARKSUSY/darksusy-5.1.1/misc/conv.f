      program conv
      implicit none
      real*8 u,f

      open(unit=13,file='dat/gauss.txt',
     &  status='old',form='formatted')
      open(unit=14,file='dat/vdfearth-sdgauss.dat',
     &  status='unknown',form='formatted')

 10   read(13,*,end=100) u,f
      u=u/1000.0d0  ! m/s -> km/s
      f=f*1000.0d0  ! s/m -> s/km
      if (u.gt.0.0d0) then 
        f=f/u
      else
        f=0.0d0
      endif

      write(14,'(2(1x,E14.8))') u,f
      goto 10

 100  close(13)
      close(14)

      end

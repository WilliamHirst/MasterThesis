**********************************************************************
*** subroutine dsgalproptable
*** input: force - if file found, compute anyway
*** 
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine dsgalprop_gettable(force)
      implicit none

      include 'dsgalpropcom.h'

      integer m,n,i,j,force
      real *8 keinput,kemin,kemax,kefactor,y,dummy(5)

      gpgffile=dsinstall
      call dscharadd(gpgffile,'share/DarkSUSY/FITS/GF_50p_')
      call dscharadd(gpgffile,gpmodtag)
      call dscharadd(gpgffile,'.dat')
      gpgfunit=97
      if (force.ne.0) then
         write (6,*) 'New Green''s functions requested'
         write (6,*) 'Calling GALPROP to construct Green''s functions'
         goto 35
      endif

      open(unit=gpgfunit,file=gpgffile,status='unknown',
     +     form='formatted',err=10)
      do i=1,gpnumin
         do j=1,gpnumout
            read(gpgfunit,*,end=20) kegpgf(i,j),kepgpgf(i,j),
     +           epgpgf(i,j),pbgpgf(i,j)
         enddo
      enddo
      read(gpgfunit,*,end=5)
     +     dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)

      close(gpgfunit)
      write(6,*) 'Warning: old GALPROP Green''s function file ',
     +     'was longer than needed'
      gpgfread=.true.
      return

 5    close(gpgfunit)
      gpgfread=.true.
      return

 10   close(gpgfunit)
      write (6,*) 'Calling GALPROP to construct Green''s functions'
      goto 30

 20   close(gpgfunit)
      write (6,*) 'Green''s function table incomplete, calling GALPROP'

 30   write (6,*) 'This may take a while'
 35   call dsgalprop_maketable()

      open(unit=gpgfunit,file=gpgffile,status='unknown',
     +     form='formatted',err=40)
      do i=1,gpnumin
         do j=1,gpnumout
            read(gpgfunit,*,end=40) kegpgf(i,j),kepgpgf(i,j),
     +           epgpgf(i,j),pbgpgf(i,j)
         enddo
      enddo
      read(gpgfunit,*,end=45)
     +     dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)
      close(gpgfunit)

      write (6,*) 'Warning: newly computed GALPROP Green''s function ',
     +     'file was longer than needed'
      gpgfread=.true.
      return
      
 45   close(gpgfunit)
      write (6,*) 'GALPROP Green''s functions completed'
      gpgfread=.true.
      return

 40   write (6,*) 'Error reading GALPROP Green''s function table'
      stop
      end

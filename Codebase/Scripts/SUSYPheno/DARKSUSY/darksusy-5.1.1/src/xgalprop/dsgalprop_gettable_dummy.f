**********************************************************************
*** subroutine dsgalprop_gettable_dummy
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

 20   close(gpgfunit)
      write(6,*) 'GALPROP file incomplete, exiting'
      write (6,*) 'To compute Green''s functions, compile GALPROP'
      stop

 10   close(gpgfunit)
      write (6,*) 'Can''t open Green''s function file, exiting'
      write (6,*) 'To compute Green''s functions, compile GALPROP'
      stop
      end

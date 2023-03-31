      subroutine dsgalprop_write_galdef(ekin,dekin)
      implicit none
      integer i, flag
      real*8 ekin,dekin

      include 'dsgalpropcom.h'

      open (unit=66,file=gpgaldeffile,status='unknown',
     +     form='formatted',err=101)
      open (unit=67,file=gpgaldeftmpfile,status='unknown',
     +     form='formatted')

c...  copy the two line header
      do i=1,2
         call dsgalprop_copy_line(flag)
         if (flag.eq.0) then
            write (6,*) 'error reading galdef file, stopping'
            stop
         endif
      enddo

      write (67,'(a22,1x,a1)') 'max_Z                =', '1'
      write (67,'(a22,1x,a1)') 'DM_positrons         =', '1'
      write (67,'(a22,1x,a1)') 'DM_electrons         =', '0'
      write (67,'(a22,1x,a1)') 'DM_antiprotons       =', '1'
      write (67,'(a22,1x,a1)') 'DM_gammas            =', '0'
      write (67,'(a22,1x,e14.8)') 'DM_double3           =', ekin
      write (67,'(a22,1x,e14.8)') 'Ekin_factor          =', dekin

c...  append the contents of the initial galdef file
 100  call dsgalprop_copy_line(flag)
      if (flag.eq.1) goto 100

      close(66)
      close(67)
      return

 101  write (6,*) 'Cannot open file ', gpgaldeffile
      end

      subroutine dsgalprop_copy_line(flag)
      implicit none
      character*1000 line
      character*80 formatstring
      integer i, numchars, flag
      
      do i=1,1000
         line(i:i)='\000'
      enddo
      read (66,'(a)',end=300) line
      do i=1000,1,-1
         if (.not.(line(i:i).eq.'')) then
            numchars=i
            goto 200
         endif
      enddo
      numchars=0
 200  write (formatstring,1000) numchars
      write (67,formatstring) line
      flag=1
      return
 300  flag=0
      return
 1000 format('(a',I4,')')
      end

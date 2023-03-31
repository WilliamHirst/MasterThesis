**********************************************************************
*** subroutine dsgalproptable
*** input: n number of points at which to compute the Green function
*** 
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine dsgalprop_maketable()
      implicit none
      character*200 header,fitsfile,strnum
      integer i,j
      real*8 ekin
      
      include 'dsgalpropcom.h'

      gpgaldeffile=dsinstall
      call dscharadd(gpgaldeffile,
     +     'share/DarkSUSY/GALDEF/galdef_50p_')
      call dscharadd(gpgaldeffile,gpmodtag)
      call dscharadd(gpgaldeffile,'_template')
      gpgfunit=97
      open(unit=gpgfunit,file=gpgffile,form='Formatted')
      do i=1,gpnumin
         ekin=10.d0**((i-1)/dble(gpnumdecade)+2.d0)
         write(strnum,'(i4.4)') i
         header=gpmodtag
         call dscharadd(header,strnum)
         gpgaldeftmpfile=dsinstall
         call dscharadd(gpgaldeftmpfile,
     +        'share/DarkSUSY/GALDEF/galdef_50p_')
         call dscharadd(gpgaldeftmpfile,header)
         fitsfile=dsinstall
         call dscharadd(fitsfile,
     +        'share/DarkSUSY/FITS/nuclei_50p_')
         call dscharadd(fitsfile,header)
         call dsgalprop_write_galdef(ekin,10.d0**(1.d0/gpnumdecade))
         call dscharadd(header,'\000')
 100     call dsgalprop_run_galprop(header)
         call dsgalprop_read_fits(fitsfile,ekin)
      enddo

      close(gpgfunit)
      end

**********************************************************************
*** subroutine dsgalproptable
*** input: n number of points at which to compute the Green function
*** 
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine dsgalprop_maketable_one(i)
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

      ekin=10.d0**((i-1)/dble(gpnumdecade)+2.d0)
      write(strnum,'(i4.4)') i

      gpgffile=dsinstall
      call dscharadd(gpgffile,'share/DarkSUSY/FITS/GF_50p_')
      call dscharadd(gpgffile,gpmodtag)
      call dscharadd(gpgffile,strnum)
      call dscharadd(gpgffile,'.dat')
      gpgfunit=97
      open(unit=gpgfunit,file=gpgffile,form='Formatted')
      
      header=gpmodtag
      call dscharadd(header,strnum)
      gpgaldeftmpfile=dsinstall
      call dscharadd(gpgaldeftmpfile,
     +     'share/DarkSUSY/GALDEF/galdef_50p_')
      call dscharadd(gpgaldeftmpfile,header)
      fitsfile=dsinstall
      call dscharadd(fitsfile,
     +     'share/DarkSUSY/FITS/nuclei_50p_')
      call dscharadd(fitsfile,header)
      call dsgalprop_write_galdef(ekin,10.d0**(1.d0/gpnumdecade))
      call dscharadd(header,'\000')
 100  call dsgalprop_run_galprop(header)
      call dsgalprop_read_fits(fitsfile,ekin)

      close(gpgfunit)
      end

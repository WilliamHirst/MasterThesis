      subroutine dsreadnuclides
      implicit none
      include 'dsdirver.h'
      include 'dsnuclides.h'
      character*200 filename,scratch*2
      character*300 msg
      integer dsi_trim,itmp,i
      real*8 atomicmassunit,tmp,etmp,mpu,mnu
      parameter (atomicmassunit=0.931494028d0)
      parameter (mpu=1.00727646688d0)
      parameter (mnu=1.0086649156d0)
      filename=dsinstall(:dsi_trim(dsinstall))//
     &     'share/DarkSUSY/'//'nuclides.dat'
      write (msg,'(a,a)') 'dsreadnuclides: reading ',filename
      call dswrite(0,0,msg)
      open (unit=13,file=filename,status='unknown',
     &     form='formatted')
      read(13,'(a)',err=2000) scratch
      do i=1,nnucld
         read(13,*,err=2000) nucldsym(i),nucldz(i),nuclda(i),nucldn(i),
     &        tmp,itmp,nucldnatab(i),nucldj(i),nucldmu(i),etmp
         nucldam(i)=tmp*atomicmassunit
         nucldm(i)=(nucldz(i)*mpu+nucldn(i)*mnu)*atomicmassunit
         if (etmp.ne.9999.d0) then
            nucldm(i)=nucldm(i)-etmp*nuclda(i)*1.d-3
         endif
         if (itmp.eq.0) nucldstable(i)=.true.
         if (itmp.eq.1) nucldstable(i)=.false.
      enddo
      close(13)
      inucld=1
ccc begin dbg
ccc      do i=1,nnucld
ccc         write(*,*) nucldsym(i),nucldz(i),nuclda(i),nucldn(i),nucldm(i),
ccc     &        nucldstable(i),nucldnatab(i),nucldj(i),nucldmu(i)
ccc      enddo
ccc end dbg
      return
 2000 continue
      close(13)
      write (*,*) 'dsreadnuclides: error while reading ',filename
      write (*,*) 'DarkSUSY will stop'
      stop
      end

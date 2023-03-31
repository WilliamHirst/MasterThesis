      subroutine dsntctabread(wh,i,file)

***********************************************************************
*** Reads in tabulated capture rates (apart from cross section)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to read
***        file = file name to read
*** Author: Joakim Edsjo
*** Date: 2003-11-27
*** Modified: 2004-02-01
***********************************************************************

      implicit none
      include 'dsntcap.h'
      include 'dshmcom.h'

      integer i,fl,l,m,ifile

      character*2 wh,num,whfile
      character*12 vdf,vdfc

      character*200 file,scr
      logical bad

      real*8 tmp1,tmp2,tmp3

c...Now read in the data file

c...Make sure 'num' and 'numc' is treated consistenly     
      vdfc=veldf
      if (vdfc.eq.'numc') vdfc='num'

      write(*,*) 'dsntctabread: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='old',err=200)


      read(13,101) whfile,vdf
 101  format(1x,1x,A2,1x,A12,1x,A)


c...Check consistency of files
      bad=.false.
      if (wh.eq.'su'.or.wh.eq.'SU') then
        if (whfile.ne.'su'.and.whfile.ne.'SU') bad=.true.
        if (vdf.ne.vdfc) bad=.true.
        if (bad) then
          write(*,*) 'ERROR in dsntctabread: ',
     &      'File type mismatch.'
          write(*,*) 'Tried to read wh=',wh,' and vdfc=',
     &      vdfc,'.'
          write(*,*) 'but found wh=',whfile,' and vdfc=',vdfc,
     &      'in file. Stopping.'
          stop
        endif
      else  ! if (wh.eq.'ea'.or.wh.eq.'EA') then
        if (whfile.ne.'ea'.and.whfile.ne.'EA') bad=.true.
        if (vdf.ne.veldfearth) bad=.true.
        if (bad) then
          write(*,*) 'ERROR in dsntctabread: ',
     &      'File type mismatch.'
          write(*,*) 'Tried to read wh=',wh,' and veldfearth=',
     &      veldfearth,'.'
          write(*,*) 'but found wh=',whfile,' and veldfearth=',vdf,
     &      'in file. Stopping.'
          stop
        endif
      endif

      do l=0,nc
        if (wh.eq.'su'.or.wh.eq.'SU') then
          read(13,*) ctabsusi(l,i),ctabsusd(l,i)
        else
          read(13,*) ctabea(l,i)
        endif
      enddo

      close(13)
      
      return

 200  continue  ! we get here if file didn't exist

      close(13)

c...Create file
      write(*,*) 'WARNING in dsntctabread:'
      write(*,*) 'The requested file ',file,' does not exist.'
      write(*,*) 'I will create it for you',
     &  ' (this only needs to be done once),'
      write(*,*)  'but it could take several hours.',
     &  ' Be patient please...'
      call dsntctabcreate(wh,i)
      call dsntctabwrite(wh,i,file)

      return

      end


      


      subroutine dsntctabwrite(wh,i,file)

***********************************************************************
*** Writes out tabulated capture rates (apart from cross section)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to write
***        file = file to write to
*** Author: Joakim Edsjo
*** Date: 2003-11-27
*** Modified: 2004-02-01
***********************************************************************

      implicit none
      include 'dsntcap.h'
      include 'dshmcom.h'
      include 'dsdirver.h'

      integer i,fl,l,m

      character*2 wh,num

      character*200 file,scr


c...Now write the data file
     
      write(*,*) 'dsntctabwrite: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='unknown')

      if (wh.eq.'su'.or.wh.eq.'SU') then
        write(13,101) wh,veldf,dsversion
      else
        write(13,101) wh,veldfearth,dsversion
      endif
 101  format('#',1x,A2,1x,A12,1x,A)

      do l=0,nc
        if (wh.eq.'su'.or.wh.eq.'SU') then
          write(13,'(2(1x,E14.8))') ctabsusi(l,i),ctabsusd(l,i)
        else
          write(13,'(1x,E14.8)') ctabea(l,i)
        endif
      enddo
      
      close(13)

      return

      end


      


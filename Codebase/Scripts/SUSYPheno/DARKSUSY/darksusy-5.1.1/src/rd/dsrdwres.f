************************************************************************
                      subroutine dsrdwres
c_______________________________________________________________________
c  write out dsrdtab and check the interpolation routine
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdtab,dsrdwintp.f
c  author: joakim edsjo (edsjo@physto.se)
c  date: 96-03-26
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      external dsrdwintp
      real*8 wpmax,dp,w,p,dsrdwintp,pmin
      integer uout,i,nmax
      character*25 filetab,fileint
c      character*1 filenr1
c      character*2 filenr2
      integer nmb
      data nmb /0/
      save nmb
      parameter (uout=14)
c------------------------------------------------------------------------
      nmb=nmb+1
c      if (nmb.lt.10) then
c        write(filenr1,1001) nmb
c        filetab='tab'//filenr1//'.dat'
c        fileint='int'//filenr1//'.dat'
c      elseif (nmb.ge.10.and.nmb.lt.100) then
c        write(filenr2,1002) nmb
c        filetab='tab'//filenr2//'.dat'
c        fileint='int'//filenr2//'.dat'
c      else
c        return
c      endif

      filetab='tab'//rdtag//'.dat'
      fileint='int'//rdtag//'.dat'
      open (unit=uout,file=filetab,status='unknown',form='formatted')
c      write (uout,*) '#...p...   ...w...'
      do i=1,nr
        write(uout,*) pp(indx(i)),exp(yy(indx(i)))*1.0d10,
     &    yy2(indx(i))*1.0d10
      enddo
      close(uout)

      open (unit=uout,file=fileint,status='unknown',form='formatted')
c      write (uout,*) '#...p...   ...w...'
      nmax=nr*50
      nmax=5000
      pmin=0.0d0
      wpmax=pp(indx(nr))*0.9999d0
c      pmin=31.0d0
c      wpmax=33.0d0
      dp=(wpmax-pmin)/dble(nmax)
      do i=1,nmax
        p=dble(i)*dp+pmin
        w=dsrdwintp(p)
        write(uout,*) p,w*1.0d10
      enddo
      close(uout)

 1001 format(i1)
 1002 format(i2)
      end









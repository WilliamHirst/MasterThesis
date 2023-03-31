      subroutine dsntctabcreate(wh,i)

***********************************************************************
*** Creates tabulated capture rates (apart from cross section)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to store the results in
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      implicit none

      include 'dsntcap.h'
      include 'dshmcom.h'
      real*8 dsntcapearthnum,dsntcapsunnum,
     &  sigsi,sigsd,mx,rhotmp,tmp1,tmp2
      integer i,j,n
     
      character*2 wh

      sigsi=0.0d0
      sigsd=0.0d0
      rhotmp=rhox
      rhox=0.3d0
      do j=0,nc
        mx=10**(dble(j)*5.0d0/dble(nc))
        write(*,*) '  Taking care of mass number ',j,'/',nc,
     &    ' (',mx,' GeV)'
        
        if (wh.eq.'su'.or.wh.eq.'SU') then
          sigsi=1.0d-40
          sigsd=0.0d0
          tmp1=dsntcapsunnum(mx,sigsi,sigsd)
          ctabsusi(j,i)=tmp1
          sigsi=0.0d0
          sigsd=1.0d-40
          tmp2=dsntcapsunnum(mx,sigsi,sigsd)
          ctabsusd(j,i)=tmp2
          write(*,*) '    result: ',tmp1,tmp2
        else
          sigsi=1.0d-40
          tmp1=dsntcapearthnum(mx,sigsi)
          ctabea(j,i)=tmp1
          write(*,*) '    result: ',tmp1
        endif
      enddo

      rhox=rhotmp
      return

      end


      


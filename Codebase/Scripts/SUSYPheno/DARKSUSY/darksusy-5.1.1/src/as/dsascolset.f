      subroutine dsascolset(iloctype)
      implicit none
      include 'dsascom.h'
      integer iloctype
      integer j1,j2
c
      if(iloctype.eq.1) then
        do j1=1,3
        do j2=1,3   
           colfactor(j1,j2)=0.d0
        enddo
        enddo
        do j1=4,6
        do j2=1,3   
           colfactor(j1,j2)=0.d0
           colfactor(j2,j1)=0.d0
        enddo
        enddo
        do j1=3,6
        do j2=3,6   
           colfactor(j1,j2)=1.d0
        enddo
        enddo
      elseif(iloctype.eq.2) then
        colfactor(1,1)=2.d0
        colfactor(1,2)=-2.d0/3.d0
        colfactor(1,3)=0.d0
        colfactor(1,4)=0.d0
        colfactor(1,5)=4.d0
        colfactor(1,6)=0.d0
        colfactor(2,1)=-2.d0/3.d0
        colfactor(2,2)=2.d0
        colfactor(2,3)=-2.d0/3.d0
        colfactor(2,4)=4.d0
        colfactor(2,5)=0.d0
        colfactor(2,6)=4.d0
        colfactor(3,1)=0.d0
        colfactor(3,2)=-2.d0/3.d0
        colfactor(3,3)=2.d0
        colfactor(3,4)=0.d0
        colfactor(3,5)=4.d0
        colfactor(3,6)=0.d0
        colfactor(4,1)=0.d0
        colfactor(4,2)=4.d0
        colfactor(4,3)=0.d0
        colfactor(4,4)=9.d0
        colfactor(4,5)=3.d0
        colfactor(4,6)=0.d0
        colfactor(5,1)=4.d0
        colfactor(5,2)=0.d0
        colfactor(5,3)=4.d0
        colfactor(5,4)=3.d0
        colfactor(5,5)=9.d0
        colfactor(5,6)=3.d0
        colfactor(6,1)=0.d0
        colfactor(6,2)=4.d0
        colfactor(6,3)=0.d0
        colfactor(6,4)=0.d0
        colfactor(6,5)=3.d0
        colfactor(6,6)=9.d0
      else
        write(*,*) 'DS: dsascolset called with wrong type = ',iloctype
        write(*,*) 'DS: rather than 1=leptons or 2=quarks'
        write(*,*) 'DS: program stopped'
        stop
      endif  
      return
      end

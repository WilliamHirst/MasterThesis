************************************************************************
*** dsrdstart stores coannihilation, resonance and threshold information
*** in common blocks and sort them
*** author: joakim edsjo (edsjo@physto.se)
*** date: 03-march-98
*** modifed: 08-may-98
***   08-may-98: bug with mdof not being sorted correctly fixed.
***   27-feb-02: bug with allocation of threshold array fixed.
***              increased number of possible coannihilations
************************************************************************


      subroutine dsrdstart(npart,kpart,mgev,dof,nrs,rm,rw,nt,tm)
      implicit none

      include 'dsrdcom.h'
      integer npart,nrs,nt,i,j
      real*8 mgev(npart),dof(npart),rm(nrs),rw(nrs),tm(nt),tmp
      integer kpart(npart),ktmp
      real*8 tmall(tharsi)

c-----------------------------------------------------------------------

c...copy to common blocks
      nco=npart
      do i=1,npart
         mco(i)=abs(mgev(i))
         mdof(i)=dof(i)
         kcoann(i)=kpart(i)
      enddo

      nres=nrs
      do i=1,nrs
        rgev(i)=rm(i)
        rwid(i)=rw(i)
      enddo
      
      nth=nt
      do i=1,nth
        tmall(i)=tm(i)
      enddo

c...When we get here, nth holds the number of thresholds supplied by the User
c...Now add coannihilation-thresholds
      if (npart.gt.1) then
        do i=1,npart
          do j=max(2,i),npart
            nth=nth+1
            if (nth.ge.tharsi) then
              write(rdluerr,*) 'DarkSUSY error in dsrdstart.',
     &              ' Increase tharsi in include/dsrdcom.h to ',
     &              nth,' or larger.'
              rderr=ibset(rderr,9)
              return
            endif
            tmall(nth)=mco(i)+mco(j)
          enddo
        enddo
      endif

c...sort
      if (nco.ge.2) then
        do i=1,nco-1
          do j=nco-1,i,-1
            if (mco(j).gt.mco(j+1)) then
              tmp=mco(j+1)
              mco(j+1)=mco(j)
              mco(j)=tmp
              tmp=mdof(j+1)
              mdof(j+1)=mdof(j)
              mdof(j)=tmp
              ktmp=kcoann(j+1)
              kcoann(j+1)=kcoann(j)
              kcoann(j)=ktmp
            endif
          enddo
        enddo
      endif

      if (nres.ge.2) then
        do i=1,nres-1
          do j=nres-1,i,-1
            if (rgev(j).gt.rgev(j+1)) then
              tmp=rgev(j+1)
              rgev(j+1)=rgev(j)
              rgev(j)=tmp
              tmp=rwid(j+1)
              rwid(j+1)=rwid(j)
              rwid(j)=tmp
            endif
          enddo
        enddo
      endif

      if (nth.ge.2) then
        do i=1,nth-1
          do j=nth-1,i,-1
            if (tmall(j).gt.tmall(j+1)) then
              tmp=tmall(j+1)
              tmall(j+1)=tmall(j)
              tmall(j)=tmp
            endif
          enddo
        enddo
      endif


c...convert thresholds to effective momentum and store in common block
      pth(0)=0.0d0
      incth(0)=1
      do i=1,nth
        pth(i)=sqrt(tmall(i)**2/4.0d0-mco(1)**2)
        incth(i)=1
      enddo

      return

      end


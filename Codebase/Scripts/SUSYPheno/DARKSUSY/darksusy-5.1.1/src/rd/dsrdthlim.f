      subroutine dsrdthlim
c_______________________________________________________________________
c  determine which limits in p_eff (or rather u) to use when
c  integrating for the thermal average
c  author: joakim edsjo (edsjo@physto.se)
c  date: 98-04-30
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 tmp,intlim(0:2*tharsi),p2,pres,dpth
      integer i,j,tmpi,limtype(0:2*tharsi)

c-----------------------------------------------------------------------

c...thresholds
      nlim=0
      intlim(0)=0.0d0
      limtype(0)=0
      do i=1,nth+1
        if (incth(i).eq.1) then
          nlim=nlim+1
          intlim(nlim)=pth(i)
          limtype(nlim)=1
        endif
      enddo

c...resonances
      do i=1,nres
        p2=(rgev(i)-2.0d0*rwid(i))**2/4.0d0-mco(1)**2
        if (p2.gt.0.0d0) then
          pres=sqrt(p2)
          nlim=nlim+1
          intlim(nlim)=pres
          limtype(nlim)=0
        endif
        p2=(rgev(i)+2.0d0*rwid(i))**2/4.0d0-mco(1)**2
        if (p2.gt.0.0d0) then
          pres=sqrt(p2)
          nlim=nlim+1
          intlim(nlim)=pres
          limtype(nlim)=0
        endif
      enddo

c...add one limit at pdivr
      nlim=nlim+1
      intlim(nlim)=pdivr*mco(1)
      limtype(nlim)=0

c...sort
      if (nlim.ge.2) then
        do i=1,nlim-1
          do j=nlim-1,i,-1
            if (intlim(j).gt.intlim(j+1)) then
              tmp=intlim(j+1)
              intlim(j+1)=intlim(j)
              intlim(j)=tmp
              tmpi=limtype(j+1)
              limtype(j+1)=limtype(j)
              limtype(j)=tmpi
            endif
          enddo
        enddo
      endif

c...make limits of it
      dpth=dpthr*mco(1)
      do i=0,nlim-1
        if (limtype(i).eq.0) then
          plow(i+1)=intlim(i)
        else
          plow(i+1)=intlim(i)+dpth
        endif
        if (limtype(i+1).eq.0) then
          phigh(i+1)=intlim(i+1)
        else
          phigh(i+1)=intlim(i+1)-dpth
        endif
      enddo

c...check if some limits are mixed due to close thresholds
 10   do i=1,nlim
        if (plow(i).ge.phigh(i)) then
          do j=i,nlim-1
            plow(j)=plow(j+1)
            phigh(j)=phigh(j+1)
          enddo
          nlim=nlim-1
          goto 10
        endif
      enddo

      return

      end

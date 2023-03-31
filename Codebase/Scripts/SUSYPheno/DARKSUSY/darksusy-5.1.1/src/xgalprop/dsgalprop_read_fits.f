**********************************************************************
*** subroutine dsgalprop_read_fits
*** 
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine dsgalprop_read_fits(filename,ekin)
      implicit none
      character*200 filename,comment
      integer naxis,naxes(10),unit,status,blocksize,bitpix,nelements,
     +     fpixel(10),i,j,k,l,anynul,group,rwmode
      real array(1000,1,1000,10),buffer(1000000),nullval,
     +     gfe(1000),gfp(1000),ek0,dek
      real r0,dr,frac
      real*8 ekin

      include 'dsgalpropcom.h'

      status=0
      blocksize=1
      unit=131
      group=0
      rwmode=0
      nullval=0.0
      
      do i=1,10
         fpixel(i)=1
      enddo

      call ftopen(unit,filename,rwmode,blocksize,status)
      call ftgipr(unit,4,bitpix,naxis,naxes,status)
      call ftgkye(unit,'CRVAL1',r0,comment,status)
      call ftgkye(unit,'CDELT1',dr,comment,status)
      call ftgkye(unit,'CRVAL3',ek0,comment,status)
      call ftgkye(unit,'CDELT3',dek,comment,status)
      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
      call ftgpve(unit,group,fpixel,nelements,nullval,
     +     buffer,anynul,status)
      call ftclos(unit,status)


      do i=1,naxes(1)
         do j=1,naxes(2)
            do k=1,naxes(3)
               do l=1,naxes(4)
                  array(i,j,k,l)=
     +                 buffer(i+
     +                 naxes(1)*(j-1+
     +                 naxes(2)*(k-1+
     +                 naxes(3)*(l-1))))
               enddo
            enddo
         enddo
      enddo
      
      k=int((8.5d0-r0)/dr)
      frac=(8.5d0-r0)/dr-k

      do i=1,naxes(3)
         gfe(i)=(1.-frac)*array(k+1,1,i,1)+frac*array(k+2,1,i,1)
         gfp(i)=(1.-frac)*array(k+1,1,i,2)+frac*array(k+2,1,i,2)
         if (gfe(i).lt.0.d0) gfe(i)=0.d0
         if (gfp(i).lt.0.d0) gfp(i)=0.d0
         write (gpgfunit,'(2(1x,e10.4),2(1x,e14.8))')
     +        ekin,10.d0**((i-1)/dble(gpnumdecade)), gfe(i), gfp(i)
      enddo

      end

      subroutine dsrdnormlz(x,y,nx,ny)
c_______________________________________________________________________
c  find the unit vector (nx,ny) in the same direction as (x,y).
c  input:
c    x,y - coordinates of the vector (real)
c  output:
c    nx,ny - coordinates of the versor (real)
c  common:
c    'dsrdcom.h' - included common blocks
c  called by dsrdtab.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 x,y,absx,absy,nx,ny
c-----------------------------------------------------------------------
      absx=abs(x)
      absy=abs(y)
      if(absx.gt.absy)then
        nx=sign(1.0d0/sqrt(1.0d0+(absy/absx)**2),x)
        ny=(y/x)*nx
      else
        if(absy.eq.0.0d0)then
          write (rdluerr,*) 'error in dsrdtab.dsrdnormlz: zero vector'
          write (rdluerr,*) '  for model ',rdtag
          rderr=ibset(rderr,1)
        else
          ny=sign(1.0d0/sqrt(1.0d0+(absx/absy)**2),y)
          nx=(x/y)*ny
        endif
      endif
      return
      end

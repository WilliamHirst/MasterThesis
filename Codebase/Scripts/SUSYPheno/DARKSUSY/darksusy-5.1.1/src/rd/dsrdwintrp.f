************************************************************************
                      subroutine dsrdwintrp(wrate,unit)
c_______________________________________________________________________
c  write out a table of
c    initial cm momentum p
c    invariant annihilation rate w
c    integration variable u
c    integrand f
c    interpolated integrand g
c    interpolation relative error f/g-1
c  input:
c    unit - logical unit to write on (integer)
c    wrate - invariant annihilation rate (real, external)
c  common:
c    'dsmssm.h' - file with susy common blocks
c    'dsrdcom.h' - included common blocks
c  uses dsrdtab,dsrdfunc
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsrdcom.h'
      real*8 u,du,utop,p,w,wrate,f,x,dsrdfunc,u2,
     &  y,g,dsrdwintp
      external wrate,dsrdwintp
      integer unit,i,n
      parameter (n=500)
c------------------------------------------------------------------------
      if (unit.le.0) return
      x=2.0d0
      mco(1)=abs(mass(lsp))
      utop=5.0d0
      du=utop/n
      call dsrdtab(wrate,x)
      nlo=1
      nhi=2
      write (unit,*) '  p  w  u  f  g  f/g-1'
      do 10 i=0,n
        u=i*du
        u2=u*u
        y=1.0d0+0.5d0*u2/x
        p=mco(1)*u*sqrt(0.5d0*(y+1.0d0)/x)
        w=wrate(p)
        f=dsrdfunc(u,x,wrate)
        g=dsrdfunc(u,x,dsrdwintp)
        write (unit,1000) p,w,u,f,g,g/f-1.0d0
   10 continue
 1000 format (6(2x,e12.6))
      end









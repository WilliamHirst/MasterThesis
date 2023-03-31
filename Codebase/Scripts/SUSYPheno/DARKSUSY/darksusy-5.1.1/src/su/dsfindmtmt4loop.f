      subroutine dsfindmtmt4loop
      implicit none
      include 'dsmssm.h'
      real*8 c,x,delta,a,dsralph34loop
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      mtmt=0.9999d0*mass(kt) ! any value < mass(kt) to have nf=6 in alph3
      a=dsralph34loop(mass(kt))/pi
      x=mbmb/mass(kt)*(1.d0+dsralph34loop(mbmb)/pi*4.d0/3.d0)
      delta=x*(pi**2/8.d0+x*(-0.597d0+x*0.230d0))
      c=1.d0-4.d0/3.d0*a+a**2*(-6.478739534870565d0-4.d0/3.d0*delta)
     &     +a**3*(-60.27d0)
      mtmt=mass(kt)*c
c      write(*,*) 'for alph3 and mass(kt) = ',alph3,mass(kt)
c      write(*,*) 'mtmt (4 loops) = ',mtmt
      return
      end  

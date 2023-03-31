      real*8 function dsntdkyf(xi,xf)
c y_f according to damour (gamma^{tot}_a = \half\alpha y_f^2 )
c lb 1998-02-24
      implicit none
      real*8 xf,dsai,dsbi,dsaip,dsbip,num,den,xi
c      write(*,*) 'xi=',xi,'  xf=',xf
      num=dsaip(xf)*dsbip(xi)-dsbip(xf)*dsaip(xi)
      den=dsai(xf)*dsbip(xi)-dsbi(xf)*dsaip(xi)
      if (den.eq.0.d0) then
         write(*,*) 'zero denominator in dsntdkyf! '
         stop
      else
            dsntdkyf=num/den
      endif
      return
      end

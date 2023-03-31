      subroutine dsneusct
c_______________________________________________________________________
c  neutralino masses and mixings. base is b-ino, w3-ino, h1-ino, h2-ino.
c  uses quartic.
c  called by susyin or mrkin.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c  history:
c     940528 readability improvement (pg)
c     950316 order by increasing mass (pg)
c     951110 positive mass convention (pg)
c     970211 loop corrections via switch neuloop (joakim edsjo)
c     990724 drop v1,v2 and change mass scale to max(mz,m1,m2,mu) (pg)
c     080606 use T. Hahn Diag routines (pg)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      real*8 mx(4),a3,a2,a1,a0,x,w,mm
      real*8 s1
      real*8 neuzmx(4,4),delta_33,delta_44,qscale,hb,ht,sin2thb,sin2tht
      complex*16 dsb0loop,m1m,m2m,mum,mzm,m(4,4)
      integer i, j, k(4), tmpk
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      external dsb0loop

c-----------------------------------------------------------------------
      mm = max(mass(kz),abs(m1),abs(m2),abs(mu))
      m1m = dcmplx(m1/mm,0.d0)
      m2m = dcmplx(m2/mm,0.d0)
      mum = dcmplx(mu/mm,0.d0)
      mzm = dcmplx(mass(kz)/mm,0.d0)
c mass matrix
      m(1,1)=m1m
      m(1,2)=dcmplx(0.0d0,0.0d0)
      m(1,3)=-mzm*sinthw*cosbe
      m(1,4)=+mzm*sinthw*sinbe
      m(2,1)=m(1,2)
      m(2,2)=m2m
      m(2,3)=+mzm*costhw*cosbe
      m(2,4)=-mzm*costhw*sinbe
      m(3,1)=m(1,3)
      m(3,2)=m(2,3)
      m(3,3)=dcmplx(0.d0,0.d0)
      m(3,4)=-mum
      m(4,1)=m(1,4)
      m(4,2)=m(2,4)
      m(4,3)=m(3,4)
      m(4,4)=dcmplx(0.d0,0.d0)

c...drees, nojiri, roy and yamada loop corrections, hep-ph/9701219
c...je addition 97-02-11
      delta_33=0.d0
      delta_44=0.d0
      if (neuloop.ge.1) then
        qscale = abs(mu)
c        hb=g2weak*mass(kb)/(sqrt(2.0d0)*mass(kw)*cosbe)
c        ht=g2weak*mass(kt)/(sqrt(2.0d0)*mass(kw)*sinbe)
c... change pu 02-07-02
        hb=yukawa(kb)
        ht=yukawa(kt)
        sin2thb = -2.0d0*dreal(sqdlmx(3,3)*sqdlmx(6,3))
        sin2tht = -2.0d0*dreal(squlmx(3,3)*squlmx(6,3))
c... check what scale are mass(kb) and mass(kt) below
        delta_33 = -3.0d0/(16.0d0*pi**2)*hb**2*mass(kb)*sin2thb*
     &    dreal(dsb0loop(qscale,mass(kb),mass(ksqd(3)))-
     &    dsb0loop(qscale,mass(kb),mass(ksqd(6))))
        delta_44 = -3.0d0/(16.0d0*pi**2)*ht**2*mass(kt)*sin2tht*
     &    dreal(dsb0loop(qscale,mass(kt),mass(ksqu(3)))-
     &    dsb0loop(qscale,mass(kt),mass(ksqu(6))))
        m(3,3) = dcmplx(delta_33/mm)
        m(4,4) = dcmplx(delta_44/mm)
      endif

      if (prtlevel.ge.2) then
         write (*, '('' neutralino mass matrix ='',4(f10.2))')
     &        (mm*m(1,i),i=1,4)
         do j=2,4
            write (*, '(''                         '',4(f10.2))')
     &           (mm*m(j,i),i=1,4)
         enddo
      endif

c eigenvalues and eigenvectors
      call TakagiFactor(4,m,4,mx,neunmx,4,1)
      kln=1
      do i=1,4
         mass(kn(i))=mm*mx(i)
      enddo

      end

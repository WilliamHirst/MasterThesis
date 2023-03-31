      subroutine dsddffsi(q,a,z,ff)
c_______________________________________________________________________
c  Spin-independent form factor for direct detection.
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(2*M*E) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    ff : |F(q)|^2, the square of the form factor
c  author: paolo gondolo (paolo@physics.utah.edu) 2004-2008
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      include 'dsge.h'

      real*8 q
      integer a,z
      real*8 ff
      real*8 x,x2,fermiGeV,s,y,f,ha,cpa,hc,r,r1,mni,r1tmp
      real*8 atomicmassunit,mpu,mnu
      parameter (fermiGeV = 1.d0/0.1973269602d0)
      parameter (atomicmassunit=0.931494028d0)
      parameter (mpu=1.00727646688d0)
      parameter (mnu=1.0086649156d0)
      
      if (a.eq.1) then
         ff = 1.d0

c best available form factor

      else if (ddffsi.eq.'best') then
         call dsddfbff(q,a,z,ff)
         if (ff.eq.-1.d0) then
            if (prtlevel.gt.0) write (*,*) 
     &           'dsddffsi: switching to SOG....'
            call dsddsogff(q,a,z,ff)
            if (ff.eq.-1.d0) then
               if (prtlevel.gt.0) write (*,*) 
     &              'dsddffsi: switching to Fermi....'
               call dsddfermi(q,a,z,ff)
               if (ff.eq.-1.d0) then
                  if (prtlevel.gt.0) write (*,*) 
     &                 'dsddffsi: switching to Lewin-Smith....'
                  call dsddlsff(q,a,z,ff)
               endif
            endif
         endif
         
c no form factor

      else if (ddffsi.eq.'no-ff') then
         ff=1.d0

c Helm form factor with Lewin-Smith parameters

      else if (ddffsi.eq.'L-S') then
         call dsddlsff(q,a,z,ff)

c Helm form factor as in DarkSUSY 4.1 (with expansion of f fixed, PG 20080216)

      else if (ddffsi.eq.'ds4.1') then
         mni = (z*mpu+(a-z)*mnu)*atomicmassunit
         r = (0.91d0*exp(log(mni)/3.d0)+0.3d0)
         s = 1.d0
         r1tmp = r*r-5.d0*s*s
         if (r1tmp.gt.0.d0) then
            r1 = sqrt(r1tmp)
         else
            r1 = r
         endif 
         x = dabs(q*r1*fermiGeV)
         y = q*s*fermiGeV
         if (x.gt.5.d-8) then
            f = (3.d0*(sin(x)-x*cos(x))/x**3)
         else
            x2 = x*x
            f = 1.d0+x2*(-0.1d0+x2*(1.d0/280.d0+x2*(-1.d0/15120.d0+
     &           x2/1330560.d0)))
         endif
         ff = f**2*exp(-y**2)

c gaussian form factor

      else if (ddffsi.eq.'gauss') then
         r = 0.89d0*exp(log(dble(a))/3.d0)+0.3d0
         x = q*r*fermiGeV
         ff = exp(-x*x/3.d0)

cc Sum of Gaussian form factor                                cc

      else if (ddffsi.eq.'SOG') then
         call dsddsogff(q,a,z,ff)

cc Fourier-Bessel form factor                                 cc

      else if (ddffsi.eq.'FB') then
         call dsddfbff(q,a,z,ff)

cc Fermi Integration form factor                              cc

      else if (ddffsi.eq.'Fermi') then
         call dsddfermi(q,a,z,ff)

c unrecognized form factor label

      else
         write (*,*) 
     &        'Unrecognized spin-independent form factor',ddffsi
         write (*,*) 'DarkSUSY stops'
         stop
      endif
      
      return
      end

      subroutine dsddffsd(q,a,z,l2jjpp,l2jjnn,l2jjpn)
c_______________________________________________________________________
c  Spin-dependent structure functions for direct detection.
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(M*E/2/mu^2) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    l2jjnn, l2jjpp, l2jjpn : the factors \lambda^2 J (J+1) made functions
c       of the momentum transfer q, for pp, nn, and pn; they are
c       related to Spp, Snn, and Spn through
c           S = \lambda^2 J (J+1) (2J+1)/\pi
c       In turn, Spp, Snn, and Spn are related to the spin-dependent 
c       structure functions S_{00}(q), S_{01}(q), S_{11}(q) defined in
c       Engel (PLB264,114,1991) via
c           Spp = S00+S11+S01
c           Snn = S00+S11-S01
c           Spn = 2*(S00-S11)
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c  modified: pg 040605 switched s01 and s11 in Na-23
c  modified: pg 080217 changed to \lambda^2 J (J+1)
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      include 'dsge.h'

      real*8 q
      integer a,z
      real*8 l2jjpp,l2jjnn,l2jjpn
      real*8 fermiGeV,y,j
      parameter (fermiGeV = 1.d0/0.1973269602d0)
      
      if (a.eq.1.and.z.eq.1) then ! proton
         l2jjpp = 0.75d0
         l2jjnn = 0.d0
         l2jjpn = 0.d0

      else if (a.eq.1.and.z.eq.0) then ! neutron
         l2jjpp = 0.d0
         l2jjnn = 0.75d0
         l2jjpn = 0.d0

c best available form factor

      else if (ddffsd(2:).eq.'best') then
         if (ddffsd(1:1).eq.'1') then
            call dsddismff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
         else
            call dsddismff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
         endif
         if (l2jjpp.eq.-1.d0) then
            if (prtlevel.gt.0) write (*,*) 
     &           'dsddffsd: switching to OddG....'
            if (ddffsd(1:1).eq.'1') then
               call dsddoddgff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
            else
               call dsddoddgff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
            endif
            if (l2jjpp.eq.-1.d0) then
               if (prtlevel.gt.0) write (*,*) 
     &              'dsddffsd: switching to SPSM....'
               if (ddffsd(1:1).eq.'1') then
                  call dsddspsmff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
               else
                  call dsddspsmff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
               endif
               if (l2jjpp.eq.-1.d0) then
                  if (prtlevel.gt.0) write (*,*) 
     &                 'dsddffsd: giving up.... sigma set to zero'
                  l2jjpp=0.d0
                  l2jjnn=0.d0
                  l2jjpn=0.d0
               endif
            endif
         endif
         
c interacting shell model

      else if (ddffsd(2:).eq.'ISM'.or.ddffsd(2:).eq.'ISMR') then
         if (ddffsd(1:1).eq.'1') then ! with form factor = 1
            call dsddismff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
         else
            call dsddismff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
         endif

c single particle shell model

      else if (ddffsd(2:).eq.'SPSM') then
         if (ddffsd(1:1).eq.'1') then ! with form factor = 1
            call dsddspsmff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
         else
            call dsddspsmff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
         endif

c odd-group model

      else if (ddffsd(2:).eq.'OddG') then
         if (ddffsd(1:1).eq.'1') then ! with form factor = 1
            call dsddoddgff(0.d0,a,z,l2jjpp,l2jjnn,l2jjpn)
         else
            call dsddoddgff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
         endif

c unrecognized form factor label

      else
         write (*,*) 
     &        'Unrecognized spin-dependent spin or form factor',ddffsd
         write (*,*) 'DarkSUSY stops'
         stop
      endif
      
      return
      end

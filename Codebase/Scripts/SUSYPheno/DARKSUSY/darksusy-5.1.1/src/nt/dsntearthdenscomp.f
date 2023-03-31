

***********************************************************************
*** dsntearthdenscomp gives the number density of nucleons of mass
*** number a per cm^3
*** input: radius - in meters
***        mass number - 16 for o
***                      24 for mg
***                      28 for si
***                      56 for fe JE UPDATE!
*** the radius should be given in m and the density is returned in
*** g/cm^3
*** author: joakim edsjo
*** date: april 6, 1999
*** Updated with new values by J. Edsjo, 2003-11-21
***********************************************************************

      real*8 function dsntearthdenscomp(r,a)
      implicit none
      include 'dsmpconst.h'
      real*8 r,m_a,dsntearthdens,rcore,fraction
      integer a,index
      parameter (rcore=3483.0d3)           ! division between core and mantle

      real*8 mfrac(2,11) ! core and mantle

c...Best values from Encyclopedia Britannica
c      data mfrac/0.05d0  ,0.4366d0,   ! O
c     &           0.0d0   ,0.2100d0,   ! Si
c     &           0.0d0   ,0.2304d0,   ! Mg
c     &           0.90d0  ,0.0591d0,   ! Fe
c     &           0.0d0   ,0.0d0,      ! Ca
c     &           0.0d0   ,0.0d0,      ! P
c     &           0.0d0   ,0.0d0,      ! Na
c     &           0.0d0   ,0.0d0,      ! S
c     &           0.0d0   ,0.0d0,      ! Ni
c     &           0.0d0   ,0.0d0,      ! Al
c     &           0.0d0   ,0.0d0/      ! Cr

c...Values from McDonough, Treatise on Geochemistry, Vol 2, Elsevier, 2003.
c...These values are very close to those in The Encyclopedia of
c...Geochemistry, Eds. Marshall and Fairbridge, Klower Acadmic Publ, 1998.
      data mfrac/0.00d0  ,0.44d0,     ! O
     &           0.060d0 ,0.2100d0,   ! Si
     &           0.0d0   ,0.228d0,    ! Mg
     &           0.855d0 ,0.0626d0,   ! Fe
     &           0.0d0   ,0.0253d0,   ! Ca
     &           0.0020d0,0.00009d0,  ! P
     &           0.0d0   ,0.0027d0,   ! Na
     &           0.019d0 ,0.00025d0,  ! S
     &           0.052d0 ,0.00196d0,  ! Ni
     &           0.0d0   ,0.0235d0,   ! Al
     &           0.009d0 ,0.0026d0/      ! Cr

c...Values from jkg translated to core and mantle values normalized to
c...give the same total mass fractions. With these values, dsntcapearthfull
c...and dsntcapearthnumi (with the same gaussian) are expeted to give
c...roughly the same result, and they do within about 0.5%
c      data mfrac/0.0d0   ,0.4444d0,   ! O
c     &           0.0d0   ,0.2222d0,   ! Si
c     &           0.0d0   ,0.2074d0,   ! Mg
c     &           0.9230d0,0.0d0,      ! Fe
c     &           0.0d0   ,0.0222d0,   ! Ca
c     &           0.0d0   ,0.0163d0,   ! P
c     &           0.0d0   ,0.0059d0,   ! Na
c     &           0.1539d0,0.0d0,      ! S
c     &           0.0923d0,0.0d0,      ! Ni
c     &           0.0d0   ,0.0d0,      ! Al
c     &           0.0d0   ,0.0d0/      ! Cr

      if (a.eq.16) then ! O
        index=1
        m_a=16.0d0
      elseif (a.eq.28) then ! Si
        index=2
        m_a=28.1d0
      elseif (a.eq.24) then ! Mg
        index=3
        m_a=24.3d0
      elseif (a.eq.56) then ! Fe
        index=4
        m_a=55.85d0
      elseif (a.eq.40) then ! Ca
        index=5
        m_a=40.0d0
      elseif (a.eq.30) then ! P
        index=6
        m_a=30.0d0
      elseif (a.eq.23) then ! Na
        index=7
        m_a=23.0d0
      elseif (a.eq.32) then ! S
        index=8
        m_a=32.0d0
      elseif (a.eq.59) then ! Ni
        index=9
        m_a=59.0d0
      elseif (a.eq.27) then ! Al
        index=10
        m_a=27.0d0
      elseif (a.eq.52) then ! Cr
        index=11
        m_a=52.0d0
      else
        write(*,*) 'warning in dsntearthdenscomp:'
        write(*,*) '  element ',a,' requested - not available!'
        write(*,*) '  please update dsntearthdenscomp'
        dsntearthdenscomp=0.0d0
        return
      endif

c...m_a is now in atomic units u which is what we want

      if (r.lt.rcore) then
        fraction=mfrac(1,index)
      else
        fraction=mfrac(2,index)
      endif

      dsntearthdenscomp=fraction*dsntearthdens(r)/m_a*n_avogadro ! number/cm^3

      return
      end

      function dsandwdcosopt(p,costheta)
c_______________________________________________________________________
c  annihilation differential invariant rate.
c  input:
c    p - initial cm momentum (real) for lsp annihilations
c    costheta - cosine of c.m. annihilation angle
c  common:
c    'dssusy.h' - file with susy common blocks
c  uses dsandwdcosnn, dsandwdcoscn and dsandwdcoscc
c  called by wx when dwopt=true
c  author: joakim edsjo (edsjo@physto.se)
c  date: 97-09-09
c  modified: 98-05-01
c=======================================================================
      implicit none
      include 'dssusy.h'
      include 'dsandwcom.h'
      real*8 dsandwdcosopt,dsandwdcosnn,dsandwdcoscn,dsandwdcoscc,mx
      real*8 w,p,costheta,pnew,s,mp1,mp2,f,tmp,partmax
      integer i,j,k,kp1(36),kp2(36),nch,ch,kpn(4),nn,kpcha(2),ncha
      integer imax,jmax,kmax

      do i=1,6
        do j=1,6
          do k=1,54
            parts(i,j,k)=0.0d0
            partincl(i,j,k)=.true.
            tur(i,j,k)=.false.        ! t- or u-channel resonance
          enddo
        enddo
      enddo

c... check for t- or u-channel resonances
c... neutralino-neutralino
      do i=1,4
        do j=i,4
          call dsantunn(p,i,j)
        enddo
      enddo

c... chargino-neutralino
      do i=5,6
        do j=1,4
          call dsantucn(p,i,j)
        enddo
      enddo

c... chargino-chargino
      do i=5,6
        do j=i,6
          call dsantucc(p,i,j)
        enddo
      enddo

c      write(*,*) 'dsandwdcosopt called with p = ',p
c      write(*,*) '            costheta = ',costheta
c...initialize akinvar variables
      call dsankinvar(0.0d0,0.0d0,0,0,0,0,0)

      if (coann.ne.0) then
c...decide which neutralinos to include
        kpn(1)=kn(1)
        nn=1
        mx = mass(kn(1))
        do i=2,4
          if (mass(kn(i))/mx.le.mcofr) then
            nn=nn+1
            kpn(nn)=kn(i)
          endif
        enddo
c...decide which charginos to include
        ncha=0
        do i=1,2
          if (mass(kcha(i))/mx.le.mcofr) then
            ncha=ncha+1
            kpcha(ncha)=kcha(i)
          endif
        enddo

c---------------------------------- neutralino-neutralino annihilation
        nch=0
        do i=1,nn
          do j=i,nn
            nch=nch+1
            kp1(nch)=kpn(i)
            kp2(nch)=kpn(j)
          enddo
        enddo

        w = dsandwdcosnn(p,costheta,kp1(1),kp2(1))
        do k=1,28
          parts(1,1,k)=prtial(k)
        enddo

c...convert given p for lsp:s to pnew and call dsandwdcosnn with pnew
        s = 4*(mass(kn(1))**2+p**2)
        do ch=2,nch
          mp1 = mass(kp1(ch))
          mp2 = mass(kp2(ch))
          tmp=(s-(mp1-mp2)**2)*(s-(mp1+mp2)**2)/(4.0d0*s)
          if (tmp.gt.0.0d0) then
            pnew = sqrt(tmp)
            f=1.0d0
            if (kp1(ch).ne.kp2(ch)) f=2.0d0*f
            w=w+f*pnew/p*dsandwdcosnn(pnew,costheta,kp1(ch),kp2(ch))
            do k=1,29
              parts(kp1(ch)-kn(1)+1,kp2(ch)-kn(1)+1,k)=prtial(k)
              parts(kp2(ch)-kn(1)+1,kp1(ch)-kn(1)+1,k)=prtial(k)
            enddo
          endif
        enddo

        if (ncha.gt.0) then
c------------------------------------ neutralino-chargino annihilation
          nch=0
          do i=1,ncha
            do j=1,nn      ! i->1 je correction 97-05-12
              nch=nch+1
              kp1(nch)=kpcha(i)
              kp2(nch)=kpn(j)
            enddo
          enddo
c...convert given p for lsp:s to pnew and call dsandwdcosii with pnew
          s = 4*(mass(kn(1))**2+p**2)
          do ch=1,nch
            mp1 = mass(kp1(ch))
            mp2 = mass(kp2(ch))
            tmp=(s-(mp1-mp2)**2)*(s-(mp1+mp2)**2)/(4.0d0*s)
            if (tmp.gt.0.0d0) then
              pnew = sqrt(tmp)
              f=2.0d0     ! degrees of freedom for chargino=4
              if (kp1(ch).ne.kp2(ch)) f=2.0d0*f
              w=w+f*pnew/p*dsandwdcoscn(pnew,costheta,kp1(ch),kp2(ch))
              do k=1,28
                parts(kp2(ch)-kn(1)+1,kp1(ch)-kcha(1)+5,k)=prtial(k)
                parts(kp1(ch)-kcha(1)+5,kp2(ch)-kn(1)+1,k)=prtial(k)
              enddo
            endif
          enddo

c-------------------------------------- chargino-chargino annihilation
          nch=0
          do i=1,ncha
            do j=i,ncha
              nch=nch+1
              kp1(nch)=kpcha(i)
              kp2(nch)=kpcha(j)
            enddo
          enddo
c...convert given p for lsp:s to pnew and call dsandwdcoscc with pnew
          s = 4*(mass(kn(1))**2+p**2)
          do ch=1,nch
            mp1 = mass(kp1(ch))
            mp2 = mass(kp2(ch))
            tmp=(s-(mp1-mp2)**2)*(s-(mp1+mp2)**2)/(4.0d0*s)
            if (tmp.gt.0.0d0) then
              pnew = sqrt(tmp)
              f=4.0d0     ! degrees of freedom for chargino=4
              if (kp1(ch).ne.kp2(ch)) f=2.0d0*f
              w=w+f*pnew/p*dsandwdcoscc(pnew,costheta,kp1(ch),kp2(ch))
              do k=1,54
                parts(kp1(ch)-kcha(1)+5,kp2(ch)-kcha(1)+5,k)=prtial(k)
                parts(kp2(ch)-kcha(1)+5,kp1(ch)-kcha(1)+5,k)=prtial(k)
              enddo
            endif
          enddo
        endif      ! ncha.gt.0

      else       ! no coannihilation
        w = dsandwdcosnn(p,costheta,kn(1),kn(1))
        do k=1,29
          parts(1,1,k)=prtial(k)
        enddo
      endif

      dsandwdcosopt=w
c      write(33,*) costheta,w    ! je test

c... optimization stuff
      partmax=0.0d0
      do i=1,6
        do j=1,6
          do k=1,54
            if ((.not.tur(i,j,k)).and.parts(i,j,k).gt.partmax) then
              partmax=parts(i,j,k)
              imax=i
              jmax=j
              kmax=k
            endif
          enddo
        enddo
      enddo

      do i=1,6
        do j=1,6
          do k=1,54
            if ((.not.tur(i,j,k)).and.parts(i,j,k)/partmax.lt.brmin)
     &        partincl(i,j,k)=.false.
          enddo
        enddo
      enddo

c      if (p.gt.3*mass(kn(1))) then
c        write(*,*) 'partmax = ',partmax
c        do k=1,28
c          write(*,*) '  parts(1,1,',k,') = ',parts(1,1,k),
c     &      partincl(1,1,k)
c        enddo
c        stop
c      endif

c       write(*,*)
c       write(*,*) 'dsandwdcosopt called with p = ',p
c       write(*,*) '  maximum partial = ',partmax,' for:'
c       write(*,*) '  i=',imax,' j=',jmax,' k=',kmax

      end








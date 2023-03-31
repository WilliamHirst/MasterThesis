      subroutine dsaswcomp(p,costheta,kp1,kp2,dwbsmax,dmbsmax)
c_______________________________________________________________________
c  Routine to compare annihilation cross sections to find
c  out how big they are
c    p - initial cm momentum (real) for lsp annihilations
c    costheta - cosine of c.m. annihilation angle
c  uses dsandwdcosnn, dsandwdcoscn and dsandwdcoscc
c  author: joakim edsjo (edsjo@physto.se)
c  date: 02-05-23
c  modified: Joakim Edsjo, 02-10-22
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      real*8 dsandwdcos,dsandwdcosnn,dsandwdcoscn,dsandwdcoscc,
     &  dsasdwdcossfsf,dsasdwdcossfchi
      real*8 dw(30,30),dw0,dwmax,dwbsmax,dmbsmax
      real*8 mx,bsup,peff,sqrts
      real*8 p,costheta,pnew,s,mp1,mp2,f,tmp
      integer kpi(30),type(30),imax,jmax,ibsmax,jbsmax
      integer i,j,kp1,kp2

c      write(*,*) 'dsandwdcos called with p = ',p
c      write(*,*) '            costheta = ',costheta
c...initialize dsankinvar variables
      call dsankinvar(0.0d0,0.0d0,0,0,0,0,0)

      mx = mass(kn(1))

c---------------------------------------------------------------------
c----- OK, now calculate dwdcos for all combinations of initial states
c---------------------------------------------------------------------

      kpi(1)=kn(1)
      kpi(2)=kn(2)
      kpi(3)=kn(3)
      kpi(4)=kn(4)
      do i=1,4
        type(i)=1  ! neutralino
      enddo

      kpi(5)=kcha(1)
      kpi(6)=kcha(2)
      do i=5,6
        type(i)=2  ! chargino
      enddo

      kpi(7)=ksl(1)
      kpi(8)=ksl(2)
      kpi(9)=ksl(3)
      kpi(10)=ksl(4)
      kpi(11)=ksl(5)
      kpi(12)=ksl(6)

      kpi(13)=ksnu(1)
      kpi(14)=ksnu(2)
      kpi(15)=ksnu(3)

      kpi(16)=ksqu(1)
      kpi(17)=ksqu(2)
      kpi(18)=ksqu(3)
      kpi(19)=ksqu(4)
      kpi(20)=ksqu(5)
      kpi(21)=ksqu(6)

      kpi(22)=ksqd(1)
      kpi(23)=ksqd(2)
      kpi(24)=ksqd(3)
      kpi(25)=ksqd(4)
      kpi(26)=ksqd(5)
      kpi(27)=ksqd(6)
      do i=7,27
        type(i)=3  ! sfermion
      enddo

c...Go on and calculate invariant rates
      dw0=0.0d0
      dwmax=0.0d0
      imax=0
      jmax=0
      dwbsmax=0.0d0
      ibsmax=0
      jbsmax=0
      dmbsmax=0.0d0
      do i=1,27
        do j=i,27
          if (type(i).eq.1.and.type(j).eq.1) then
            dw(i,j)=dsandwdcosnn(p,costheta,kpi(i),kpi(j))
          elseif (type(i).eq.1.and.type(j).eq.2) then
            dw(i,j)=dsandwdcoscn(p,costheta,kpi(j),kpi(i))
          elseif (type(i).eq.1.and.type(j).eq.3) then
            dw(i,j)=dsasdwdcossfchi(p,costheta,kpi(j),kpi(i))
          elseif (type(i).eq.2.and.type(j).eq.2) then
            dw(i,j)=dsandwdcoscc(p,costheta,kpi(j),kpi(i))
          elseif (type(i).eq.2.and.type(j).eq.3) then
            dw(i,j)=dsasdwdcossfchi(p,costheta,kpi(j),kpi(i))
          elseif (type(i).eq.3.and.type(j).eq.3) then
            dw(i,j)=dsasdwdcossfsf(p,costheta,kpi(i),kpi(j))
          else
            write(*,*) 'ERROR in dsaswcomp: Shouldn''t get here...'
            stop
          endif
          dw(j,i)=dw(i,j)
          if (i.eq.1.and.j.eq.1) dw0=dw(1,1)
          if (dw(i,j).gt.dwmax) then
            dwmax=dw(i,j)
            imax=i
            jmax=j
          endif
          if (i.ne.1.and.j.ne.1) then
c...Now include the Boltzmann suppression and factors of p in 
c...thermal average. No degrees of freedom included
            bsup=exp(-20.0d0*
     &        (mass(kpi(i))+mass(kpi(j))-2.0d0*mx)/mx)
            peff=sqrt((sqrt(mass(kpi(i))**2+p**2)+
     &           sqrt(mass(kpi(j))**2+p**2))**2/4-mx**2)
            bsup=bsup*peff/p
            if (dw(i,j)*bsup.gt.dwbsmax) then
              dwbsmax=dw(i,j)*bsup
              ibsmax=i
              jbsmax=j
              dmbsmax=max(mass(kpi(i))/mx,mass(kpi(j))/mx)
c              write(*,*) 'i=',i,'  j=',j
c              write(*,*) 'kp1=',kpi(i),'  kp2=',kpi(j)
c              write(*,*) 'mx=',mx
c              write(*,*) 'm1=',mass(kpi(i))
c              write(*,*) 'm2=',mass(kpi(j))
c              write(*,*) 'dmbsmax=',dmbsmax
c              write(*,*) 'dwbsmax=',dwbsmax
            endif
          endif
        enddo
      enddo

c      write(*,*) ' '
c      write(*,*) 'Model: ',idtag
c      write(*,*) 'dwmax/dw0=',dwmax/dw0,' i=',imax,' j=',jmax
c      write(*,*) '    kp1=',kpi(imax),' kp2=',kpi(jmax)
c      write(*,*) 'dwbsmax/dw0=',dwbsmax/dw0,' i=',ibsmax,' j=',jbsmax
c      write(*,*) '    kp1=',kpi(ibsmax),' kp2=',kpi(jbsmax),
c     &  ' dmbsmax=',dmbsmax
c      write(*,*) 'm1=',mass(kpi(ibsmax)),'  m2=',mass(kpi(jbsmax))

      kp1=kpi(ibsmax)
      kp2=kpi(jbsmax)
      dwbsmax=dwbsmax/dw0

      end



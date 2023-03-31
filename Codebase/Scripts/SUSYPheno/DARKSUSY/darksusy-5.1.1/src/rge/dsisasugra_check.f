      subroutine dsisasugra_check(valid)
c=======================================================================
c  This routine checks that the neutralino, chargino and neutralino
c  mixing matrices extracted from ISASUGRA are consistent with the
c  DarkSUSY convention. It does this by checking that they really
c  do diagonlize the mass matrices. If any inconsistencies (apart
c  from small numerical differences are found), an error message
c  is written.
c  Author: J. Edsjo and M. Schelke, 2002-12-03
c=======================================================================
      implicit none

      include 'dsmssm.h'
      include 'dsidtag.h'
      integer valid
      integer i,n,j,k
      real*8 p,w,dsanwx,mz,dsabsq,mw
      complex*16 um(2,2),mdiag(2,2)
      complex*16 neudiag(4,4),tmp(4,4)
      real*8 sfm(2,2),R(2,2),stmp(2,2),m2ll,m2lr,m2rr,d,mf,merr

      include 'dsisasugra.h'

c======================================== Check chargino matrices
      mw=amw  ! Use ISASUGRA m_W value for this check

      um(1,1)=chaumx(1,1)*m2
     &        +chaumx(1,2)*sqrt(2.d0)*mw*cosbe
      um(1,2)=chaumx(1,1)*sqrt(2.d0)*mw*sinbe 
     &        +chaumx(1,2)*mu
      um(2,1)=chaumx(2,1)*m2
     &        +chaumx(2,2)*sqrt(2.d0)*mw*cosbe
      um(2,2)=chaumx(2,1)*sqrt(2.d0)*mw*sinbe 
     &        +chaumx(2,2)*mu

      mdiag(1,1)=chavmx(1,1)*um(1,1)+chavmx(1,2)*um(1,2)
      mdiag(1,2)=chavmx(2,1)*um(1,1)+chavmx(2,2)*um(1,2)
      mdiag(2,1)=chavmx(1,1)*um(2,1)+chavmx(1,2)*um(2,2)
      mdiag(2,2)=chavmx(2,1)*um(2,1)+chavmx(2,2)*um(2,2)           

      if(dabs(dble(mdiag(1,1))-mass(kcha1)).gt.1.d0) then
        write(*,*) 'ERROR in dsisasugra_check: ',
     &   'chargino mass eigenvalue #1 wrong'
        write(*,*) '  Model: ',idtag
        write(*,*) '  Mass eigenvalue from ISASUGRA: ',mass(kcha1)
        write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &   dble(mdiag(1,1))
c        valid=10
      endif

      if(dabs(dble(mdiag(2,2))-mass(kcha2)).gt.1.d0) then
        write(*,*) 'ERROR in dsisasugra_check: ',
     &   'chargino mass eigenvalue #2 wrong'
        write(*,*) '  Model: ',idtag
        write(*,*) '  Mass eigenvalue from ISASUGRA: ',mass(kcha2)
        write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &   dble(mdiag(2,2))
c        valid=10
      endif

      if(dabs(dble(mdiag(1,2))).gt.1.d0) then
        write(*,*) 'ERROR in dsisasugra_check: ',
     &   'chargino off-diagonal element non-zero'
        write(*,*) '  Model: ',idtag
        write(*,*) '  Entry (1,2) from diagonalized matrix: ',
     &   dble(mdiag(1,2))
c        valid=10
      endif

      if(dabs(dble(mdiag(2,1))).gt.1.d0) then
        write(*,*) 'ERROR in dsisasugra_check: ',
     &   'chargino off-diagonal element non-zero'
        write(*,*) '  Model: ',idtag
        write(*,*) '  Entry (2,1) from diagonalized matrix: ',
     &   dble(mdiag(2,1))
c        valid=10
      endif

c...Write out matrix (uncomment if needed)
c      do i=1,2
c        do j=1,2
c          write(*,*) 'cha i=',i,' j=',j,' mdiag=',mdiag(i,j)
c        enddo
c      enddo

c======================================== Check neutralino matrices

c...Define neutralino mass matrix
      mz=mass(kz)
      do i=1,4
        do j=1,4
          neudiag(i,j)=dcmplx(0.0d0,0.0d0)
        enddo
      enddo
      neudiag(1,1)=m1
      neudiag(2,2)=m2
      neudiag(3,4)=-mu
      neudiag(4,3)=-mu
      neudiag(1,3)=-mz*cosbe*sinthw
      neudiag(3,1)=neudiag(1,3)
      neudiag(1,4)=mz*sinbe*sinthw
      neudiag(4,1)=neudiag(1,4)
      neudiag(2,3)=mz*cosbe*costhw
      neudiag(3,2)=neudiag(2,3)
      neudiag(2,4)=-mz*sinbe*costhw
      neudiag(4,2)=neudiag(2,4)

c...Now perform N* . neudiag . N^dagger
c...Start with N* . neudiag
      do i=1,4
        do j=1,4
          tmp(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,4
            tmp(i,j)=tmp(i,j)+conjg(neunmx(i,k))*neudiag(k,j)
          enddo
        enddo
      enddo

c...Then do tmp . N^dagger
      do i=1,4
        do j=1,4
          neudiag(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,4
            neudiag(i,j)=neudiag(i,j)+
     &        tmp(i,k)*conjg(neunmx(j,k))
          enddo
        enddo
      enddo

c...Check matrix
      do i=1,4
        do j=1,4
          if (i.eq.j) then
            if (dabs(dble(neudiag(i,j))-mass(kn(i))).gt.1.d0) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &         'neutralino mass eigenvalue #',i,' wrong'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Mass eigenvalue from ISASUGRA: ',mass(kn(i))
              write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &          dabs(dble(neudiag(i,j))) 
c              valid=11
           endif
          else
            if (sqrt(dsabsq(neudiag(i,j))).gt.1.0d0) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &          'neutralino off-diagonal element non-zero'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Entry (',i,',',j,
     &          ') from diagonalized matrix: ',
     &          neudiag(i,j)
c              valid=11
            endif
          endif
c          write(*,*) 'neu i=',i,' j=',j,' neudiag=',neudiag(i,j)
        enddo
      enddo

c======================================== Check up-squark matrices

      d = mass(kz)**2*(1.d0-tanbe*tanbe)/(1.d0+tanbe*tanbe) ! mz^2 cos(2beta)

c...The matrix that would diagnolize the mass matrix is R and
c...R^T M R = diag(m_1,m_2)
c...Define R
      R(1,1)=squlmx(3,3)
      R(1,2)=squlmx(6,3)
      R(2,1)=squrmx(3,3)
      R(2,2)=squrmx(6,3)

c...Define the mass matrix
      mf=mtq  ! top mass at weak scale, used by isasgura
      m2ll = mass2q(3)+mf**2+d*(0.5d0-2.d0/3.d0*s2thw)
      m2rr = mass2u(3)+mf**2+d*(2.d0/3.d0*s2thw)
      m2lr = mf*(asoftu(3)-mu/tanbe)
      sfm(1,1)=m2ll
      sfm(1,2)=m2lr
      sfm(2,1)=m2lr
      sfm(2,2)=m2rr

c...Now perform R^T . sfm . R
c...Start with R^T . sfm
      do i=1,2
        do j=1,2
          stmp(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            stmp(i,j)=stmp(i,j)+R(k,i)*sfm(k,j)
          enddo
        enddo
      enddo
c...Then do stmp . R
      do i=1,2
        do j=1,2
          sfm(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            sfm(i,j)=sfm(i,j)+
     &        stmp(i,k)*R(k,j)
          enddo
        enddo
      enddo

c...Check matrix
c...Maximal error allowed
      merr=min(mass(ksqu(3))**2,mass(ksqu(6))**2)*0.05d0
      merr=max(20.0d0,merr)
      do i=1,2
        do j=1,2
          if (i.eq.j) then
            if (dabs(dble(sfm(i,j))-mass(ksqu(3*i))**2).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &         'top-squark mass eigenvalue #',i,' wrong'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Mass eigenvalue from ISASUGRA: ',
     &          mass(ksqu(3*i))
              write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=12
            endif
          else
            if (dabs(sfm(i,j)).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &          'top-squark off-diagonal element non-zero'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Sqrt of entry (',i,',',j,
     &          ') from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=12
            endif
          endif
c          write(*,*) 'squ i=',i,' j=',j,' sfm=',sfm(i,j)
        enddo
      enddo

c======================================== Check down-squark matrices

c...The matrix that would diagnolize the mass matrix is R and
c...R^T M R = diag(m_1,m_2)
c...Define R
      R(1,1)=sqdlmx(3,3)
      R(1,2)=sqdlmx(6,3)
      R(2,1)=sqdrmx(3,3)
      R(2,2)=sqdrmx(6,3)

c...Define the mass matrix
      mf=mbq  ! bottom mass at weak scale, used by isasgura
      m2ll = mass2q(3)+mf**2+d*(-0.5d0+1.d0/3.d0*s2thw)
      m2rr = mass2d(3)+mf**2+d*(-1.d0/3.d0*s2thw)
      m2lr = mf*(asoftd(3)-mu*tanbe)
      sfm(1,1)=m2ll
      sfm(1,2)=m2lr
      sfm(2,1)=m2lr
      sfm(2,2)=m2rr

c...Now perform R^T . sfm . R
c...Start with R^T . sfm
      do i=1,2
        do j=1,2
          stmp(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            stmp(i,j)=stmp(i,j)+R(k,i)*sfm(k,j)
          enddo
        enddo
      enddo
c...Then do stmp . R
      do i=1,2
        do j=1,2
          sfm(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            sfm(i,j)=sfm(i,j)+
     &        stmp(i,k)*R(k,j)
          enddo
        enddo
      enddo

c...Check matrix
c...Maximal error allowed
      merr=min(mass(ksqd(3))**2,mass(ksqd(6))**2)*0.05d0
      merr=max(20.0d0,merr)
      do i=1,2
        do j=1,2
          if (i.eq.j) then
            if (dabs(dble(sfm(i,j))-mass(ksqd(3*i))**2).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &         'bottom-squark mass eigenvalue #',i,' wrong'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Mass eigenvalue from ISASUGRA: ',
     &          mass(ksqd(3*i))
              write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=13
            endif
          else
            if (dabs(sfm(i,j)).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &          'bottom-squark off-diagonal element non-zero'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Sqrt of entry (',i,',',j,
     &          ') from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=13
            endif
          endif
c          write(*,*) 'sqd i=',i,' j=',j,' sfm=',sfm(i,j)
        enddo
      enddo

c======================================== Check slepton matrices

c...The matrix that would diagnolize the mass matrix is R and
c...R^T M R = diag(m_1,m_2)
c...Define R
      R(1,1)=sldlmx(3,3)
      R(1,2)=sldlmx(6,3)
      R(2,1)=sldrmx(3,3)
      R(2,2)=sldrmx(6,3)

c...Define the mass matrix
      mf=mlq  ! tau mass at the weak scale
      m2ll = mass2l(3)+mf**2+d*(-0.5d0+s2thw)
      m2rr = mass2e(3)+mf**2+d*(-s2thw)
      m2lr = mf*(asofte(3)-mu*tanbe)
      sfm(1,1)=m2ll
      sfm(1,2)=m2lr
      sfm(2,1)=m2lr
      sfm(2,2)=m2rr

c...Now perform R^T . sfm . R
c...Start with R^T . sfm
      do i=1,2
        do j=1,2
          stmp(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            stmp(i,j)=stmp(i,j)+R(k,i)*sfm(k,j)
          enddo
        enddo
      enddo
c...Then do stmp . R
      do i=1,2
        do j=1,2
          sfm(i,j)=dcmplx(0.0d0,0.0d0)
          do k=1,2
            sfm(i,j)=sfm(i,j)+
     &        stmp(i,k)*R(k,j)
          enddo
        enddo
      enddo

c...Check matrix
c...Maximal error allowed
      merr=min(mass(ksl(3))**2,mass(ksl(6))**2)*0.05d0
      merr=max(20.0d0,merr)
      do i=1,2
        do j=1,2
          if (i.eq.j) then
            if (dabs(dble(sfm(i,j))-mass(ksl(3*i))**2).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &         'charged slepton mass eigenvalue #',i,' wrong'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Mass eigenvalue from ISASUGRA: ',
     &          mass(ksl(3*i))
              write(*,*) '  Mass eigenvalue from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=14
            endif
          else
            if (dabs(sfm(i,j)).gt.merr) then
              write(*,*) 'ERROR in dsisasugra_check: ',
     &          'charged slepton off-diagonal element non-zero'
              write(*,*) '  Model: ',idtag
              write(*,*) '  Sqrt of entry (',i,',',j,
     &          ') from diagonalized matrix: ',
     &          sqrt(sfm(i,j))
c              valid=14
            endif
          endif
c          write(*,*) 'sl i=',i,' j=',j,' sfm=',sfm(i,j)
        enddo
      enddo

c-------------------------------------------------- finished checking

      return
      end


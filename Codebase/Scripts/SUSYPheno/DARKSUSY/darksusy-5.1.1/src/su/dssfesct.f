      subroutine dssfesct(unphys)
c_______________________________________________________________________
c  sfermion masses and mixings.
c  options:
c    msquarks=0  : squark masses and mixing from mass matrix
c    msquarks>0  : all squark masses equal to msquarks, no mixing
c    msquarks<0  : all squark masses = max(lsp,-msquarks), no mixing
c    msleptons=0  : squark masses and mixing from mass matrix
c    msleptons>0  : all squark masses equal to msleptons, no mixing
c    msleptons<0  : all squark masses = max(lsp,-msleptons), no mixing
c  called by susyin.
c  needs neusct.
c  author: paolo gondolo 1994-1999
c     981000 correction to diagonalization of mass matrices
c     941100 addition of generation-mixing for squarks
c     950100 correction to charged slepton mass matrix
c     990715 pg options msquarks and msleptons added
c     990724 pg drop chankowski constants
c     020219 pg better reporting of unphys
c     020405 pg matrix diagonalization rewritten to handle degenerate case
c     080611 pg 6x6 mass matrices
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      real*8 auxa,d,m2ll,m2rr,m2lr,msq1,msq2,theta,dsabsq,cmax
      complex*16 mm(6,6),mix(6,6)
      real*8 m(6),c(6,3)
      integer i,j,g,unphys,k,imax(6)
      character*80 message
      logical aux
      do g=1,3
        do i=1,3
          slulmx(i,g) = dcmplx(0.d0,0.d0)
        enddo
        do i=1,6
          sldlmx(i,g) = dcmplx(0.d0,0.d0)
          sldrmx(i,g) = dcmplx(0.d0,0.d0)
          squlmx(i,g) = dcmplx(0.d0,0.d0)
          squrmx(i,g) = dcmplx(0.d0,0.d0)
          sqdlmx(i,g) = dcmplx(0.d0,0.d0)
          sqdrmx(i,g) = dcmplx(0.d0,0.d0)
        enddo
        slulmx(g,g) = dcmplx(1.d0,0.d0)
        sldlmx(g,g) = dcmplx(1.d0,0.d0)
        sldrmx(g+3,g) = dcmplx(1.d0,0.d0)
        squlmx(g,g) = dcmplx(1.d0,0.d0)
        squrmx(g+3,g) = dcmplx(1.d0,0.d0)
        sqdlmx(g,g) = dcmplx(1.d0,0.d0)
        sqdrmx(g+3,g) = dcmplx(1.d0,0.d0)
      enddo
c------------------------------------- special options

      d = mass(kz)**2*(1.d0-tanbe*tanbe)/(1.d0+tanbe*tanbe) ! mz^2 cos(2beta)

      if (msleptons.gt.0.d0) then
         do i=1,6
            mass(ksnu(i)) = msleptons
            mass(ksl(i)) = msleptons
         enddo
         do i=1,3
            asofte(i) = mu*tanbe
         enddo
      else if (msleptons.lt.0.d0) then
         if (neuloop.ne.0) then
            call dswrite(0,1,
     &    'susfesct: neuloop set to zero since msfermions=mlsp etc')
            neuloop=0
         endif
         call dsneusct
         auxa = abs(mass(kn(kln)))
         if (auxa.lt.-msleptons) auxa = -msleptons
         do i=1,6
            mass(ksnu(i)) = auxa
            mass(ksl(i)) = auxa
         enddo
         do i=1,3
            asofte(i) = mu*tanbe
         enddo
      else
c-------------------------------------------------------------sneutrinos
         do i=1,3
            msq1 = mass2l(i)+d*0.5d0
            if (prtlevel.ge.2) then
               write (*,'('' sneutrino('',i1,'') mass^2 = '',f12.2)')
     &              i,msq1
            endif
            if (msq1.lt.0.0d0) then
               write (message,*) 'susfesct: negative ',pname(ksnu(i)),
     &              ' mass squared (',real(msq1),')'
               call dswrite(1,1,message)
               unphys = -10-i
            endif
            mass(ksnu(i)) = sqrt(abs(msq1))
            slulmx(i,i)=dcmplx(1.d0,0.d0)
         enddo
c--------------------------------------------------------charged sleptons
         do i=1,6
            do j=1,6
               mm(i,j)=dcmplx(0.d0,0.d0)
            enddo
         enddo
         do i=1,3
            m2ll = mass2l(i)+mass(kl(i))**2+d*(-0.5d0+s2thw)
            m2rr = mass2e(i)+mass(kl(i))**2+d*(-s2thw)
            m2lr = mass(kl(i))*(asofte(i)-mu*tanbe)
            if (prtlevel.ge.2) then
               write (*,
     &            '('' slepton('',i1,'') mass^2 matrix = '',2(f12.2))')
     &              i,m2ll,m2lr
               write (*,
     &              '(''                            '',2(f12.2))')
     &              m2lr,m2rr
            endif
            mm(2*i-1,2*i-1)=dcmplx(m2ll,0.d0)
            mm(2*i-1,2*i)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i-1)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i)=dcmplx(m2rr,0.d0)
         enddo
         call HEigensystem(6,mm,6,m,mix,6,1)
c one could simply order by increasing mass, as in the
c following lines (commented out), but for absurd conventional
c reasons we will order by flavor first and then by mass
c within each flavor
c         do i=1,6
c            if (m(i).lt.0.0d0) then
c               write (message,*) 'susfesct: negative ',pname(ksl(i)),
c     &              ' mass squared (',real(m(i)),')'
c               call dswrite(1,1,message)
c               unphys = -10-3-i
c            endif
c            mass(ksl(i)) = dsqrt(abs(m(i)))
c            do j=1,3
c               sldlmx(i,j) = mix(i,2*j-1)
c               sldrmx(i,j) = mix(i,2*j)
c            enddo
c         enddo
         do i=1,6
            do j=1,3
               c(i,j) = dsabsq(mix(i,2*j-1))+dsabsq(mix(i,2*j))
            enddo
         enddo
         do i=1,6
            imax(i)=0
         enddo
         do j=1,6
            cmax=0.d0
            do i=1,6
               aux=.true.
               if (j.gt.1) then
                  do k=1,j-1
                     aux=aux.and.i.ne.imax(k)
                  enddo
               endif
               if (aux.and.c(i,1+mod(j-1,3)).gt.cmax) then
                  cmax=c(i,1+mod(j-1,3))
                  imax(j)=i
               endif
            enddo
         enddo
         do i=1,6
            if (m(imax(i)).lt.0.0d0) then
               write (message,*) 'susfesct: negative ',pname(ksl(i)),
     &              ' mass squared (',real(m(imax(i))),')'
               call dswrite(1,1,message)
               unphys = -10-3-i
            endif
            mass(ksl(i)) = dsqrt(abs(m(imax(i))))
            do j=1,3
               sldlmx(i,j) = mix(imax(i),2*j-1)
               sldrmx(i,j) = mix(imax(i),2*j)
            enddo
         enddo
      endif

      if (msquarks.gt.0.d0) then
         do i=1,6
            mass(ksqu(i)) = msquarks
            mass(ksqd(i)) = msquarks
         enddo
         do i=1,3
            asoftu(i) = mu/tanbe
            asoftd(i) = mu*tanbe
         enddo
      else if (msquarks.lt.0.d0) then
         if (neuloop.ne.0) then
            call dswrite(0,1,
     &       'susfesct: neuloop set to zero since msfermions=mlsp etc')
            neuloop=0
         endif
         call dsneusct
         auxa = abs(mass(kn(kln)))
         if (auxa.lt.-msquarks) auxa = -msquarks
         do i=1,6
            mass(ksqu(i)) = auxa
            mass(ksqd(i)) = auxa
         enddo
         do i=1,3
            asoftu(i) = mu/tanbe
            asoftd(i) = mu*tanbe
         enddo
      else
c-------------------------------------------------------up-type squarks
         do i=1,6
            do j=1,6
               mm(i,j)=dcmplx(0.d0,0.d0)
            enddo
         enddo
         do i=1,3
            m2ll = mass2q(i)+mass(kqu(i))**2+d*(0.5d0-2.d0/3.d0*s2thw)
            m2rr = mass2u(i)+mass(kqu(i))**2+d*(2.d0/3.d0*s2thw)
            m2lr = mass(kqu(i))*(asoftu(i)-mu/tanbe)
            if (prtlevel.ge.2) then
               write (*,
     &            '('' squarku('',i1,'') mass^2 matrix = '',2(f12.2))')
     &              i,m2ll,m2lr
               write (*,
     &              '(''                            '',2(f12.2))')
     &              m2lr,m2rr
            endif
            mm(2*i-1,2*i-1)=dcmplx(m2ll,0.d0)
            mm(2*i-1,2*i)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i-1)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i)=dcmplx(m2rr,0.d0)
         enddo
         call HEigensystem(6,mm,6,m,mix,6,1)
c one could simply order by increasing mass, as in the
c following lines (commented out), but for absurd conventional
c reasons we will order by flavor first and then by mass
c within each flavor
c         do i=1,6
c            if (m(i).lt.0.0d0) then
c               write (message,*) 'susfesct: negative ',pname(ksqu(i)),
c     &              ' mass squared (',real(m(i)),')'
c               call dswrite(1,1,message)
c               unphys = -10-3-i
c            endif
c            mass(ksqu(i)) = dsqrt(abs(m(i)))
c            do j=1,3
c               squlmx(i,j) = mix(i,2*j-1)
c               squrmx(i,j) = mix(i,2*j)
c            enddo
c         enddo
         do i=1,6
            do j=1,3
               c(i,j) = dsabsq(mix(i,2*j-1))+dsabsq(mix(i,2*j))
            enddo
         enddo
         do i=1,6
            imax(i)=0
         enddo
         do j=1,6
            cmax=0.d0
            do i=1,6
               aux=.true.
               if (j.gt.1) then
                  do k=1,j-1
                     aux=aux.and.i.ne.imax(k)
                  enddo
               endif
               if (aux.and.c(i,1+mod(j-1,3)).gt.cmax) then
                  cmax=c(i,1+mod(j-1,3))
                  imax(j)=i
               endif
            enddo
         enddo
         do i=1,6
            if (m(imax(i)).lt.0.0d0) then
               write (message,*) 'susfesct: negative ',pname(ksqu(i)),
     &              ' mass squared (',real(m(imax(i))),')'
               call dswrite(1,1,message)
               unphys = -10-3-i
            endif
            mass(ksqu(i)) = dsqrt(abs(m(imax(i))))
            do j=1,3
               squlmx(i,j) = mix(imax(i),2*j-1)
               squrmx(i,j) = mix(imax(i),2*j)
            enddo
         enddo
c------------------------------------------------------down-type squarks
         do i=1,3
            m2ll = mass2q(i)+mass(kqd(i))**2+d*(-0.5d0+1.d0/3.d0*s2thw)
            m2rr = mass2d(i)+mass(kqd(i))**2+d*(-1.d0/3.d0*s2thw)
            m2lr = mass(kqd(i))*(asoftd(i)-mu*tanbe)
            if (prtlevel.ge.2) then
               write (*,
     &            '('' squarkd('',i1,'') mass^2 matrix = '',2(f12.2))')
     &              i,m2ll,m2lr
               write (*,
     &              '(''                            '',2(f12.2))')
     &              m2lr,m2rr
            endif
            mm(2*i-1,2*i-1)=dcmplx(m2ll,0.d0)
            mm(2*i-1,2*i)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i-1)=dcmplx(m2lr,0.d0)
            mm(2*i,2*i)=dcmplx(m2rr,0.d0)
         enddo
         call HEigensystem(6,mm,6,m,mix,6,1)
c one could simply order by increasing mass, as in the
c following lines (commented out), but for absurd conventional
c reasons we will order by flavor first and then by mass
c within each flavor
c         do i=1,6
c            if (m(i).lt.0.0d0) then
c               write (message,*) 'susfesct: negative ',pname(ksqd(i)),
c     &              ' mass squared (',real(m(i)),')'
c               call dswrite(1,1,message)
c               unphys = -10-3-i
c            endif
c            mass(ksqd(i)) = dsqrt(abs(m(i)))
c            do j=1,3
c               sqdlmx(i,j) = mix(i,2*j-1)
c               sqdrmx(i,j) = mix(i,2*j)
c            enddo
c         enddo
         do i=1,6
            do j=1,3
               c(i,j) = dsabsq(mix(i,2*j-1))+dsabsq(mix(i,2*j))
            enddo
         enddo
         do i=1,6
            imax(i)=0
         enddo
         do j=1,6
            cmax=0.d0
            do i=1,6
               aux=.true.
               if (j.gt.1) then
                  do k=1,j-1
                     aux=aux.and.i.ne.imax(k)
                  enddo
               endif
               if (aux.and.c(i,1+mod(j-1,3)).gt.cmax) then
                  cmax=c(i,1+mod(j-1,3))
                  imax(j)=i
               endif
            enddo
         enddo
         do i=1,6
            if (m(imax(i)).lt.0.0d0) then
               write (message,*) 'susfesct: negative ',pname(ksqd(i)),
     &              ' mass squared (',real(m(imax(i))),')'
               call dswrite(1,1,message)
               unphys = -10-3-i
            endif
            mass(ksqd(i)) = dsqrt(abs(m(imax(i))))
            do j=1,3
               sqdlmx(i,j) = mix(imax(i),2*j-1)
               sqdrmx(i,j) = mix(imax(i),2*j)
            enddo
         enddo
      endif

      if (unphys.lt.0) then
         write (message,*) 'unphysical sfermion sector (',unphys,')'
         call dswrite(1,1,message)
      endif

      end


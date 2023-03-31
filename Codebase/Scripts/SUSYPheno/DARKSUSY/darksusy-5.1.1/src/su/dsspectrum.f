      subroutine dsspectrum(unphys,hwarning)
c_______________________________________________________________________
c  particle spectrum and mixing matrices
c  uses sconst, sfesct, chasct, higsct, neusct
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1999
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsprep.h'

      integer i,unphys,hwarning
      real*8 zgm,dsabsq
      character*80 message
      integer j

c----------------------------------------------------------- gluino mass
      mass(kgluin) = abs(m3)

c------------------------------------------- sfermion masses and mixings
      call dssfesct(unphys)

c------------------------------------------ neutralino masses and mixing
      call dsneusct

c------------------------------------------- chargino masses and mixings
      call dschasct
      
c---------------------------------------- higgs boson masses and mixings
      call dshigsct(unphys,hwarning) ! must follow sfesct if higloop = 1

      if (unphys.lt.0) then
         write (message,*) 'unphysical higgs sector (',unphys,')'
         call dswrite(1,1,message)
         return
      endif

c-------------------------------------- lightest supersymmetric particle
      lsp = kn(kln)
      do i=21,48
         if (abs(mass(i)).lt.abs(mass(lsp))) lsp = i
      enddo
      if (lsp.lt.42.or.lsp.gt.45) then
         unphys=1
         write (message,*)
     &        'lsp is ',pname(lsp),
     &        ' m_lsp=',real(mass(lsp)),
     &        ' (m_x0=',real(mass(kn(kln))),')'
         call dswrite(1,1,message)
      endif
      zg = dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
      zgm = dsabsq(neunmx(1,3))+dsabsq(neunmx(1,4))
      lgzg = log10(zg)-log10(zgm)

      mx=mass(kn(1))

      end

************************************************************************
      subroutine dswwidth(unit)
c_______________________________________________________________________
c  write out a table of Higgs decay widths
c  input:
c    unit - logical unit to write on (integer)
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsmssm.h'
      integer unit,i,j,kh(4)
      character*3 hname(4)
      character*14 d0name(32),dpname(32)
      data hname/'H10','H20','A0 ','H+ '/
      data d0name/
     &  'H1 H1','H1 H2','H2 H2','H3 H3','H1 H3',
     &  'H2 H3','H+ H-','Z H1','Z H2','Z H3',
     &  'W+ H- & W- H+','ZZ','W+ W-','nu_e nu_e~','e+ e-',
     &  'nu_mu nu_mu~','mu+ mu-','nu_tau nu_tau~','tau+ tau-',
     &  'u u~',
     &  'd d~','c c~','s s~','t t~','b b~',
     &  'g g','q q g','ga ga','Z ga','sfermions',
     &  'neutralinos','charginos'/
      data dpname/
     &  'u d~','u s~','u b~','c d~','c s~',
     &  'c b~','t d~','t s~','t b~','nu_e e+',
     &  'nu_mu mu+','nu_tau tau+','W+ H1','W+ H2','W+ H3',
     &  4*' ',
     &  'sfermions','neu and cha',11*' '/

c...Set up Higgses
      kh(1)=kh1
      kh(2)=kh2
      kh(3)=kh3
      kh(4)=khc


      if (unit.le.0) return
      write(unit,*) '*** Higgs decay widths ***'
      do i=1,3
         write(unit,*) ' '
         write(unit,2000) hname(i),mass(kh(i)),width(kh(i))
         do j=1,32
           write(unit,1000) i,j,hname(i),d0name(j),hdwidth(j,i),
     &        hdwidth(j,i)/width(kh(i))
         enddo
      enddo

      do i=4,4
         write(unit,*) ' '
         write(unit,2000) hname(i),mass(kh(i)),width(kh(i))
         do j=1,21
           write(unit,1000) i,j,hname(i),dpname(j),hdwidth(j,i),
     &        hdwidth(j,i)/width(kh(i))
         enddo
      enddo


 1000 format (1x,I1,1x,I2,1x,A3,1x,'->',1x,A14,1x,E14.8,1x,F12.8)
 2000 format (1x,'Higgs: ',A3,1x,'Mass:',1x,E12.6,1x,'Width:',1x,E14.8)

      return
      end

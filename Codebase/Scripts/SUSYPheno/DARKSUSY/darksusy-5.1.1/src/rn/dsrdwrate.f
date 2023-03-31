************************************************************************
      subroutine dsrdwrate(unit1,unit2,ich)
c_______________________________________________________________________
c  write out a table of
c    initial cm momentum p
c    invariant annihilation rate w
c  input:
c    unit1 - logical unit to write total rate to (integer)
c    unit2 - logical unit to write partial differ. rates to (integer)
c    ich   - what initial channel to look at:
c              ich=1 nn-ann.  ich=2 cn-ann.  ich=3  cc-ann
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c  changes by je to include prtial in partials and coannihilation routines
c=======================================================================
      implicit none
      include 'dsandwcom.h'
      include 'dsmssm.h'
      real*8 du,umax,p,w,x,m,s,dsabsq
c     real*8 u,u2,y
      external dsandwdcosnn,dsandwdcoscn,dsandwdcoscc
      real*8 dsandwdcosnn,dsandwdcoscn,dsandwdcoscc
      integer unit1,unit2,i,n,c,ich,kchal,kstaul,chmax(5)
      character*12 chname(54,5)
c      logical doprtial
      parameter (n=10)!100)
      data chname /'h1 h1       ','h1 h2       ',
     &  'h2 h2       ','a a         ','h1 a        ',
     &  'h2 a        ','h+ h-       ','h1 z        ',
     &  'h2 z        ','a z         ','h+ w-       ',
     &  'z z         ','w w         ','nue nue     ',
     &  'e e         ','numu numu   ','mu mu       ',
     &  'nutau nutau ','tau tau     ','u u-bar     ',
     &  'd d-bar     ','c c-bar     ','s s-bar     ',
     &  't t-bar     ','b b-bar     ',29*'            ',
     &  'h+ h1       ','h+ h2       ','h+ a        ',
     &  'w+ h1       ','w+ h2       ','w+ a        ',
     &  'z h+        ','gamma h+    ','z w+        ',
     &  'gamma w+    ','nu_e e-bar  ','nu_e mu-bar ',
     &  'nu_e tau-bar',
     &  'nu_mu e-bar ','nu_mu mu-bar','nu_mu tau-ba',
     &  'nu_tau e-bar','nu_tau mu-ba','nu_tau tau-b',
     &  'u d-bar     ','u s-bar     ','u b-bar     ',
     &  'c d-bar     ','c s-bar     ','c b-bar     ',
     &  't d-bar     ','t s-bar     ','t b-bar     ',
     &  26*'            ',
     &  'h1 h1       ','h1 h2       ','h2 h2       ',
     &  'a a         ','h1 a        ','h2 a        ',
     &  'h+ h-       ','h1 z        ','h2 z        ',
     &  'a z         ','h+ w-       ','z z         ',
     &  'w+ w-       ','gamma gamma ','z0 gamma    ',
     &  'nu_e nu_e-ba',
     &  'nu_e nu_mu-b','nu_e nu_taub','e e-bar     ',
     &  'e mu-bar    ','e tau-bar   ','nu_mu nu_e-b',
     &  'nu_mu nu_m-b',
     &  'nu_mu nu_tab','mu e-bar    ','mu mu-bar   ',
     &  'mu tau-bar  ','nu_tau nu_eb','nu_tau nu_mb',
     &  'nu_tau nu_tb','tau e-bar   ','tau mu-bar  ',
     &  'tau tau-bar ','u u-bar     ','u c-bar     ',
     &  'u t-bar     ','d d-bar     ','d s-bar     ',
     &  'd b-bar     ','c u-bar     ','c c-bar     ',
     &  'c t-bar     ','s d-bar     ','s s-bar     ',
     &  's b-bar     ','t u-bar     ','t c-bar     ',
     &  't t-bar     ','b d-bar     ','b s-bar     ',
     &  'b b-bar     ',
     &  'h+ h+       ','h+ w+       ','w+ w+       ',
     &  'z tau       ','gamma tau   ','tau h       ',
     &  'tau H       ',50*'            ',
     &  'w+ w-       ','z z         ','z gamma     ',
     &  'gamma gamma ','z h         ','z H         ',
     &  'gamma h     ','gamma H     ','z a         ',
     &  'nu_e nu_e-ba','e e-bar     ','nu_mu nu_m-b',
     &  'mu mu-bar   ','nu_tau nu_tb',
     &  'tau tau-bar ','u u-bar     ','d d-bar     ',
     &  'c c-bar     ','s s-bar     ','t t-bar     ',
     &  'b b-bar     ','h h         ','h a         ',
     &  'H a         ','w+ h-       ','a a         ',
     &  'h H         ','H H         ','h+ h-       ',
     &  'tau tau     ',24*'            '/
      data chmax /25,28,54,4,30/

      if (unit1.le.0.and.unit2.le.0) return
      kchal=kcha(1)
      if (mass(kcha(2)).lt.mass(kcha(1))) kchal=kcha(2)
      kstaul=kstau(1)
      if (mass(kstau(2)).lt.mass(kstau(1))) kstaul=kstau(2)
      x=2.0d0
      m=abs(mass(kn(kln)))
      umax=5.0d0
      umax=90.0d0
      du=umax/n
      write (unit1,*)
      if (ich.eq.1) then
        write(unit1,*) 'neutralino-neutralino annihilation'
      elseif (ich.eq.2) then
        write(unit1,*) 'chargino-neutralino annihilation'
      elseif (ich.eq.3) then
        write(unit1,*) 'chargino-chargino annihilation'
      elseif (ich.eq.4) then
        write(unit1,*) 'stau-neutralino annihilation'
      elseif (ich.eq.5) then
        write(unit1,*) 'stau-stau annihilation'
      endif
      do 10 i=0,n
         p=i/real(n)*(m/10.)*200.0
         s=4*(m**2+p**2)
         if (ich.eq.1) then
           w=dsandwdcosnn(p,.1d0,kn(1),kn(1))
         elseif (ich.eq.2) then
           w=dsandwdcoscn(p,.1d0,kchal,kn(1))
         elseif (ich.eq.3) then
           w=dsandwdcoscc(p,.1d0,kchal,kchal)
         endif
         if (unit1.gt.0) then
           write (unit1,*)
           write(unit1,*) 'neutralino masses and gaugino fraction: '
           do c=1,4
             write(unit1,*) '    ',mass(kn(c)),'  ',
     &       dsabsq(neunmx(c,1)) + dsabsq(neunmx(c,2))
           enddo
           write(unit1,*) 'chargino masses: '
           do c=1,2
             write(unit1,*) '    ',mass(kcha(c)),'  '
           enddo
           write(unit1,*) 'stau masses: '
           do c=1,2
             write(unit1,*) '    ',mass(kstau(c)),'  '
           enddo
           write (unit1,*) ' p = ',p
         endif
         write (unit2,*) ' sqrt(s) = ',sqrt(s)
         write (unit2,*) ' dsandwdcos = ',w
         do c=1,chmax(ich)
           write(unit2,*) chname(c,ich),prtial(c)
         enddo
 10   continue
 1000 format (4(2x,e12.6))
      end

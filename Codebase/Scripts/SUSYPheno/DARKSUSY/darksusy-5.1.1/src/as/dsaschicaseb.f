c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.34, October 8, 2001, edsjo@physto.se)
c....Template file for dsaschicaseb begins here

**************************************************************
*** SUBROUTINE dsaschicaseb                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** anti-sfermion(i) + neutralino(j)/chargino^+(j)         *** 
*** -> gauge-boson + anti-fermion                          ***
***                                                        ***
*** The anti-sfermion must be the first mentioned          ***
*** particle (kp1) and the neutralino/chargino             *** 
*** the other (kp2) -- not the opposite.                   ***
*** For the final state the gauge boson must be mentioned  ***
*** first (i.e. kp3) and next the anti-fermion (kp4) --    ***
*** not the opposite.                                      ***
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-03                                         ***
*** QCD included: 02-03-21                                 ***
*** comment added by Piero Ullio, 02-07-01                 ***
*** added flavour changing charged exchange:               ***
*** for the 1st,2nd, 3rd and 4th case                      ***
*** added by Mia Schelke January 2007                      ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

***** Note that it is assumed that coupling constants that do 
***** not exist have already been set to zero!!!!!
***** Thus, many of the coefficients defined in this code 
***** simplify when the diagrams contain sneutrinos 
***** or neutrinos.

      subroutine dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro

      complex*16 dsaschicasebas,tmp1,tmp2
      real*8  msfer2,mchi2,mgb2,mfer2
      real*8  massorga
      complex*16  zfbtt1l,zfbtt1r,zfbtt2l,zfbtt2r
      complex*16  zfbuu1l,zfbuu1r,zfbuu2l,zfbuu2r
      complex*16  zfbuu3l,zfbuu3r,zfbuu4l,zfbuu4r,zfbuu5l,zfbuu5r
      complex*16  zfbuu6l,zfbuu6r,zfbuu7l,zfbuu7r,zfbuu8l,zfbuu8r
      complex*16  zfbss,zfbss1l,zfbss1r,zfbss2l,zfbss2r 
      complex*16  zfbss3l,zfbss3r,zfbss4l,zfbss4r,zfbss5l,zfbss5r
      complex*16  zfbss6l,zfbss6r,zfbss7l,zfbss7r,zfbss8l,zfbss8r
      complex*16  zfbtu1l,zfbtu1r,zfbtu2l,zfbtu2r
      complex*16  zfbtu3l,zfbtu3r,zfbtu4l,zfbtu4r
      complex*16  zfbts,zfbts1l,zfbts1r,zfbts2l,zfbts2r
      complex*16  zfbts3l,zfbts3r,zfbts4l,zfbts4r
      complex*16  zfbsu,zfbsu1l,zfbsu1r,zfbsu2l,zfbsu2r 
      complex*16  zfbsu3l,zfbsu3r,zfbsu4l,zfbsu4r,zfbsu5l,zfbsu5r
      complex*16  zfbsu6l,zfbsu6r,zfbsu7l,zfbsu7r,zfbsu8l,zfbsu8r
      real*8      ctt,cuu,css,ctu,cts,csu  
      integer i,k,l
      integer kgb,kfer
      integer ksfert(6),kchiu(4)
      integer kfers(3)
      integer nsfert,nchiu,nfers
      real*8 pi
      parameter (pi=3.141592653589793238d0)

c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      

***** set particles in final state:
      kgb=kp3
      kfer=kp4

***** masses in final state:  
      mass3=mass(kgb)  
      mass4=mass(kfer)  
***** define the kinematic variables  
      call dsaskinset2
   
c....symmetry factor, for non-identical final state particles
      s34=1.d0

c....mass symbols used in form-code  
      msfer2=mass(kp1)**2
      mchi2=mass(kp2)**2
      mgb2=mass(kgb)**2
      mfer2=mass(kfer)**2
      s=Svar  
      t=Tvar  
      u=Uvar  

***** now set the color factors to 1.d0 for a slepton in
***** the initial state and to 3.d0 for a squark in 
***** the initial state
***** if the gauge boson in the final state is a gluon, then the 
***** color factor is changed -- see case six below 
***** the color factors enter the form generated expression for 
***** the squared amplitude

      ctt=1.d0
      cuu=1.d0
      css=1.d0
      ctu=1.d0
      cts=1.d0
      csu=1.d0

      if(iifam(1).ge.ivfam(ksu1).and.iifam(1).le.ivfam(ksb2)) then
        ctt=3.d0
        cuu=3.d0
        css=3.d0
        ctu=3.d0
        cts=3.d0
        csu=3.d0
      endif

***** initially setting array of exchange empty
      do i=1,6
         ksfert(i)=0
      enddo
      do i=1,3
         kfers(i)=0
      enddo

*****
*****
***** the first case    
***** anti-up-type-sfermion(i) + \chi^+(j) 
***** -> Z + anti-down-type-fermion
      if(icase.eq.1) then
        massorga=1.d0/mass(kgb)**2
***** now set particles in the intermediate states 
***** in the t-channel (1 or 2 sfermions):  
        nsfert=ncsfert
        do i=1,nsfert
          ksfert(i)=kcsfertn(i)
        enddo
***** in the u-channel (2 charginos):
        do k=1,2
          kchiu(k)=kcha(k)
          kchiu(k+2)=0
        enddo
        nchiu=2    
***** in the s-channel (1 fermion):
        kfers(1)=kfer
        nfers=1
        goto 200
      endif
***** 
***** the second case
***** anti-up-type-sfermion(i) + neutralino(j)
***** -> W^- + anti-down-type-fermion
      if(icase.eq.2) then
        massorga=1.d0/mass(kgb)**2 
***** set particles in the intermediate states 
***** in the t-channel (2 sfermions):
        nsfert=ncsfertc
        do i=1,nsfert ! JE Correction 2009-04-03
          ksfert(i)=kcsfertc(i)
        enddo
***** in the u-channel (2 charginos):
        do k=1,2
          kchiu(k)=kcha(k)
          kchiu(k+2)=0
        enddo
        nchiu=2
***** in the s-channel (1 fermion): 
        kfers(1)=kcfers
        nfers=1
        goto 200
      endif
*****
***** the third case
***** anti-up-type-sfermion(i) + chargino^+(j) 
***** -> W^+ + anti-up-type-fermion
      if(icase.eq.3) then
        massorga=1.d0/mass(kgb)**2
***** no t-channel:  
        nsfert=0
***** in the u-channel (0 or 4 neutralinos):
        if(iifam(1).eq.ivfam(kfer)) then
          do k=1,4
            kchiu(k)=kn(k)
          enddo
          nchiu=4
        else
          nchiu=0
        endif
***** in the s-channel (1 lepton or 3 quarks):
        nfers=ncferd
        do i=1,nfers
          kfers(i)=kcferd(i)
        enddo
        goto 200
      endif
*****
***** the fourth case
***** anti-down-type-sfermion(i) + chargino^+(j) 
***** -> W^+ + anti-down-type-fermion
      if(icase.eq.4) then
        massorga=1.d0/mass(kgb)**2
***** in the t-channel (1 or 6 sfermions):
        nsfert=ncsfertc
        do i=1,nsfert
          ksfert(i)=kcsfertc(i)
        enddo
***** in the u-channel (0 or 4 neutralinos):
        if(iifam(1).eq.ivfam(kfer)) then
          do k=1,4
            kchiu(k)=kn(k)
          enddo
          nchiu=4
        else
          nchiu=0
        endif
***** no s-channel:
        nfers=0
        goto 200
      endif
      
      write(*,*) 'DS: dsaschicaseb called with wrong icase : ',icase   
      write(*,*) 'DS: initial or final states : '  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  
      
 200  continue

c....specify the coefficients used in the form-code
c....
c....note that when writing down the Feynman rules for the diagrams 
c....we have started with the anti-fermion instead of ending there
c....(the two results are identical, which can be seen by 
c....transposing one of them)
c....in this way we got a result that looked very similar to the 
c....one in dsaschicasea
c....
c....first for the t-channel(sfermion exchange)
      zfbtt1l=dcmplx(0.d0,0.d0)
      zfbtt1r=dcmplx(0.d0,0.d0)
      zfbtt2l=dcmplx(0.d0,0.d0)
      zfbtt2r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
      do k=1,nsfert
        do l=1,nsfert 
         zfbtt1l=zfbtt1l
     &        +dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,kp1,ksfert(k))*conjg(gl(kgb,kp1,ksfert(l)))
     &        *gr(ksfert(k),kp2,kfer)*conjg(gr(ksfert(l),kp2,kfer))
         zfbtt1r=zfbtt1r
     &        +dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,kp1,ksfert(k))*conjg(gl(kgb,kp1,ksfert(l)))
     &        *gl(ksfert(k),kp2,kfer)*conjg(gl(ksfert(l),kp2,kfer))
         zfbtt2l=zfbtt2l+mass(kfer)*mass(kp2)
     &        *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,kp1,ksfert(k))*conjg(gl(kgb,kp1,ksfert(l)))
     &        *gl(ksfert(k),kp2,kfer)*conjg(gr(ksfert(l),kp2,kfer))
         zfbtt2r=zfbtt2r+mass(kfer)*mass(kp2)
     &        *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &        *gl(kgb,kp1,ksfert(k))*conjg(gl(kgb,kp1,ksfert(l)))
     &        *gr(ksfert(k),kp2,kfer)*conjg(gl(ksfert(l),kp2,kfer))
        enddo
      enddo
      endif


c....then the u-channel(chargino/neutralino exchange)
      zfbuu1l=dcmplx(0.d0,0.d0)
      zfbuu1r=dcmplx(0.d0,0.d0)
      zfbuu2l=dcmplx(0.d0,0.d0)
      zfbuu2r=dcmplx(0.d0,0.d0)
      zfbuu3l=dcmplx(0.d0,0.d0)
      zfbuu3r=dcmplx(0.d0,0.d0)
      zfbuu4l=dcmplx(0.d0,0.d0)
      zfbuu4r=dcmplx(0.d0,0.d0)
      zfbuu5l=dcmplx(0.d0,0.d0)
      zfbuu5r=dcmplx(0.d0,0.d0)
      zfbuu6l=dcmplx(0.d0,0.d0)
      zfbuu6r=dcmplx(0.d0,0.d0)
      zfbuu7l=dcmplx(0.d0,0.d0)
      zfbuu7r=dcmplx(0.d0,0.d0)
      zfbuu8l=dcmplx(0.d0,0.d0)
      zfbuu8r=dcmplx(0.d0,0.d0)
      if(nchiu.gt.0) then
      do k=1,nchiu
        do l=1,nchiu  
         zfbuu1l=zfbuu1l
     &    +dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu1r=zfbuu1r
     &    +dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu2l=zfbuu2l+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu2r=zfbuu2r+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu3l=zfbuu3l+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu3r=zfbuu3r+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu4l=zfbuu4l+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu4r=zfbuu4r+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu5l=zfbuu5l+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu5r=zfbuu5r+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu6l=zfbuu6l+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         zfbuu6r=zfbuu6r+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu7l=zfbuu7l+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu7r=zfbuu7r+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbuu8l=zfbuu8l+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(kgb,kchiu(k),kp2)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbuu8r=zfbuu8r+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(kgb,kchiu(k),kp2)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
        enddo
      enddo
      endif

c....then the s-channel(fermion-exchange)
      zfbss=dcmplx(0.d0,0.d0)
      zfbss1l=dcmplx(0.d0,0.d0)
      zfbss1r=dcmplx(0.d0,0.d0)
      zfbss2l=dcmplx(0.d0,0.d0)
      zfbss2r=dcmplx(0.d0,0.d0) 
      zfbss3l=dcmplx(0.d0,0.d0)
      zfbss3r=dcmplx(0.d0,0.d0)
      zfbss4l=dcmplx(0.d0,0.d0)
      zfbss4r=dcmplx(0.d0,0.d0) 
      zfbss5l=dcmplx(0.d0,0.d0)
      zfbss5r=dcmplx(0.d0,0.d0)
      zfbss6l=dcmplx(0.d0,0.d0)
      zfbss6r=dcmplx(0.d0,0.d0) 
      zfbss7l=dcmplx(0.d0,0.d0)
      zfbss7r=dcmplx(0.d0,0.d0)
      zfbss8l=dcmplx(0.d0,0.d0)
      zfbss8r=dcmplx(0.d0,0.d0)
      if(nfers.gt.0) then
***** the expression below has ben changed Jan 2007
***** before we had only kfers now kfers(l),kfers(k)
      do k=1,nfers
        do l=1,nfers
         zfbss=dsasdepro(s,kfers(k))*conjg(dsasdepro(s,kfers(l)))
         zfbss1l=zfbss1l+zfbss*gr(kgb,kfers(k),kfer)
     &      *gr(kp1,kp2,kfers(k))*conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer)) 
         zfbss1r=zfbss1r+zfbss*gl(kgb,kfers(k),kfer)
     &      *gl(kp1,kp2,kfers(k))*conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer)) 
         zfbss2l=zfbss2l+zfbss*mass(kp2)*mass(kfers(l))
     &      *gr(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer)) 
         zfbss2r=zfbss2r+zfbss*mass(kp2)*mass(kfers(l))
     &      *gl(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer))  
         zfbss3l=zfbss3l+zfbss*mass(kfers(k))*mass(kfers(l))
     &      *gr(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer))  
         zfbss3r=zfbss3r+zfbss*mass(kfers(k))*mass(kfers(l))
     &      *gl(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer)) 
         zfbss4l=zfbss4l+zfbss*mass(kfers(k))*mass(kp2)
     &      *gr(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer)) 
         zfbss4r=zfbss4r+zfbss*mass(kfers(k))*mass(kp2)
     &      *gl(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer))  
         zfbss5l=zfbss5l+zfbss*mass(kfer)*mass(kfers(k))
     &      *mass(kp2)*mass(kfers(l))
     &      *gl(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer))  
         zfbss5r=zfbss5r+zfbss*mass(kfer)*mass(kfers(k))
     &      *mass(kp2)*mass(kfers(l))
     &      *gr(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer))
         zfbss6l=zfbss6l+zfbss*mass(kfer)*mass(kfers(k))
     &      *gl(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer))
         zfbss6r=zfbss6r+zfbss*mass(kfer)*mass(kfers(k))
     &      *gr(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer)) 
         zfbss7l=zfbss7l+zfbss*mass(kfer)*mass(kp2)
     &      *gl(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer)) 
         zfbss7r=zfbss7r+zfbss*mass(kfer)*mass(kp2)
     &      *gr(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer)) 
         zfbss8l=zfbss8l+zfbss*mass(kfer)*mass(kfers(l))
     &      *gl(kgb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))
     &      *conjg(gr(kgb,kfers(l),kfer)) 
         zfbss8r=zfbss8r+zfbss*mass(kfer)*mass(kfers(l))
     &      *gr(kgb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))
     &      *conjg(gl(kgb,kfers(l),kfer))
        enddo
      enddo
***** we must redefine the term zfbss, because it will be 
***** multiplied on every ss term in the form expression 
***** this was ok before when just one exchange particle
***** but not now with 3 exchange particles (k,l)
***** where it has to be introduced on the terms here above
      zfbss=dcmplx(1.d0,0.d0)
      endif

c....then t-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      zfbtu1l=dcmplx(0.d0,0.d0)
      zfbtu1r=dcmplx(0.d0,0.d0)
      zfbtu2l=dcmplx(0.d0,0.d0)
      zfbtu2r=dcmplx(0.d0,0.d0)
      zfbtu3l=dcmplx(0.d0,0.d0)
      zfbtu3r=dcmplx(0.d0,0.d0)
      zfbtu4l=dcmplx(0.d0,0.d0)
      zfbtu4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nchiu.gt.0) then
      do k=1,nsfert
        do l=1,nchiu
         zfbtu1l=zfbtu1l
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kgb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         zfbtu1r=zfbtu1r
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kgb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbtu2l=zfbtu2l
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kp2)*mass(kchiu(l))
     &    *gl(kgb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         zfbtu2r=zfbtu2r
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kp2)*mass(kchiu(l))
     &    *gl(kgb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbtu3l=zfbtu3l
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kchiu(l))
     &    *gl(kgb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         zfbtu3r=zfbtu3r
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kchiu(l))
     &    *gl(kgb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         zfbtu4l=zfbtu4l
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)
     &    *gl(kgb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gr(kgb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         zfbtu4r=zfbtu4r
     &    +dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)
     &    *gl(kgb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gl(kgb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
        enddo
      enddo
      endif
  

c....s-channel amplitude multiplied by
c....hermitian conjugated t-channel amplitude 
      zfbts=dcmplx(0.d0,0.d0)
      zfbts1l=dcmplx(0.d0,0.d0)
      zfbts1r=dcmplx(0.d0,0.d0)
      zfbts2l=dcmplx(0.d0,0.d0)
      zfbts2r=dcmplx(0.d0,0.d0)
      zfbts3l=dcmplx(0.d0,0.d0)
      zfbts3r=dcmplx(0.d0,0.d0)
      zfbts4l=dcmplx(0.d0,0.d0)
      zfbts4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nfers.gt.0) then
*       zfbts=dsasdepro(s,kfers)
***** changed Jan 2007 to have kfers(l) not just kfers
***** must then include zfbts in all terms from the start
***** instaed of in the form result, therefore  
***** set zfbts=1 at the end
       do k=1,nsfert
         do l=1,nfers
         zfbts=dsasdepro(s,kfers(l))         
         zfbts1l=zfbts1l+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gr(kgb,kfers(l),kfer)
     &       *gr(kp1,kp2,kfers(l))*conjg(gr(ksfert(k),kp2,kfer))
         zfbts1r=zfbts1r+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gl(kgb,kfers(l),kfer)
     &       *gl(kp1,kp2,kfers(l))*conjg(gl(ksfert(k),kp2,kfer))
         zfbts2l=zfbts2l+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfers(l))*mass(kp2)
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gr(kgb,kfers(l),kfer)
     &       *gl(kp1,kp2,kfers(l))*conjg(gr(ksfert(k),kp2,kfer))
         zfbts2r=zfbts2r+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfers(l))*mass(kp2)
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gl(kgb,kfers(l),kfer)
     &       *gr(kp1,kp2,kfers(l))*conjg(gl(ksfert(k),kp2,kfer))
         zfbts3l=zfbts3l+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfer)*mass(kp2)
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gl(kgb,kfers(l),kfer)
     &       *gl(kp1,kp2,kfers(l))*conjg(gr(ksfert(k),kp2,kfer))
         zfbts3r=zfbts3r+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfer)*mass(kp2)
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gr(kgb,kfers(l),kfer)
     &       *gr(kp1,kp2,kfers(l))*conjg(gl(ksfert(k),kp2,kfer))
         zfbts4l=zfbts4l+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfer)*mass(kfers(l))
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gl(kgb,kfers(l),kfer)
     &       *gr(kp1,kp2,kfers(l))*conjg(gr(ksfert(k),kp2,kfer))
         zfbts4r=zfbts4r+zfbts*conjg(dsasdepro(t,ksfert(k)))
     &       *mass(kfer)*mass(kfers(l))
     &       *conjg(gl(kgb,kp1,ksfert(k)))*gr(kgb,kfers(l),kfer)
     &       *gl(kp1,kp2,kfers(l))*conjg(gl(ksfert(k),kp2,kfer))
         enddo
       enddo
       zfbts=dcmplx(1.d0,0.d0)
      endif
 

c....s-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      zfbsu=dcmplx(0.d0,0.d0) 
      zfbsu1l=dcmplx(0.d0,0.d0)
      zfbsu1r=dcmplx(0.d0,0.d0)
      zfbsu2l=dcmplx(0.d0,0.d0)
      zfbsu2r=dcmplx(0.d0,0.d0) 
      zfbsu3l=dcmplx(0.d0,0.d0)
      zfbsu3r=dcmplx(0.d0,0.d0)
      zfbsu4l=dcmplx(0.d0,0.d0)
      zfbsu4r=dcmplx(0.d0,0.d0)
      zfbsu5l=dcmplx(0.d0,0.d0)
      zfbsu5r=dcmplx(0.d0,0.d0)
      zfbsu6l=dcmplx(0.d0,0.d0)
      zfbsu6r=dcmplx(0.d0,0.d0)
      zfbsu7l=dcmplx(0.d0,0.d0)
      zfbsu7r=dcmplx(0.d0,0.d0)
      zfbsu8l=dcmplx(0.d0,0.d0)
      zfbsu8r=dcmplx(0.d0,0.d0) 
      if(nfers.gt.0.and.nchiu.gt.0) then
***** changed Jan 2007, comments see st interference terms
*       zfbsu=dsasdepro(s,kfers)
       do k=1,nchiu 
         do l=1,nfers
          zfbsu=dsasdepro(s,kfers(l))
          zfbsu1l=zfbsu1l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *gr(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu1r=zfbsu1r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *gl(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu2l=zfbsu2l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gr(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu2r=zfbsu2r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gl(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer)) 
          zfbsu3l=zfbsu3l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kchiu(k))
     &           *gr(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu3r=zfbsu3r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kchiu(k))
     &           *gl(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu4l=zfbsu4l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kp2)
     &           *gr(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu4r=zfbsu4r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfers(l))*mass(kp2)
     &           *gl(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu5l=zfbsu5l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gl(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu5r=zfbsu5r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *mass(kp2)*mass(kchiu(k))
     &           *gr(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu6l=zfbsu6l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *gl(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu6r=zfbsu6r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kfers(l))
     &           *gr(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu7l=zfbsu7l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kp2)
     &           *gl(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu7r=zfbsu7r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kp2)
     &           *gr(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))
          zfbsu8l=zfbsu8l+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kchiu(k))
     &           *gl(kgb,kfers(l),kfer)*gl(kp1,kp2,kfers(l))
     &           *conjg(gl(kgb,kchiu(k),kp2))
     &           *conjg(gr(kp1,kchiu(k),kfer))
          zfbsu8r=zfbsu8r+zfbsu*conjg(dsasdepro(u,kchiu(k)))
     &           *mass(kfer)*mass(kchiu(k))
     &           *gr(kgb,kfers(l),kfer)*gr(kp1,kp2,kfers(l))
     &           *conjg(gr(kgb,kchiu(k),kp2))
     &           *conjg(gl(kp1,kchiu(k),kfer))  
         enddo  
       enddo
       zfbsu=dcmplx(1.d0,0.d0)
      endif

***** After the template file follows the form expression of the 
***** ('complex') amplitude squared: dsaschicasebas
***** which is (M_tM_t^\dagger + M_sM_s^\dagger + M_uM_u^\dagger
***** + 2M_tM_u^\dagger + 2M_sM_t^\dagger + 2M_sM_u^\dagger)
***** As M_tM_u^\dagger+M_uM_t^\dagger=2Re(M_tM_u^\dagger) 
***** we have to take the real part of the interference terms,
***** in order to obtain the true amplitude squared.
***** This is most easily done by taking the real part of the
***** whole expression dsaschicasebas.
***** This is done at the very end of the code.

c....Template file for dsaschicaseb ends here









      tmp1=
     &  + msfer2*mchi2*mgb2*massorga*zfbuu1l*cuu
     &  + msfer2*mchi2*mgb2*massorga*zfbuu1r*cuu
     &  + msfer2*mchi2*mfer2*massorga*zfbuu1l*cuu
     &  + msfer2*mchi2*mfer2*massorga*zfbuu1r*cuu
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfbts*zfbts1l*cts
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfbts*zfbts1r*cts
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*msfer2*mchi2*mfer2*massorga*zfbsu*zfbsu1r*csu
     &  - msfer2*mchi2*s*massorga*zfbuu1l*cuu
     &  - msfer2*mchi2*s*massorga*zfbuu1r*cuu
     &  + 2.d0*msfer2*mchi2*s*massorga*zfbts*zfbts1l*cts
     &  + 2.d0*msfer2*mchi2*s*massorga*zfbts*zfbts1r*cts
     &  + 2.d0*msfer2*mchi2*s*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*msfer2*mchi2*s*massorga*zfbsu*zfbsu1r*csu
     &  - 2.d0*msfer2*mchi2*t*massorga*zfbtt1l*ctt
     &  - 2.d0*msfer2*mchi2*t*massorga*zfbtt1r*ctt
     &  - 2.d0*msfer2*mchi2*t*massorga*zfbuu1l*cuu
     &  - 2.d0*msfer2*mchi2*t*massorga*zfbuu1r*cuu
     &  - 4.d0*msfer2*mchi2*t*massorga*zfbtu1l*ctu
     &  - 4.d0*msfer2*mchi2*t*massorga*zfbtu1r*ctu
     &  - 2.d0*msfer2*mchi2*zfbtt1l*ctt
     &  - 2.d0*msfer2*mchi2*zfbtt1r*ctt
     &  - 2.d0*msfer2*mchi2*zfbtu1l*ctu
     &  - 2.d0*msfer2*mchi2*zfbtu1r*ctu
     &  - 4.d0*msfer2*mchi2*zfbts*zfbts1l*cts
     &  - 4.d0*msfer2*mchi2*zfbts*zfbts1r*cts
     &  - 4.d0*msfer2*mchi2*zfbsu*zfbsu1l*csu
     &  - 4.d0*msfer2*mchi2*zfbsu*zfbsu1r*csu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbuu1l*cuu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbuu1r*cuu
     &  + msfer2*mgb2*mfer2*massorga*zfbss*zfbss1l*css
     &  + msfer2*mgb2*mfer2*massorga*zfbss*zfbss1r*css
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbtu1l*ctu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbtu1r*ctu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*msfer2*mgb2*mfer2*massorga*zfbsu*zfbsu1r*csu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfbuu1l*cuu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfbuu1r*cuu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfbtu1l*ctu
     &  - 2.d0*msfer2*mgb2*t*massorga*zfbtu1r*ctu
     &  + msfer2*mgb2*massorga*zfbuu3l*cuu
     &  + msfer2*mgb2*massorga*zfbuu3r*cuu
     &  + msfer2*mgb2*massorga*zfbuu6l*cuu
     &  + msfer2*mgb2*massorga*zfbuu6r*cuu
     &  + 2.d0*msfer2*mgb2*massorga*zfbuu7l*cuu
     &  + 2.d0*msfer2*mgb2*massorga*zfbuu7r*cuu
     &  + msfer2*mgb2*massorga*zfbuu8l*cuu
     &  + msfer2*mgb2*massorga*zfbuu8r*cuu
     &  + msfer2*mgb2*massorga*zfbss*zfbss3l*css
     &  + msfer2*mgb2*massorga*zfbss*zfbss3r*css
     &  - msfer2*mgb2*massorga*zfbss*zfbss6l*css
     &  - msfer2*mgb2*massorga*zfbss*zfbss6r*css
     &  - msfer2*mgb2*massorga*zfbss*zfbss8l*css
     &  - msfer2*mgb2*massorga*zfbss*zfbss8r*css
     &  + 2.d0*msfer2*mgb2*massorga*zfbtu2l*ctu
     &  + 2.d0*msfer2*mgb2*massorga*zfbtu2r*ctu
     &  + 4.d0*msfer2*mgb2*massorga*zfbtu4l*ctu
     &  + 4.d0*msfer2*mgb2*massorga*zfbtu4r*ctu
     &  - 2.d0*msfer2*mgb2*massorga*zfbts*zfbts2l*cts
     &  - 2.d0*msfer2*mgb2*massorga*zfbts*zfbts2r*cts
     &  + 2.d0*msfer2*mgb2*massorga*zfbts*zfbts3l*cts
     &  + 2.d0*msfer2*mgb2*massorga*zfbts*zfbts3r*cts
     &  - 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu3l*csu
     &  - 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu3r*csu
     &  - 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu8l*csu
     &  + 2.d0*msfer2*mgb2*massorga*zfbsu*zfbsu8r*csu
     &  + 2.d0*msfer2*mgb2*zfbuu1l*cuu
     &  + 2.d0*msfer2*mgb2*zfbuu1r*cuu
     &  + 2.d0*msfer2*mgb2*zfbss*zfbss1l*css
     &  + 2.d0*msfer2*mgb2*zfbss*zfbss1r*css
     &  - 4.d0*msfer2*mgb2*zfbsu*zfbsu1l*csu
     &  - 4.d0*msfer2*mgb2*zfbsu*zfbsu1r*csu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbuu1l*cuu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbuu1r*cuu
     &  + msfer2*mfer2*s*massorga*zfbss*zfbss1l*css
     &  + msfer2*mfer2*s*massorga*zfbss*zfbss1r*css
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbtu1l*ctu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbtu1r*ctu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*msfer2*mfer2*s*massorga*zfbsu*zfbsu1r*csu
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbtt1l*ctt
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbtt1r*ctt
     &  - 4.d0*msfer2*mfer2*t*massorga*zfbuu1l*cuu
     &  - 4.d0*msfer2*mfer2*t*massorga*zfbuu1r*cuu
     &  - 6.d0*msfer2*mfer2*t*massorga*zfbtu1l*ctu
     &  - 6.d0*msfer2*mfer2*t*massorga*zfbtu1r*ctu
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbts*zfbts1l*cts
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbts*zfbts1r*cts
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*msfer2*mfer2*t*massorga*zfbsu*zfbsu1r*csu
     &  + msfer2*mfer2*massorga*zfbuu3l*cuu
     &  + msfer2*mfer2*massorga*zfbuu3r*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfbuu6l*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfbuu6r*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfbuu8l*cuu
     &  + 2.d0*msfer2*mfer2*massorga*zfbuu8r*cuu
     &  + msfer2*mfer2*massorga*zfbss*zfbss3l*css
     &  + msfer2*mfer2*massorga*zfbss*zfbss3r*css
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu2l*ctu
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu2r*ctu
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu3l*ctu
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu3r*ctu
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu4l*ctu
     &  + 2.d0*msfer2*mfer2*massorga*zfbtu4r*ctu
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts2l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts2r*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts3l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts3r*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts4l*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbts*zfbts4r*cts
     &  - 2.d0*msfer2*mfer2*massorga*zfbsu*zfbsu3l*csu
     &  - 2.d0*msfer2*mfer2*massorga*zfbsu*zfbsu3r*csu
     &  - 4.d0*msfer2*mfer2*massorga*zfbsu*zfbsu6l*csu
     &  - 4.d0*msfer2*mfer2*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*msfer2*mfer2*massorga*zfbsu*zfbsu7l*csu
     &  - 2.d0*msfer2*mfer2*massorga*zfbsu*zfbsu7r*csu
     &  - 2.d0*msfer2*mfer2*zfbtt1l*ctt
     &  - 2.d0*msfer2*mfer2*zfbtt1r*ctt
     &  + 2.d0*msfer2*mfer2*zfbuu1l*cuu
     &  + 2.d0*msfer2*mfer2*zfbuu1r*cuu
     &  - 2.d0*msfer2*mfer2*zfbss*zfbss1l*css
     &  - 2.d0*msfer2*mfer2*zfbss*zfbss1r*css
     &  - 6.d0*msfer2*mfer2*zfbtu1l*ctu
     &  - 6.d0*msfer2*mfer2*zfbtu1r*ctu
     &  - 2.d0*msfer2*mfer2*zfbts*zfbts1l*cts
     &  - 2.d0*msfer2*mfer2*zfbts*zfbts1r*cts
     &  - 8.d0*msfer2*mfer2*zfbsu*zfbsu1l*csu
     &  - 8.d0*msfer2*mfer2*zfbsu*zfbsu1r*csu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbuu1l*cuu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbuu1r*cuu
     &  - msfer2*mfer2**2*massorga*zfbss*zfbss1l*css
     &  - msfer2*mfer2**2*massorga*zfbss*zfbss1r*css
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbtu1l*ctu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbtu1r*ctu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*msfer2*mfer2**2*massorga*zfbsu*zfbsu1r*csu
     &  + 2.d0*msfer2*s*t*massorga*zfbuu1l*cuu
     &  + 2.d0*msfer2*s*t*massorga*zfbuu1r*cuu
     &  + 2.d0*msfer2*s*t*massorga*zfbtu1l*ctu
     &  + 2.d0*msfer2*s*t*massorga*zfbtu1r*ctu
     &  - 2.d0*msfer2*s*t*massorga*zfbts*zfbts1l*cts
     &  - 2.d0*msfer2*s*t*massorga*zfbts*zfbts1r*cts
     &  - 2.d0*msfer2*s*t*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*msfer2*s*t*massorga*zfbsu*zfbsu1r*csu
     &  - msfer2*s*massorga*zfbuu3l*cuu
     &  - msfer2*s*massorga*zfbuu3r*cuu
     &  - 2.d0*msfer2*s*massorga*zfbuu6l*cuu
     &  - 2.d0*msfer2*s*massorga*zfbuu6r*cuu
     &  - 2.d0*msfer2*s*massorga*zfbuu8l*cuu
     &  - 2.d0*msfer2*s*massorga*zfbuu8r*cuu
     &  - msfer2*s*massorga*zfbss*zfbss3l*css
     &  - msfer2*s*massorga*zfbss*zfbss3r*css
     &  - 2.d0*msfer2*s*massorga*zfbtu2l*ctu
     &  - 2.d0*msfer2*s*massorga*zfbtu2r*ctu
     &  - 2.d0*msfer2*s*massorga*zfbtu3l*ctu
     &  - 2.d0*msfer2*s*massorga*zfbtu3r*ctu
     &  - 2.d0*msfer2*s*massorga*zfbtu4l*ctu
     &  - 2.d0*msfer2*s*massorga*zfbtu4r*ctu
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts2l*cts
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts2r*cts
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts3l*cts
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts3r*cts
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts4l*cts
     &  + 2.d0*msfer2*s*massorga*zfbts*zfbts4r*cts
     &  + 2.d0*msfer2*s*massorga*zfbsu*zfbsu3l*csu
     &  + 2.d0*msfer2*s*massorga*zfbsu*zfbsu3r*csu
     &  + 4.d0*msfer2*s*massorga*zfbsu*zfbsu6l*csu
     &  + 4.d0*msfer2*s*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*msfer2*s*massorga*zfbsu*zfbsu7l*csu
     &  + 2.d0*msfer2*s*massorga*zfbsu*zfbsu7r*csu
     &  - 2.d0*msfer2*s*zfbuu1l*cuu
     &  - 2.d0*msfer2*s*zfbuu1r*cuu
     &  - 2.d0*msfer2*s*zfbss*zfbss1l*css
     &  - 2.d0*msfer2*s*zfbss*zfbss1r*css
     &  + 4.d0*msfer2*s*zfbsu*zfbsu1l*csu
     &  + 4.d0*msfer2*s*zfbsu*zfbsu1r*csu
     &  - 4.d0*msfer2*t*massorga*zfbtt2l*ctt
     &  - 4.d0*msfer2*t*massorga*zfbtt2r*ctt
     &  - 2.d0*msfer2*t*massorga*zfbuu6l*cuu
     &  - 2.d0*msfer2*t*massorga*zfbuu6r*cuu
     &  - 2.d0*msfer2*t*massorga*zfbuu8l*cuu
     &  - 2.d0*msfer2*t*massorga*zfbuu8r*cuu
     &  - 4.d0*msfer2*t*massorga*zfbtu3l*ctu
     &  - 4.d0*msfer2*t*massorga*zfbtu3r*ctu
     &  - 4.d0*msfer2*t*massorga*zfbtu4l*ctu
     &  - 4.d0*msfer2*t*massorga*zfbtu4r*ctu
     &  + 4.d0*msfer2*t*massorga*zfbts*zfbts4l*cts
     &  + 4.d0*msfer2*t*massorga*zfbts*zfbts4r*cts
     &  + 4.d0*msfer2*t*massorga*zfbsu*zfbsu6l*csu
     &  + 4.d0*msfer2*t*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*msfer2*t*zfbtt1l*ctt
     &  + 2.d0*msfer2*t*zfbtt1r*ctt
     &  + 2.d0*msfer2*t*zfbtu1l*ctu
     &  + 2.d0*msfer2*t*zfbtu1r*ctu
     &  + 4.d0*msfer2*t*zfbts*zfbts1l*cts
     &  + 4.d0*msfer2*t*zfbts*zfbts1r*cts
     &  + 4.d0*msfer2*t*zfbsu*zfbsu1l*csu
     &  + 4.d0*msfer2*t*zfbsu*zfbsu1r*csu
     &  + 2.d0*msfer2*t**2*massorga*zfbtt1l*ctt
     &  + 2.d0*msfer2*t**2*massorga*zfbtt1r*ctt
     &  + 2.d0*msfer2*t**2*massorga*zfbuu1l*cuu
     &  + 2.d0*msfer2*t**2*massorga*zfbuu1r*cuu
     &  + 4.d0*msfer2*t**2*massorga*zfbtu1l*ctu
     &  + 4.d0*msfer2*t**2*massorga*zfbtu1r*ctu
     &  - 4.d0*msfer2*zfbtt2l*ctt
     &  - 4.d0*msfer2*zfbtt2r*ctt
     &  + 2.d0*msfer2*zfbuu6l*cuu
     &  + 2.d0*msfer2*zfbuu6r*cuu
     &  - 8.d0*msfer2*zfbuu7l*cuu
     &  - 8.d0*msfer2*zfbuu7r*cuu
     &  + 2.d0*msfer2*zfbuu8l*cuu
     &  + 2.d0*msfer2*zfbuu8r*cuu
     &  + 4.d0*msfer2*zfbss*zfbss6l*css
     &  + 4.d0*msfer2*zfbss*zfbss6r*css
     &  + 4.d0*msfer2*zfbss*zfbss8l*css
     &  + 4.d0*msfer2*zfbss*zfbss8r*css
     &  - 2.d0*msfer2*zfbtu3l*ctu
     &  - 2.d0*msfer2*zfbtu3r*ctu
     &  - 6.d0*msfer2*zfbtu4l*ctu
     &  - 6.d0*msfer2*zfbtu4r*ctu
     &  - 4.d0*msfer2*zfbts*zfbts3l*cts
     &  - 4.d0*msfer2*zfbts*zfbts3r*cts
     &  + 2.d0*msfer2*zfbts*zfbts4l*cts
     &  + 2.d0*msfer2*zfbts*zfbts4r*cts
     &  - 4.d0*msfer2*zfbsu*zfbsu6l*csu
     &  - 4.d0*msfer2*zfbsu*zfbsu6r*csu
     &  + 4.d0*msfer2*zfbsu*zfbsu7l*csu
     &  + 4.d0*msfer2*zfbsu*zfbsu7r*csu
     &  - 8.d0*msfer2*zfbsu*zfbsu8l*csu
     &  - 8.d0*msfer2*zfbsu*zfbsu8r*csu
     &  + msfer2**2*mchi2*massorga*zfbtt1l*ctt
     &  + msfer2**2*mchi2*massorga*zfbtt1r*ctt
     &  + msfer2**2*mchi2*massorga*zfbuu1l*cuu
     &  + msfer2**2*mchi2*massorga*zfbuu1r*cuu
     &  + 2.d0*msfer2**2*mchi2*massorga*zfbtu1l*ctu
     &  + 2.d0*msfer2**2*mchi2*massorga*zfbtu1r*ctu
     &  + msfer2**2*mfer2*massorga*zfbtt1l*ctt
     &  + msfer2**2*mfer2*massorga*zfbtt1r*ctt
     &  + msfer2**2*mfer2*massorga*zfbuu1l*cuu
     &  + msfer2**2*mfer2*massorga*zfbuu1r*cuu
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbtu1l*ctu
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbtu1r*ctu
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbts*zfbts1l*cts
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbts*zfbts1r*cts
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*msfer2**2*mfer2*massorga*zfbsu*zfbsu1r*csu
     &  - msfer2**2*t*massorga*zfbtt1l*ctt
     &  - msfer2**2*t*massorga*zfbtt1r*ctt
     &  - msfer2**2*t*massorga*zfbuu1l*cuu
     &  - msfer2**2*t*massorga*zfbuu1r*cuu
     &  - 2.d0*msfer2**2*t*massorga*zfbtu1l*ctu
     &  - 2.d0*msfer2**2*t*massorga*zfbtu1r*ctu
     &  + 2.d0*msfer2**2*massorga*zfbtt2l*ctt
     &  + 2.d0*msfer2**2*massorga*zfbtt2r*ctt
     &  + msfer2**2*massorga*zfbuu6l*cuu
     &  + msfer2**2*massorga*zfbuu6r*cuu
     &  + msfer2**2*massorga*zfbuu8l*cuu
     &  + msfer2**2*massorga*zfbuu8r*cuu
     &  + 2.d0*msfer2**2*massorga*zfbtu3l*ctu
     &  + 2.d0*msfer2**2*massorga*zfbtu3r*ctu
     &  + 2.d0*msfer2**2*massorga*zfbtu4l*ctu
     &  + 2.d0*msfer2**2*massorga*zfbtu4r*ctu
     &  - 2.d0*msfer2**2*massorga*zfbts*zfbts4l*cts
     &  - 2.d0*msfer2**2*massorga*zfbts*zfbts4r*cts
     &  - 2.d0*msfer2**2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*msfer2**2*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*mchi2*mgb2*mfer2*massorga*zfbuu1l*cuu
     &  - 2.d0*mchi2*mgb2*mfer2*massorga*zfbuu1r*cuu
     &  - mchi2*mgb2*mfer2*massorga*zfbss*zfbss1l*css
     &  - mchi2*mgb2*mfer2*massorga*zfbss*zfbss1r*css
     &  - 4.d0*mchi2*mgb2*mfer2*massorga*zfbsu*zfbsu1l*csu
     &  - 4.d0*mchi2*mgb2*mfer2*massorga*zfbsu*zfbsu1r*csu
     &  + mchi2*mgb2*s*massorga*zfbuu1l*cuu
     &  + mchi2*mgb2*s*massorga*zfbuu1r*cuu
     &  + mchi2*mgb2*massorga*zfbuu2l*cuu
     &  + mchi2*mgb2*massorga*zfbuu2r*cuu
     &  - mchi2*mgb2*massorga*zfbuu3l*cuu
     &  - mchi2*mgb2*massorga*zfbuu3r*cuu
     &  + mchi2*mgb2*massorga*zfbuu4l*cuu
     &  + mchi2*mgb2*massorga*zfbuu4r*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfbuu6l*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfbuu6r*cuu
     &  + 2.d0*mchi2*mgb2*massorga*zfbuu7l*cuu
     &  + 2.d0*mchi2*mgb2*massorga*zfbuu7r*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfbuu8l*cuu
     &  - 2.d0*mchi2*mgb2*massorga*zfbuu8r*cuu
     &  - mchi2*mgb2*massorga*zfbss*zfbss3l*css
     &  - mchi2*mgb2*massorga*zfbss*zfbss3r*css
     &  + mchi2*mgb2*massorga*zfbss*zfbss6l*css
     &  + mchi2*mgb2*massorga*zfbss*zfbss6r*css
     &  + mchi2*mgb2*massorga*zfbss*zfbss8l*css
     &  + mchi2*mgb2*massorga*zfbss*zfbss8r*css
     &  + 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu3l*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu3r*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu4l*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu4r*csu
     &  + 4.d0*mchi2*mgb2*massorga*zfbsu*zfbsu6l*csu
     &  + 4.d0*mchi2*mgb2*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu7l*csu
     &  + 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu7r*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu8l*csu
     &  - 2.d0*mchi2*mgb2*massorga*zfbsu*zfbsu8r*csu
     &  + mchi2*mgb2*zfbtt1l*ctt
     &  + mchi2*mgb2*zfbtt1r*ctt
     &  + 2.d0*mchi2*mgb2*zfbuu1l*cuu
     &  + 2.d0*mchi2*mgb2*zfbuu1r*cuu
     &  - 2.d0*mchi2*mgb2*zfbss*zfbss1l*css
     &  - 2.d0*mchi2*mgb2*zfbss*zfbss1r*css
     &  + 4.d0*mchi2*mgb2*zfbts*zfbts1l*cts
     &  + 4.d0*mchi2*mgb2*zfbts*zfbts1r*cts
     &  + 4.d0*mchi2*mgb2*zfbsu*zfbsu1l*csu
     &  + 4.d0*mchi2*mgb2*zfbsu*zfbsu1r*csu
     &  - mchi2*mgb2**2*massorga*zfbuu1l*cuu
     &  - mchi2*mgb2**2*massorga*zfbuu1r*cuu
     &  - 2.d0*mchi2*mfer2*s*massorga*zfbss*zfbss1l*css
     &  - 2.d0*mchi2*mfer2*s*massorga*zfbss*zfbss1r*css
     &  - mchi2*mfer2*t*massorga*zfbuu1l*cuu
     &  - mchi2*mfer2*t*massorga*zfbuu1r*cuu
     &  + 2.d0*mchi2*mfer2*t*massorga*zfbts*zfbts1l*cts
     &  + 2.d0*mchi2*mfer2*t*massorga*zfbts*zfbts1r*cts
     &  + 2.d0*mchi2*mfer2*t*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*mchi2*mfer2*t*massorga*zfbsu*zfbsu1r*csu
     &  + 6.d0*mchi2*mfer2*zfbuu1l*cuu
     &  + 6.d0*mchi2*mfer2*zfbuu1r*cuu
     &  + 2.d0*mchi2*mfer2*zfbss*zfbss1l*css
     &  + 2.d0*mchi2*mfer2*zfbss*zfbss1r*css
     &  - 4.d0*mchi2*mfer2*zfbtu1l*ctu
     &  - 4.d0*mchi2*mfer2*zfbtu1r*ctu
     &  - 2.d0*mchi2*mfer2*zfbts*zfbts1l*cts
     &  - 2.d0*mchi2*mfer2*zfbts*zfbts1r*cts
     &  + 8.d0*mchi2*mfer2*zfbsu*zfbsu1l*csu
     &  + 8.d0*mchi2*mfer2*zfbsu*zfbsu1r*csu
     &  + mchi2*mfer2**2*massorga*zfbss*zfbss1l*css
     &  + mchi2*mfer2**2*massorga*zfbss*zfbss1r*css
     &  + mchi2*s*t*massorga*zfbuu1l*cuu
     &  + mchi2*s*t*massorga*zfbuu1r*cuu
     &  - 2.d0*mchi2*s*t*massorga*zfbts*zfbts1l*cts
     &  - 2.d0*mchi2*s*t*massorga*zfbts*zfbts1r*cts
     &  - 2.d0*mchi2*s*t*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*mchi2*s*t*massorga*zfbsu*zfbsu1r*csu
     &  - 4.d0*mchi2*s*zfbuu1l*cuu
     &  - 4.d0*mchi2*s*zfbuu1r*cuu
     &  - 2.d0*mchi2*s*zfbts*zfbts1l*cts
     &  - 2.d0*mchi2*s*zfbts*zfbts1r*cts
     &  + mchi2*s**2*massorga*zfbss*zfbss1l*css
     &  + mchi2*s**2*massorga*zfbss*zfbss1r*css
     &  - 2.d0*mchi2*t*zfbtt1l*ctt
     &  - 2.d0*mchi2*t*zfbtt1r*ctt
     &  - 2.d0*mchi2*t*zfbuu1l*cuu
     &  - 2.d0*mchi2*t*zfbuu1r*cuu
     &  - 2.d0*mchi2*t*zfbtu1l*ctu
     &  - 2.d0*mchi2*t*zfbtu1r*ctu
     &  - 4.d0*mchi2*t*zfbts*zfbts1l*cts
     &  - 4.d0*mchi2*t*zfbts*zfbts1r*cts
     &  - 4.d0*mchi2*t*zfbsu*zfbsu1l*csu
     &  - 4.d0*mchi2*t*zfbsu*zfbsu1r*csu
     &  + mchi2*t**2*massorga*zfbtt1l*ctt
     &  + mchi2*t**2*massorga*zfbtt1r*ctt
     &  + mchi2*t**2*massorga*zfbuu1l*cuu
     &  + mchi2*t**2*massorga*zfbuu1r*cuu
     &  + 2.d0*mchi2*t**2*massorga*zfbtu1l*ctu
     &  + 2.d0*mchi2*t**2*massorga*zfbtu1r*ctu
     &  - 4.d0*mchi2*zfbuu2l*cuu
     &  - 4.d0*mchi2*zfbuu2r*cuu
     &  + 2.d0*mchi2*zfbuu3l*cuu
     &  + 2.d0*mchi2*zfbuu3r*cuu
     &  - 4.d0*mchi2*zfbuu4l*cuu
     &  - 4.d0*mchi2*zfbuu4r*cuu
     &  + 4.d0*mchi2*zfbuu6l*cuu
     &  + 4.d0*mchi2*zfbuu6r*cuu
     &  - 8.d0*mchi2*zfbuu7l*cuu
     &  - 8.d0*mchi2*zfbuu7r*cuu
     &  + 4.d0*mchi2*zfbuu8l*cuu
     &  + 4.d0*mchi2*zfbuu8r*cuu
     &  + 2.d0*mchi2*zfbss*zfbss3l*css
     &  + 2.d0*mchi2*zfbss*zfbss3r*css
     &  - 4.d0*mchi2*zfbss*zfbss6l*css
     &  - 4.d0*mchi2*zfbss*zfbss6r*css
     &  - 4.d0*mchi2*zfbss*zfbss8l*css
     &  - 4.d0*mchi2*zfbss*zfbss8r*css
     &  - 4.d0*mchi2*zfbtu2l*ctu
     &  - 4.d0*mchi2*zfbtu2r*ctu
     &  - 4.d0*mchi2*zfbtu3l*ctu
     &  - 4.d0*mchi2*zfbtu3r*ctu
     &  - 4.d0*mchi2*zfbtu4l*ctu
     &  - 4.d0*mchi2*zfbtu4r*ctu
     &  + 4.d0*mchi2*zfbts*zfbts2l*cts
     &  + 4.d0*mchi2*zfbts*zfbts2r*cts
     &  + 4.d0*mchi2*zfbts*zfbts3l*cts
     &  + 4.d0*mchi2*zfbts*zfbts3r*cts
     &  + 4.d0*mchi2*zfbts*zfbts4l*cts
     &  + 4.d0*mchi2*zfbts*zfbts4r*cts
     &  - 4.d0*mchi2*zfbsu*zfbsu3l*csu
     &  - 4.d0*mchi2*zfbsu*zfbsu3r*csu
     &  + 8.d0*mchi2*zfbsu*zfbsu4l*csu
     &  + 8.d0*mchi2*zfbsu*zfbsu4r*csu
     &  - 8.d0*mchi2*zfbsu*zfbsu6l*csu
     &  - 8.d0*mchi2*zfbsu*zfbsu6r*csu
     &  - 4.d0*mchi2*zfbsu*zfbsu7l*csu
     &  - 4.d0*mchi2*zfbsu*zfbsu7r*csu
     &  + 8.d0*mchi2*zfbsu*zfbsu8l*csu
     &  + 8.d0*mchi2*zfbsu*zfbsu8r*csu
     &  - mchi2**2*mgb2*massorga*zfbuu1l*cuu
     &  - mchi2**2*mgb2*massorga*zfbuu1r*cuu
     &  + 2.d0*mchi2**2*zfbuu1l*cuu
     &  + 2.d0*mchi2**2*zfbuu1r*cuu
     &  + 4.d0*mchi2**2*zfbts*zfbts1l*cts
     &  + 4.d0*mchi2**2*zfbts*zfbts1r*cts
     &  + 4.d0*mchi2**2*zfbsu*zfbsu1l*csu
     &  + 4.d0*mchi2**2*zfbsu*zfbsu1r*csu
     &  - 2.d0*mgb2*mfer2*s*massorga*zfbuu1l*cuu
     &  - 2.d0*mgb2*mfer2*s*massorga*zfbuu1r*cuu
     &  - mgb2*mfer2*s*massorga*zfbss*zfbss1l*css
     &  - mgb2*mfer2*s*massorga*zfbss*zfbss1r*css
     &  - 4.d0*mgb2*mfer2*t*massorga*zfbuu1l*cuu
     &  - 4.d0*mgb2*mfer2*t*massorga*zfbuu1r*cuu
     &  - 2.d0*mgb2*mfer2*t*massorga*zfbtu1l*ctu
     &  - 2.d0*mgb2*mfer2*t*massorga*zfbtu1r*ctu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu2l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu2r*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu4l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu4r*cuu
     &  + mgb2*mfer2*massorga*zfbuu6l*cuu
     &  + mgb2*mfer2*massorga*zfbuu6r*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu7l*cuu
     &  + 2.d0*mgb2*mfer2*massorga*zfbuu7r*cuu
     &  + mgb2*mfer2*massorga*zfbuu8l*cuu
     &  + mgb2*mfer2*massorga*zfbuu8r*cuu
     &  - mgb2*mfer2*massorga*zfbss*zfbss2l*css
     &  - mgb2*mfer2*massorga*zfbss*zfbss2r*css
     &  - mgb2*mfer2*massorga*zfbss*zfbss4l*css
     &  - mgb2*mfer2*massorga*zfbss*zfbss4r*css
     &  + 2.d0*mgb2*mfer2*massorga*zfbsu*zfbsu2l*csu
     &  + 2.d0*mgb2*mfer2*massorga*zfbsu*zfbsu2r*csu
     &  - 4.d0*mgb2*mfer2*massorga*zfbsu*zfbsu4l*csu
     &  - 4.d0*mgb2*mfer2*massorga*zfbsu*zfbsu4r*csu
     &  - 2.d0*mgb2*mfer2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*mgb2*mfer2*massorga*zfbsu*zfbsu6r*csu
     &  + mgb2*mfer2*zfbtt1l*ctt
     &  + mgb2*mfer2*zfbtt1r*ctt
     &  + 2.d0*mgb2*mfer2*zfbtu1l*ctu
     &  + 2.d0*mgb2*mfer2*zfbtu1r*ctu
     &  + 2.d0*mgb2*mfer2**2*massorga*zfbuu1l*cuu
     &  + 2.d0*mgb2*mfer2**2*massorga*zfbuu1r*cuu
     &  + 2.d0*mgb2*s*t*massorga*zfbuu1l*cuu
     &  + 2.d0*mgb2*s*t*massorga*zfbuu1r*cuu
     &  - 2.d0*mgb2*s*t*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*mgb2*s*t*massorga*zfbsu*zfbsu1r*csu
     &  - mgb2*s*massorga*zfbuu2l*cuu
     &  - mgb2*s*massorga*zfbuu2r*cuu
     &  - mgb2*s*massorga*zfbuu3l*cuu
     &  - mgb2*s*massorga*zfbuu3r*cuu
     &  - mgb2*s*massorga*zfbuu4l*cuu
     &  - mgb2*s*massorga*zfbuu4r*cuu
     &  - mgb2*s*massorga*zfbuu6l*cuu
     &  - mgb2*s*massorga*zfbuu6r*cuu
     &  - 2.d0*mgb2*s*massorga*zfbuu7l*cuu
     &  - 2.d0*mgb2*s*massorga*zfbuu7r*cuu
     &  - mgb2*s*massorga*zfbuu8l*cuu
     &  - mgb2*s*massorga*zfbuu8r*cuu
     &  - mgb2*s*massorga*zfbss*zfbss2l*css
     &  - mgb2*s*massorga*zfbss*zfbss2r*css
     &  - mgb2*s*massorga*zfbss*zfbss3l*css
     &  - mgb2*s*massorga*zfbss*zfbss3r*css
     &  - mgb2*s*massorga*zfbss*zfbss4l*css
     &  - mgb2*s*massorga*zfbss*zfbss4r*css
     &  + mgb2*s*massorga*zfbss*zfbss6l*css
     &  + mgb2*s*massorga*zfbss*zfbss6r*css
     &  + 2.d0*mgb2*s*massorga*zfbss*zfbss7l*css
     &  + 2.d0*mgb2*s*massorga*zfbss*zfbss7r*css
     &  + mgb2*s*massorga*zfbss*zfbss8l*css
     &  + mgb2*s*massorga*zfbss*zfbss8r*css
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu2l*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu2r*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu3l*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu3r*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu4l*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu4r*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu6l*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu7l*csu
     &  + 2.d0*mgb2*s*massorga*zfbsu*zfbsu7r*csu
     &  - 2.d0*mgb2*s*massorga*zfbsu*zfbsu8l*csu
     &  - 2.d0*mgb2*s*massorga*zfbsu*zfbsu8r*csu
     &  - 2.d0*mgb2*s*zfbuu1l*cuu
     &  - 2.d0*mgb2*s*zfbuu1r*cuu
     &  - 2.d0*mgb2*s*zfbss*zfbss1l*css
     &  - 2.d0*mgb2*s*zfbss*zfbss1r*css
     &  + 4.d0*mgb2*s*zfbsu*zfbsu1l*csu
     &  + 4.d0*mgb2*s*zfbsu*zfbsu1r*csu
      tmp2=
     &  - mgb2*t*massorga*zfbuu2l*cuu
     &  - mgb2*t*massorga*zfbuu2r*cuu
     &  - mgb2*t*massorga*zfbuu4l*cuu
     &  - mgb2*t*massorga*zfbuu4r*cuu
     &  - mgb2*t*massorga*zfbuu6l*cuu
     &  - mgb2*t*massorga*zfbuu6r*cuu
     &  - 2.d0*mgb2*t*massorga*zfbuu7l*cuu
     &  - 2.d0*mgb2*t*massorga*zfbuu7r*cuu
     &  - mgb2*t*massorga*zfbuu8l*cuu
     &  - mgb2*t*massorga*zfbuu8r*cuu
     &  - 2.d0*mgb2*t*massorga*zfbtu2l*ctu
     &  - 2.d0*mgb2*t*massorga*zfbtu2r*ctu
     &  - 4.d0*mgb2*t*massorga*zfbtu4l*ctu
     &  - 4.d0*mgb2*t*massorga*zfbtu4r*ctu
     &  + 2.d0*mgb2*t*massorga*zfbts*zfbts2l*cts
     &  + 2.d0*mgb2*t*massorga*zfbts*zfbts2r*cts
     &  - 2.d0*mgb2*t*massorga*zfbts*zfbts3l*cts
     &  - 2.d0*mgb2*t*massorga*zfbts*zfbts3r*cts
     &  + 2.d0*mgb2*t*massorga*zfbsu*zfbsu4l*csu
     &  + 2.d0*mgb2*t*massorga*zfbsu*zfbsu4r*csu
     &  + 2.d0*mgb2*t*massorga*zfbsu*zfbsu6l*csu
     &  + 2.d0*mgb2*t*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*mgb2*t*massorga*zfbsu*zfbsu7l*csu
     &  - 2.d0*mgb2*t*massorga*zfbsu*zfbsu7r*csu
     &  - mgb2*t*zfbtt1l*ctt
     &  - mgb2*t*zfbtt1r*ctt
     &  - 2.d0*mgb2*t*zfbtu1l*ctu
     &  - 2.d0*mgb2*t*zfbtu1r*ctu
     &  + 2.d0*mgb2*t**2*massorga*zfbuu1l*cuu
     &  + 2.d0*mgb2*t**2*massorga*zfbuu1r*cuu
     &  + 2.d0*mgb2*t**2*massorga*zfbtu1l*ctu
     &  + 2.d0*mgb2*t**2*massorga*zfbtu1r*ctu
     &  + 2.d0*mgb2*massorga*zfbuu5l*cuu
     &  + 2.d0*mgb2*massorga*zfbuu5r*cuu
     &  + 2.d0*mgb2*massorga*zfbss*zfbss5l*css
     &  + 2.d0*mgb2*massorga*zfbss*zfbss5r*css
     &  - 4.d0*mgb2*massorga*zfbsu*zfbsu5l*csu
     &  - 4.d0*mgb2*massorga*zfbsu*zfbsu5r*csu
     &  + 2.d0*mgb2*zfbtt2l*ctt
     &  + 2.d0*mgb2*zfbtt2r*ctt
     &  - 4.d0*mgb2*zfbuu2l*cuu
     &  - 4.d0*mgb2*zfbuu2r*cuu
     &  - 4.d0*mgb2*zfbuu4l*cuu
     &  - 4.d0*mgb2*zfbuu4r*cuu
     &  - 8.d0*mgb2*zfbuu7l*cuu
     &  - 8.d0*mgb2*zfbuu7r*cuu
     &  - 2.d0*mgb2*zfbss*zfbss2l*css
     &  - 2.d0*mgb2*zfbss*zfbss2r*css
     &  - 2.d0*mgb2*zfbss*zfbss4l*css
     &  - 2.d0*mgb2*zfbss*zfbss4r*css
     &  - 2.d0*mgb2*zfbtu2l*ctu
     &  - 2.d0*mgb2*zfbtu2r*ctu
     &  + 2.d0*mgb2*zfbts*zfbts2l*cts
     &  + 2.d0*mgb2*zfbts*zfbts2r*cts
     &  + 2.d0*mgb2*zfbts*zfbts3l*cts
     &  + 2.d0*mgb2*zfbts*zfbts3r*cts
     &  + 4.d0*mgb2*zfbsu*zfbsu2l*csu
     &  + 4.d0*mgb2*zfbsu*zfbsu2r*csu
     &  + 8.d0*mgb2*zfbsu*zfbsu4l*csu
     &  + 8.d0*mgb2*zfbsu*zfbsu4r*csu
     &  + 4.d0*mgb2*zfbsu*zfbsu7l*csu
     &  + 4.d0*mgb2*zfbsu*zfbsu7r*csu
     &  + mgb2**2*mfer2*massorga*zfbuu1l*cuu
     &  + mgb2**2*mfer2*massorga*zfbuu1r*cuu
     &  - mgb2**2*t*massorga*zfbuu1l*cuu
     &  - mgb2**2*t*massorga*zfbuu1r*cuu
     &  + mgb2**2*massorga*zfbuu2l*cuu
     &  + mgb2**2*massorga*zfbuu2r*cuu
     &  + mgb2**2*massorga*zfbuu4l*cuu
     &  + mgb2**2*massorga*zfbuu4r*cuu
     &  + 2.d0*mgb2**2*massorga*zfbuu7l*cuu
     &  + 2.d0*mgb2**2*massorga*zfbuu7r*cuu
     &  - 2.d0*mgb2**2*massorga*zfbsu*zfbsu4l*csu
     &  - 2.d0*mgb2**2*massorga*zfbsu*zfbsu4r*csu
     &  + 2.d0*mgb2**2*massorga*zfbsu*zfbsu7l*csu
     &  + 2.d0*mgb2**2*massorga*zfbsu*zfbsu7r*csu
     &  + 4.d0*mfer2*s*t*massorga*zfbuu1l*cuu
     &  + 4.d0*mfer2*s*t*massorga*zfbuu1r*cuu
     &  + mfer2*s*t*massorga*zfbss*zfbss1l*css
     &  + mfer2*s*t*massorga*zfbss*zfbss1r*css
     &  + 2.d0*mfer2*s*t*massorga*zfbtu1l*ctu
     &  + 2.d0*mfer2*s*t*massorga*zfbtu1r*ctu
     &  - 2.d0*mfer2*s*t*massorga*zfbsu*zfbsu1l*csu
     &  - 2.d0*mfer2*s*t*massorga*zfbsu*zfbsu1r*csu
     &  - 2.d0*mfer2*s*massorga*zfbuu3l*cuu
     &  - 2.d0*mfer2*s*massorga*zfbuu3r*cuu
     &  - 2.d0*mfer2*s*massorga*zfbuu6l*cuu
     &  - 2.d0*mfer2*s*massorga*zfbuu6r*cuu
     &  - 2.d0*mfer2*s*massorga*zfbuu8l*cuu
     &  - 2.d0*mfer2*s*massorga*zfbuu8r*cuu
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss2l*css
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss2r*css
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss3l*css
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss3r*css
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss4l*css
     &  - 2.d0*mfer2*s*massorga*zfbss*zfbss4r*css
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu2l*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu2r*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu3l*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu3r*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu6l*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu6r*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu7l*csu
     &  + 4.d0*mfer2*s*massorga*zfbsu*zfbsu7r*csu
     &  - 4.d0*mfer2*s*zfbuu1l*cuu
     &  - 4.d0*mfer2*s*zfbuu1r*cuu
     &  + 2.d0*mfer2*s*zfbtu1l*ctu
     &  + 2.d0*mfer2*s*zfbtu1r*ctu
     &  + 8.d0*mfer2*s*zfbsu*zfbsu1l*csu
     &  + 8.d0*mfer2*s*zfbsu*zfbsu1r*csu
     &  + mfer2*s**2*massorga*zfbuu1l*cuu
     &  + mfer2*s**2*massorga*zfbuu1r*cuu
     &  - mfer2*t*massorga*zfbuu3l*cuu
     &  - mfer2*t*massorga*zfbuu3r*cuu
     &  - 2.d0*mfer2*t*massorga*zfbuu6l*cuu
     &  - 2.d0*mfer2*t*massorga*zfbuu6r*cuu
     &  - 2.d0*mfer2*t*massorga*zfbuu8l*cuu
     &  - 2.d0*mfer2*t*massorga*zfbuu8r*cuu
     &  - mfer2*t*massorga*zfbss*zfbss3l*css
     &  - mfer2*t*massorga*zfbss*zfbss3r*css
     &  - 2.d0*mfer2*t*massorga*zfbtu2l*ctu
     &  - 2.d0*mfer2*t*massorga*zfbtu2r*ctu
     &  - 2.d0*mfer2*t*massorga*zfbtu3l*ctu
     &  - 2.d0*mfer2*t*massorga*zfbtu3r*ctu
     &  - 2.d0*mfer2*t*massorga*zfbtu4l*ctu
     &  - 2.d0*mfer2*t*massorga*zfbtu4r*ctu
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts2l*cts
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts2r*cts
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts3l*cts
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts3r*cts
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts4l*cts
     &  + 2.d0*mfer2*t*massorga*zfbts*zfbts4r*cts
     &  + 2.d0*mfer2*t*massorga*zfbsu*zfbsu3l*csu
     &  + 2.d0*mfer2*t*massorga*zfbsu*zfbsu3r*csu
     &  + 4.d0*mfer2*t*massorga*zfbsu*zfbsu6l*csu
     &  + 4.d0*mfer2*t*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*mfer2*t*massorga*zfbsu*zfbsu7l*csu
     &  + 2.d0*mfer2*t*massorga*zfbsu*zfbsu7r*csu
     &  - 2.d0*mfer2*t*zfbtt1l*ctt
     &  - 2.d0*mfer2*t*zfbtt1r*ctt
     &  - 2.d0*mfer2*t*zfbuu1l*cuu
     &  - 2.d0*mfer2*t*zfbuu1r*cuu
     &  - 4.d0*mfer2*t*zfbtu1l*ctu
     &  - 4.d0*mfer2*t*zfbtu1r*ctu
     &  + mfer2*t**2*massorga*zfbtt1l*ctt
     &  + mfer2*t**2*massorga*zfbtt1r*ctt
     &  + 3.d0*mfer2*t**2*massorga*zfbuu1l*cuu
     &  + 3.d0*mfer2*t**2*massorga*zfbuu1r*cuu
     &  + 4.d0*mfer2*t**2*massorga*zfbtu1l*ctu
     &  + 4.d0*mfer2*t**2*massorga*zfbtu1r*ctu
     &  - 8.d0*mfer2*zfbuu2l*cuu
     &  - 8.d0*mfer2*zfbuu2r*cuu
     &  + 2.d0*mfer2*zfbuu3l*cuu
     &  + 2.d0*mfer2*zfbuu3r*cuu
     &  - 8.d0*mfer2*zfbuu4l*cuu
     &  - 8.d0*mfer2*zfbuu4r*cuu
     &  + 2.d0*mfer2*zfbuu6l*cuu
     &  + 2.d0*mfer2*zfbuu6r*cuu
     &  - 8.d0*mfer2*zfbuu7l*cuu
     &  - 8.d0*mfer2*zfbuu7r*cuu
     &  + 2.d0*mfer2*zfbuu8l*cuu
     &  + 2.d0*mfer2*zfbuu8r*cuu
     &  + 2.d0*mfer2*zfbss*zfbss2l*css
     &  + 2.d0*mfer2*zfbss*zfbss2r*css
     &  + 2.d0*mfer2*zfbss*zfbss3l*css
     &  + 2.d0*mfer2*zfbss*zfbss3r*css
     &  + 2.d0*mfer2*zfbss*zfbss4l*css
     &  + 2.d0*mfer2*zfbss*zfbss4r*css
     &  + 2.d0*mfer2*zfbtu2l*ctu
     &  + 2.d0*mfer2*zfbtu2r*ctu
     &  + 2.d0*mfer2*zfbtu3l*ctu
     &  + 2.d0*mfer2*zfbtu3r*ctu
     &  + 2.d0*mfer2*zfbtu4l*ctu
     &  + 2.d0*mfer2*zfbtu4r*ctu
     &  - 2.d0*mfer2*zfbts*zfbts2l*cts
     &  - 2.d0*mfer2*zfbts*zfbts2r*cts
     &  - 2.d0*mfer2*zfbts*zfbts3l*cts
     &  - 2.d0*mfer2*zfbts*zfbts3r*cts
     &  - 2.d0*mfer2*zfbts*zfbts4l*cts
     &  - 2.d0*mfer2*zfbts*zfbts4r*cts
     &  - 4.d0*mfer2*zfbsu*zfbsu2l*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu2r*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu3l*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu3r*csu
     &  + 16.d0*mfer2*zfbsu*zfbsu4l*csu
     &  + 16.d0*mfer2*zfbsu*zfbsu4r*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu6l*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu6r*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu7l*csu
     &  - 4.d0*mfer2*zfbsu*zfbsu7r*csu
     &  - 2.d0*mfer2**2*s*massorga*zfbuu1l*cuu
     &  - 2.d0*mfer2**2*s*massorga*zfbuu1r*cuu
     &  - 3.d0*mfer2**2*t*massorga*zfbuu1l*cuu
     &  - 3.d0*mfer2**2*t*massorga*zfbuu1r*cuu
     &  - 2.d0*mfer2**2*t*massorga*zfbtu1l*ctu
     &  - 2.d0*mfer2**2*t*massorga*zfbtu1r*ctu
     &  + mfer2**2*massorga*zfbuu3l*cuu
     &  + mfer2**2*massorga*zfbuu3r*cuu
     &  + mfer2**2*massorga*zfbuu6l*cuu
     &  + mfer2**2*massorga*zfbuu6r*cuu
     &  + mfer2**2*massorga*zfbuu8l*cuu
     &  + mfer2**2*massorga*zfbuu8r*cuu
     &  + mfer2**2*massorga*zfbss*zfbss2l*css
     &  + mfer2**2*massorga*zfbss*zfbss2r*css
     &  + mfer2**2*massorga*zfbss*zfbss3l*css
     &  + mfer2**2*massorga*zfbss*zfbss3r*css
     &  + mfer2**2*massorga*zfbss*zfbss4l*css
     &  + mfer2**2*massorga*zfbss*zfbss4r*css
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu2l*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu2r*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu3l*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu3r*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu7l*csu
     &  - 2.d0*mfer2**2*massorga*zfbsu*zfbsu7r*csu
     &  + 2.d0*mfer2**2*zfbuu1l*cuu
     &  + 2.d0*mfer2**2*zfbuu1r*cuu
     &  + 2.d0*mfer2**2*zfbtu1l*ctu
     &  + 2.d0*mfer2**2*zfbtu1r*ctu
     &  + mfer2**3*massorga*zfbuu1l*cuu
     &  + mfer2**3*massorga*zfbuu1r*cuu
     &  + s*t*massorga*zfbuu3l*cuu
     &  + s*t*massorga*zfbuu3r*cuu
     &  + 2.d0*s*t*massorga*zfbuu6l*cuu
     &  + 2.d0*s*t*massorga*zfbuu6r*cuu
     &  + 2.d0*s*t*massorga*zfbuu8l*cuu
     &  + 2.d0*s*t*massorga*zfbuu8r*cuu
     &  + s*t*massorga*zfbss*zfbss3l*css
     &  + s*t*massorga*zfbss*zfbss3r*css
     &  + 2.d0*s*t*massorga*zfbtu2l*ctu
     &  + 2.d0*s*t*massorga*zfbtu2r*ctu
     &  + 2.d0*s*t*massorga*zfbtu3l*ctu
     &  + 2.d0*s*t*massorga*zfbtu3r*ctu
     &  + 2.d0*s*t*massorga*zfbtu4l*ctu
     &  + 2.d0*s*t*massorga*zfbtu4r*ctu
     &  - 2.d0*s*t*massorga*zfbts*zfbts2l*cts
     &  - 2.d0*s*t*massorga*zfbts*zfbts2r*cts
     &  - 2.d0*s*t*massorga*zfbts*zfbts3l*cts
     &  - 2.d0*s*t*massorga*zfbts*zfbts3r*cts
     &  - 2.d0*s*t*massorga*zfbts*zfbts4l*cts
     &  - 2.d0*s*t*massorga*zfbts*zfbts4r*cts
     &  - 2.d0*s*t*massorga*zfbsu*zfbsu3l*csu
     &  - 2.d0*s*t*massorga*zfbsu*zfbsu3r*csu
     &  - 4.d0*s*t*massorga*zfbsu*zfbsu6l*csu
     &  - 4.d0*s*t*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*s*t*massorga*zfbsu*zfbsu7l*csu
     &  - 2.d0*s*t*massorga*zfbsu*zfbsu7r*csu
     &  + 2.d0*s*t*zfbuu1l*cuu
     &  + 2.d0*s*t*zfbuu1r*cuu
     &  + 2.d0*s*t*zfbss*zfbss1l*css
     &  + 2.d0*s*t*zfbss*zfbss1r*css
     &  + 2.d0*s*t*zfbtu1l*ctu
     &  + 2.d0*s*t*zfbtu1r*ctu
     &  - 2.d0*s*t*zfbts*zfbts1l*cts
     &  - 2.d0*s*t*zfbts*zfbts1r*cts
     &  - 4.d0*s*t*zfbsu*zfbsu1l*csu
     &  - 4.d0*s*t*zfbsu*zfbsu1r*csu
     &  - 2.d0*s*t**2*massorga*zfbuu1l*cuu
     &  - 2.d0*s*t**2*massorga*zfbuu1r*cuu
     &  - 2.d0*s*t**2*massorga*zfbtu1l*ctu
     &  - 2.d0*s*t**2*massorga*zfbtu1r*ctu
     &  + 2.d0*s*t**2*massorga*zfbts*zfbts1l*cts
     &  + 2.d0*s*t**2*massorga*zfbts*zfbts1r*cts
     &  + 2.d0*s*t**2*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*s*t**2*massorga*zfbsu*zfbsu1r*csu
     &  + 4.d0*s*zfbuu2l*cuu
     &  + 4.d0*s*zfbuu2r*cuu
     &  + 4.d0*s*zfbuu4l*cuu
     &  + 4.d0*s*zfbuu4r*cuu
     &  - 2.d0*s*zfbuu6l*cuu
     &  - 2.d0*s*zfbuu6r*cuu
     &  + 8.d0*s*zfbuu7l*cuu
     &  + 8.d0*s*zfbuu7r*cuu
     &  - 2.d0*s*zfbuu8l*cuu
     &  - 2.d0*s*zfbuu8r*cuu
     &  + 2.d0*s*zfbss*zfbss2l*css
     &  + 2.d0*s*zfbss*zfbss2r*css
     &  + 2.d0*s*zfbss*zfbss4l*css
     &  + 2.d0*s*zfbss*zfbss4r*css
     &  - 4.d0*s*zfbss*zfbss6l*css
     &  - 4.d0*s*zfbss*zfbss6r*css
     &  - 8.d0*s*zfbss*zfbss7l*css
     &  - 8.d0*s*zfbss*zfbss7r*css
     &  - 4.d0*s*zfbss*zfbss8l*css
     &  - 4.d0*s*zfbss*zfbss8r*css
     &  + 2.d0*s*zfbtu2l*ctu
     &  + 2.d0*s*zfbtu2r*ctu
     &  + 2.d0*s*zfbtu3l*ctu
     &  + 2.d0*s*zfbtu3r*ctu
     &  + 2.d0*s*zfbtu4l*ctu
     &  + 2.d0*s*zfbtu4r*ctu
     &  - 2.d0*s*zfbts*zfbts2l*cts
     &  - 2.d0*s*zfbts*zfbts2r*cts
     &  - 2.d0*s*zfbts*zfbts3l*cts
     &  - 2.d0*s*zfbts*zfbts3r*cts
     &  - 2.d0*s*zfbts*zfbts4l*cts
     &  - 2.d0*s*zfbts*zfbts4r*cts
     &  - 4.d0*s*zfbsu*zfbsu2l*csu
     &  - 4.d0*s*zfbsu*zfbsu2r*csu
     &  - 8.d0*s*zfbsu*zfbsu4l*csu
     &  - 8.d0*s*zfbsu*zfbsu4r*csu
     &  + 4.d0*s*zfbsu*zfbsu6l*csu
     &  + 4.d0*s*zfbsu*zfbsu6r*csu
     &  + 8.d0*s*zfbsu*zfbsu8l*csu
     &  + 8.d0*s*zfbsu*zfbsu8r*csu
     &  - s**2*t*massorga*zfbuu1l*cuu
     &  - s**2*t*massorga*zfbuu1r*cuu
     &  - s**2*t*massorga*zfbss*zfbss1l*css
     &  - s**2*t*massorga*zfbss*zfbss1r*css
     &  + 2.d0*s**2*t*massorga*zfbsu*zfbsu1l*csu
     &  + 2.d0*s**2*t*massorga*zfbsu*zfbsu1r*csu
     &  + s**2*massorga*zfbuu3l*cuu
     &  + s**2*massorga*zfbuu3r*cuu
     &  + s**2*massorga*zfbuu6l*cuu
     &  + s**2*massorga*zfbuu6r*cuu
     &  + s**2*massorga*zfbuu8l*cuu
     &  + s**2*massorga*zfbuu8r*cuu
     &  + s**2*massorga*zfbss*zfbss2l*css
     &  + s**2*massorga*zfbss*zfbss2r*css
     &  + s**2*massorga*zfbss*zfbss3l*css
     &  + s**2*massorga*zfbss*zfbss3r*css
     &  + s**2*massorga*zfbss*zfbss4l*css
     &  + s**2*massorga*zfbss*zfbss4r*css
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu2l*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu2r*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu3l*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu3r*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu6r*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu7l*csu
     &  - 2.d0*s**2*massorga*zfbsu*zfbsu7r*csu
     &  + 2.d0*s**2*zfbuu1l*cuu
     &  + 2.d0*s**2*zfbuu1r*cuu
     &  + 2.d0*s**2*zfbss*zfbss1l*css
     &  + 2.d0*s**2*zfbss*zfbss1r*css
     &  - 4.d0*s**2*zfbsu*zfbsu1l*csu
     &  - 4.d0*s**2*zfbsu*zfbsu1r*csu
     &  - 4.d0*t*zfbtt2l*ctt
     &  - 4.d0*t*zfbtt2r*ctt
     &  + 4.d0*t*zfbuu2l*cuu
     &  + 4.d0*t*zfbuu2r*cuu
     &  - 2.d0*t*zfbuu3l*cuu
     &  - 2.d0*t*zfbuu3r*cuu
     &  + 4.d0*t*zfbuu4l*cuu
     &  + 4.d0*t*zfbuu4r*cuu
     &  - 2.d0*t*zfbuu6l*cuu
     &  - 2.d0*t*zfbuu6r*cuu
     &  + 8.d0*t*zfbuu7l*cuu
     &  + 8.d0*t*zfbuu7r*cuu
     &  - 2.d0*t*zfbuu8l*cuu
     &  - 2.d0*t*zfbuu8r*cuu
     &  - 2.d0*t*zfbss*zfbss3l*css
     &  - 2.d0*t*zfbss*zfbss3r*css
     &  + 4.d0*t*zfbtu2l*ctu
     &  + 4.d0*t*zfbtu2r*ctu
     &  - 2.d0*t*zfbtu3l*ctu
     &  - 2.d0*t*zfbtu3r*ctu
     &  + 2.d0*t*zfbtu4l*ctu
     &  + 2.d0*t*zfbtu4r*ctu
     &  - 4.d0*t*zfbts*zfbts2l*cts
     &  - 4.d0*t*zfbts*zfbts2r*cts
     &  + 2.d0*t*zfbts*zfbts4l*cts
     &  + 2.d0*t*zfbts*zfbts4r*cts
     &  + 4.d0*t*zfbsu*zfbsu3l*csu
     &  + 4.d0*t*zfbsu*zfbsu3r*csu
     &  - 8.d0*t*zfbsu*zfbsu4l*csu
     &  - 8.d0*t*zfbsu*zfbsu4r*csu
     &  + 4.d0*t*zfbsu*zfbsu6l*csu
     &  + 4.d0*t*zfbsu*zfbsu6r*csu
     &  + 2.d0*t**2*massorga*zfbtt2l*ctt
     &  + 2.d0*t**2*massorga*zfbtt2r*ctt
     &  + t**2*massorga*zfbuu6l*cuu
     &  + t**2*massorga*zfbuu6r*cuu
     &  + t**2*massorga*zfbuu8l*cuu
     &  + t**2*massorga*zfbuu8r*cuu
     &  + 2.d0*t**2*massorga*zfbtu3l*ctu
     &  + 2.d0*t**2*massorga*zfbtu3r*ctu
     &  + 2.d0*t**2*massorga*zfbtu4l*ctu
     &  + 2.d0*t**2*massorga*zfbtu4r*ctu
     &  - 2.d0*t**2*massorga*zfbts*zfbts4l*cts
     &  - 2.d0*t**2*massorga*zfbts*zfbts4r*cts
     &  - 2.d0*t**2*massorga*zfbsu*zfbsu6l*csu
     &  - 2.d0*t**2*massorga*zfbsu*zfbsu6r*csu
     &  + 2.d0*t**2*zfbtt1l*ctt
     &  + 2.d0*t**2*zfbtt1r*ctt
     &  + 2.d0*t**2*zfbtu1l*ctu
     &  + 2.d0*t**2*zfbtu1r*ctu
     &  - t**3*massorga*zfbtt1l*ctt
     &  - t**3*massorga*zfbtt1r*ctt
     &  - t**3*massorga*zfbuu1l*cuu
     &  - t**3*massorga*zfbuu1r*cuu
     &  - 2.d0*t**3*massorga*zfbtu1l*ctu
     &  - 2.d0*t**3*massorga*zfbtu1r*ctu
     &  - 8.d0*zfbuu5l*cuu
     &  - 8.d0*zfbuu5r*cuu
     &  - 8.d0*zfbss*zfbss5l*css
     &  - 8.d0*zfbss*zfbss5r*css
     &  + 16.d0*zfbsu*zfbsu5l*csu
     &  + 16.d0*zfbsu*zfbsu5r*csu
      dsaschicasebas = tmp1+tmp2

      if (dble(dsaschicasebas).lt.0.0d0) then     
        if(aszeroprint) then
        write(*,*) ' '
        write(*,*) 
     &    'DS: ERROR IN dsaschicaseb with negative cross section:'
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
        write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
        write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
        write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
        write(*,*) 'DS: dsaschicasebas = ',dsaschicasebas
        write(*,*) 'DS: mass1=',mass1
        write(*,*) 'DS: mass2=',mass2
        write(*,*) 'DS: mass3=',mass3
        write(*,*) 'DS: mass4=',mass4
        write(*,*) 'DS: s=',s
        write(*,*) 'DS: t=',t
        write(*,*) 'DS: u=',u
        write(*,*) 'DS: p12=',p12
        write(*,*) 'DS: costheta=',costheta
        endif
        dsaschicasebas=dcmplx(0.0d0,0.0d0)
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsaschicasebas)*k34/(8.0d0*pi*gg1c*gg2c*s34*dsqrt(s))
      return
      end

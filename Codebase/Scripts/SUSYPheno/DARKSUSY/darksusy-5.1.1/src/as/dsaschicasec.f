c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.34, October 8, 2001, edsjo@physto.se)
c....Template file for dsaschicasec begins here

**************************************************************
*** SUBROUTINE dsaschicasec                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** sfermion(i) + neutralino(j)/chargino^+(j)              ***
*** -> higgs-boson + fermion                               ***
***                                                        ***
*** The sfermion must be the first mentioned               ***
*** particle (kp1) and the neutralino/chargino             ***
*** the other (kp2) -- not the opposite.                   ***
*** For the final state the higgs boson must be mentioned  ***
*** first (i.e. kp3) and next the fermion (kp4) --         ***
*** not the opposite.                                      ***
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-10                                         ***
*** Trivial color factors included: 02-03-21               ***
*** Added flavour changing charged exchange:               ***
*** for the 3rd, 4th, 5th and 6th case                     ***
*** added by Mia Schelke January 2007                      ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

***** Note that it is assumed that coupling constants that do 
***** not exist have already been set to zero!!!!!
***** Thus, many of the coefficients defined in this code 
***** simplify when the diagrams contain sneutrinos 
***** or neutrinos.

      subroutine dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      complex*16  dsaschicasecas
      real*8      msfer2,mchi2,mhb2,mfer2
      complex*16  h1ftt,h1ftt1l,h1ftt1r,h1ftt2l,h1ftt2r
      complex*16  h1fuu1l,h1fuu1r,h1fuu2l,h1fuu2r
      complex*16  h1fuu3l,h1fuu3r,h1fuu4l,h1fuu4r
      complex*16  h1fuu5l,h1fuu5r,h1fuu6l,h1fuu6r
      complex*16  h1fuu7l,h1fuu7r,h1fuu8l,h1fuu8r
      complex*16  h1fss,h1fss1l,h1fss1r,h1fss2l
      complex*16  h1fss2r,h1fss3l,h1fss3r,h1fss4l
      complex*16  h1fss4r,h1fss5l,h1fss5r,h1fss6l
      complex*16  h1fss6r,h1fss7l,h1fss7r,h1fss8l,h1fss8r
      complex*16  h1ftu,h1ftu1l,h1ftu1r,h1ftu2l,h1ftu2r
      complex*16  h1ftu3l,h1ftu3r,h1ftu4l,h1ftu4r
      complex*16  h1fts,h1fts1l,h1fts1r,h1fts2l,h1fts2r
      complex*16  h1fts3l,h1fts3r,h1fts4l,h1fts4r
      complex*16  h1fsu,h1fsu1l,h1fsu1r,h1fsu2l,h1fsu2r 
      complex*16  h1fsu3l,h1fsu3r,h1fsu4l,h1fsu4r 
      complex*16  h1fsu5l,h1fsu5r,h1fsu6l,h1fsu6r 
      complex*16  h1fsu7l,h1fsu7r,h1fsu8l,h1fsu8r
      real*8      ctt,cuu,css,ctu,cts,csu      
      integer     i,k,l
      integer     khb,kfer
      integer     ksfert(6),kchiu(4)
      integer     kfers(3)
      integer     nsfert,nchiu,nfers
      real*8 pi
      parameter (pi=3.141592653589793238d0)


c....check if the final state is kinematically allowed
      if((Svar-(mass(kp3)+mass(kp4))**2)/Svar.lt.thstep) then
        par=0.d0
        return
      endif      

***** set particles in final state:
      khb=kp3
      kfer=kp4

***** masses in final state:  
      mass3=mass(khb)  
      mass4=mass(kfer)  
***** define the kinematic variables  
      call dsaskinset2
   
c....symmetry factor, for non-identical final state particles
      s34=1.d0

c....mass symbols used in form-code  
      msfer2=mass(kp1)**2
      mchi2=mass(kp2)**2
      mhb2=mass(khb)**2
      mfer2=mass(kfer)**2
      s=Svar  
      t=Tvar  
      u=Uvar


***** now set the color factors to 1.d0 for a slepton in
***** the initial state and to 3.d0 for a squark in 
***** the initial state
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

***** initially setting the exchange array empty
      do i=1,6
         ksfert(i)=0
      enddo
      do i=1,3
         kfers(i)=0
      enddo

*****       
***** 
***** the first case    
***** sfermion(i) + neutralino(j) -> H^0_1/H^0_2 + fermion
***** the fermion must be a supersymmetric partner to the sfermion 
      if(icase.eq.1) then
***** set particles in the intermediate states 
***** in the t-channel (1 or 2 sfermions):  
        nsfert=ncsfert
        do i=1,nsfert
          ksfert(i)=kcsfertn(i)
        enddo
***** in the u-channel (4 neutralinos):
        do k=1,4
          kchiu(k)=kn(k)
        enddo
        nchiu=4
***** in the s-channel (1 fermion or none):
        if(dabs(echarg(kfer)).lt.1.d-16) then
          nfers=0       
        else  
          kfers(1)=kfer
          nfers=1
        endif
        goto 300
      endif
*****
***** the second case
***** sfermion(i) + neutralino(j) -> H^0_3 + fermion
***** the fermion must be a supersymmetric partner to the sfermion 
      if(icase.eq.2) then
***** set particles in the intermediate states 
***** in the t-channel (2 sfermions or none):  
        nsfert=ncsfert
        if(nsfert.eq.1) then
          nsfert=0
        else
          do i=1,nsfert
            ksfert(i)=kcsfertn(i)
          enddo
        endif
***** in the u-channel (4 neutralinos):
        do k=1,4
          kchiu(k)=kn(k)
        enddo
        nchiu=4
***** in the s-channel (1 fermion or none):
        if(dabs(echarg(kfer)).lt.1.d-16) then
          nfers=0       
        else  
          kfers(1)=kfer
          nfers=1
        endif
        goto 300
      endif
*****
***** the third case
***** down-type-sfermion(i) + chargino^+(j) 
***** -> H^0_1/H^0_2/H^0_3 + up-type-fermion
      if(icase.eq.3) then
***** set particles in the intermediate states 
***** in the t-channel (2 sfermions):  
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
***** in the s-channel (1 fermion or none):
        if(dabs(echarg(kfer)).lt.1.d-16) then
          nfers=0       
        else  
          kfers(1)=kfer
          nfers=1
        endif
        goto 300
      endif
*****
***** the fourth case
***** down-type-sfermion(i) + neutralino(j) 
***** -> H^- + up-type-fermion
      if(icase.eq.4) then
***** set particles in the intermediate states 
***** in the t-channel (1 or 2 sfermions):
        nsfert=ncsfertc
        do i=1,nsfert
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
        goto 300
      endif
*****
***** the fifth case     
***** up-type-sfermion(i) + chargino^+(j) 
***** -> H^+ + up-type-fermion
      if(icase.eq.5) then
***** set particles in the intermediate states 
***** in the t-channel (2 or 6 sfermions):  
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
        goto 300
      endif
***** 
***** the sixth case    
***** down-type-sfermion(i) + chargino^+(j) 
***** -> H^+ + down-type-fermion
      if(icase.eq.6) then
***** set particles in the intermediate states 
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
        goto 300
      endif
 
      write(*,*) 'DS: dsaschicasec called with wrong icase : ',icase   
      write(*,*) 'DS: initial or final states : '  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop  
      
 300  continue


c....specify the coefficients used in the form-code
c....first for the t-channel(sfermion exchange)
      h1ftt=dcmplx(0.d0,0.d0)
      h1ftt1l=dcmplx(0.d0,0.d0)
      h1ftt1r=dcmplx(0.d0,0.d0)
      h1ftt2l=dcmplx(0.d0,0.d0)
      h1ftt2r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
      h1ftt=1.d0
      do k=1,nsfert
        do l=1,nsfert 
         h1ftt1l=h1ftt1l+
     &     dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,ksfert(k),kp1)*conjg(gl(khb,ksfert(l),kp1))
     &     *gr(ksfert(k),kfer,kp2)*conjg(gr(ksfert(l),kfer,kp2))
         h1ftt1r=h1ftt1r+
     &     dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,ksfert(k),kp1)*conjg(gl(khb,ksfert(l),kp1))
     &     *gl(ksfert(k),kfer,kp2)*conjg(gl(ksfert(l),kfer,kp2))
         h1ftt2l=h1ftt2l+mass(kfer)*mass(kp2)
     &     *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,ksfert(k),kp1)*conjg(gl(khb,ksfert(l),kp1))
     &     *gl(ksfert(k),kfer,kp2)*conjg(gr(ksfert(l),kfer,kp2))
         h1ftt2r=h1ftt2r+mass(kfer)*mass(kp2)
     &     *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,ksfert(k),kp1)*conjg(gl(khb,ksfert(l),kp1))
     &     *gr(ksfert(k),kfer,kp2)*conjg(gl(ksfert(l),kfer,kp2))
        enddo
      enddo
      endif 


c....then the u-channel(neutralino/chargino exchange)
      h1fuu1l=dcmplx(0.d0,0.d0)
      h1fuu1r=dcmplx(0.d0,0.d0)
      h1fuu2l=dcmplx(0.d0,0.d0)
      h1fuu2r=dcmplx(0.d0,0.d0)
      h1fuu3l=dcmplx(0.d0,0.d0)
      h1fuu3r=dcmplx(0.d0,0.d0)
      h1fuu4l=dcmplx(0.d0,0.d0)
      h1fuu4r=dcmplx(0.d0,0.d0)
      h1fuu5l=dcmplx(0.d0,0.d0)
      h1fuu5r=dcmplx(0.d0,0.d0)
      h1fuu6l=dcmplx(0.d0,0.d0)
      h1fuu6r=dcmplx(0.d0,0.d0)
      h1fuu7l=dcmplx(0.d0,0.d0)
      h1fuu7r=dcmplx(0.d0,0.d0)
      h1fuu8l=dcmplx(0.d0,0.d0)
      h1fuu8r=dcmplx(0.d0,0.d0)
      if(nchiu.gt.0) then
      do k=1,nchiu
        do l=1,nchiu
         h1fuu1l=h1fuu1l+dsasdepro(u,kchiu(k))
     &    *conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu1r=h1fuu1r+dsasdepro(u,kchiu(k))
     &    *conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu2l=h1fuu2l+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu2r=h1fuu2r+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu3l=h1fuu3l+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu3r=h1fuu3r+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu4l=h1fuu4l+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu4r=h1fuu4r+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu5l=h1fuu5l+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu5r=h1fuu5r+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu6l=h1fuu6l+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l))) 
         h1fuu6r=h1fuu6r+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu7l=h1fuu7l+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu7r=h1fuu7r+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1fuu8l=h1fuu8l+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1fuu8r=h1fuu8r+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
        enddo
      enddo
      endif       


c....then the s-channel(fermion-exchange)
      h1fss=dcmplx(0.d0,0.d0) 
      h1fss1l=dcmplx(0.d0,0.d0)
      h1fss1r=dcmplx(0.d0,0.d0)
      h1fss2l=dcmplx(0.d0,0.d0)
      h1fss2r=dcmplx(0.d0,0.d0) 
      h1fss3l=dcmplx(0.d0,0.d0)
      h1fss3r=dcmplx(0.d0,0.d0)
      h1fss4l=dcmplx(0.d0,0.d0)
      h1fss4r=dcmplx(0.d0,0.d0) 
      h1fss5l=dcmplx(0.d0,0.d0)
      h1fss5r=dcmplx(0.d0,0.d0)
      h1fss6l=dcmplx(0.d0,0.d0)
      h1fss6r=dcmplx(0.d0,0.d0) 
      h1fss7l=dcmplx(0.d0,0.d0)
      h1fss7r=dcmplx(0.d0,0.d0)
      h1fss8l=dcmplx(0.d0,0.d0)
      h1fss8r=dcmplx(0.d0,0.d0)
      if(nfers.gt.0) then
***** the expression below has been changed Jan 2007
***** before we had only kfers now kfers(k), kfers(l)
      do k=1,nfers
        do l=1,nfers
         h1fss=dsasdepro(s,kfers(k))*conjg(dsasdepro(s,kfers(l)))
         h1fss1l=h1fss1l+h1fss*gr(khb,kfer,kfers(k))
     &          *gl(kp1,kfers(k),kp2)*conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l))) 
         h1fss1r=h1fss1r+h1fss*gl(khb,kfer,kfers(k))
     &          *gr(kp1,kfers(k),kp2)*conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l))) 
         h1fss2l=h1fss2l+h1fss*mass(kp2)*mass(kfers(l))
     &          *gr(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l))) 
         h1fss2r=h1fss2r+h1fss*mass(kp2)*mass(kfers(l))
     &          *gl(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))  
         h1fss3l=h1fss3l+h1fss*mass(kfers(k))*mass(kfers(l))
     &          *gr(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))  
         h1fss3r=h1fss3r+h1fss*mass(kfers(k))*mass(kfers(l))
     &          *gl(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l))) 
         h1fss4l=h1fss4l+h1fss*mass(kfers(k))*mass(kp2)
     &          *gr(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l))) 
         h1fss4r=h1fss4r+h1fss*mass(kfers(k))*mass(kp2)
     &          *gl(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))  
         h1fss5l=h1fss5l+h1fss*mass(kfer)*mass(kfers(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gl(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))  
         h1fss5r=h1fss5r+h1fss*mass(kfer)*mass(kfers(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gr(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fss6l=h1fss6l+h1fss*mass(kfer)*mass(kfers(k))
     &          *gl(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fss6r=h1fss6r+h1fss*mass(kfer)*mass(kfers(k))
     &          *gr(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l))) 
         h1fss7l=h1fss7l+h1fss*mass(kfer)*mass(kp2)
     &          *gl(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l))) 
         h1fss7r=h1fss7r+h1fss*mass(kfer)*mass(kp2)
     &          *gr(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l))) 
         h1fss8l=h1fss8l+h1fss*mass(kfer)*mass(kfers(l))
     &          *gl(khb,kfer,kfers(k))*gr(kp1,kfers(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l))) 
         h1fss8r=h1fss8r+h1fss*mass(kfer)*mass(kfers(l))
     &          *gr(khb,kfer,kfers(k))*gl(kp1,kfers(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l))) 
        enddo
      enddo
***** we must redefine the term h1fss, because it will be
***** multiplied on every ss term in the form expression.
***** that was ok before when just one exchange particle,
***** but not now with 3 exchange particles (k,l)
***** where it has to be introduced on the terms here above
      h1fss=dcmplx(1.d0,0.d0)
      endif


c....then t-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      h1ftu=dcmplx(0.d0,0.d0)
      h1ftu1l=dcmplx(0.d0,0.d0)
      h1ftu1r=dcmplx(0.d0,0.d0)
      h1ftu2l=dcmplx(0.d0,0.d0)
      h1ftu2r=dcmplx(0.d0,0.d0)
      h1ftu3l=dcmplx(0.d0,0.d0)
      h1ftu3r=dcmplx(0.d0,0.d0)
      h1ftu4l=dcmplx(0.d0,0.d0)
      h1ftu4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nchiu.gt.0) then
      h1ftu=1.d0
      do k=1,nsfert
        do l=1,nchiu
         h1ftu1l=h1ftu1l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kchiu(l))
     &    *gl(khb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1ftu1r=h1ftu1r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kchiu(l))
     &    *gl(khb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1ftu2l=h1ftu2l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kp2)
     &    *gl(khb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l)))
         h1ftu2r=h1ftu2r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kp2)
     &    *gl(khb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1ftu3l=h1ftu3l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kfer)
     &    *gl(khb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kfer,kchiu(l))) 
         h1ftu3r=h1ftu3r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kfer)
     &    *gl(khb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
         h1ftu4l=h1ftu4l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)*mass(kchiu(l))
     &    *gl(khb,ksfert(k),kp1)*gl(ksfert(k),kfer,kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))
     &    *conjg(gr(kp1,kfer,kchiu(l))) 
         h1ftu4r=h1ftu4r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)*mass(kchiu(l))
     &    *gl(khb,ksfert(k),kp1)*gr(ksfert(k),kfer,kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kfer,kchiu(l)))
        enddo
      enddo
      endif


c....t-channel amplitude multiplied by
c....hermitian conjugated s-channel amplitude 
      h1fts=dcmplx(0.d0,0.d0)
      h1fts1l=dcmplx(0.d0,0.d0)
      h1fts1r=dcmplx(0.d0,0.d0)
      h1fts2l=dcmplx(0.d0,0.d0)
      h1fts2r=dcmplx(0.d0,0.d0)
      h1fts3l=dcmplx(0.d0,0.d0)
      h1fts3r=dcmplx(0.d0,0.d0)
      h1fts4l=dcmplx(0.d0,0.d0)
      h1fts4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nfers.gt.0) then
*       h1fts=conjg(dsasdepro(s,kfers)) 
***** changed Jan 2007 to have kfers(l) not just kfers
***** must then include h1fts in all terms from the start
***** instead of in the form result, therefore
***** set h1fts=1 at the end
       do k=1,nsfert
         do l=1,nfers
         h1fts=conjg(dsasdepro(s,kfers(l))) 
         h1fts1l=h1fts1l+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfers(l))*gr(ksfert(k),kfer,kp2)
     &      *conjg(gr(kp1,kfers(l),kp2))*conjg(gr(khb,kfer,kfers(l)))
         h1fts1r=h1fts1r+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfers(l))*gl(ksfert(k),kfer,kp2)
     &      *conjg(gl(kp1,kfers(l),kp2))*conjg(gl(khb,kfer,kfers(l)))
         h1fts2l=h1fts2l+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kp2)*gr(ksfert(k),kfer,kp2)
     &      *conjg(gl(kp1,kfers(l),kp2))*conjg(gr(khb,kfer,kfers(l)))
         h1fts2r=h1fts2r+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kp2)*gl(ksfert(k),kfer,kp2)
     &      *conjg(gr(kp1,kfers(l),kp2))*conjg(gl(khb,kfer,kfers(l)))
         h1fts3l=h1fts3l+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfer)*gl(ksfert(k),kfer,kp2)
     &      *conjg(gl(kp1,kfers(l),kp2))*conjg(gr(khb,kfer,kfers(l)))
         h1fts3r=h1fts3r+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfer)*gr(ksfert(k),kfer,kp2)
     &      *conjg(gr(kp1,kfers(l),kp2))*conjg(gl(khb,kfer,kfers(l)))
         h1fts4l=h1fts4l+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfer)*mass(kp2)*mass(kfers(l))
     &      *gl(ksfert(k),kfer,kp2)
     &      *conjg(gr(kp1,kfers(l),kp2))*conjg(gr(khb,kfer,kfers(l)))
         h1fts4r=h1fts4r+h1fts*dsasdepro(t,ksfert(k))
     &      *gl(khb,ksfert(k),kp1)
     &      *mass(kfer)*mass(kp2)*mass(kfers(l))
     &      *gr(ksfert(k),kfer,kp2)
     &      *conjg(gl(kp1,kfers(l),kp2))*conjg(gl(khb,kfer,kfers(l)))
         enddo
       enddo
       h1fts=dcmplx(1.d0,0.d0)
      endif  


c....u-channel amplitude multiplied by
c....hermitian conjugated s-channel amplitude
      h1fsu=dcmplx(0.d0,0.d0) 
      h1fsu1l=dcmplx(0.d0,0.d0)
      h1fsu1r=dcmplx(0.d0,0.d0)
      h1fsu2l=dcmplx(0.d0,0.d0)
      h1fsu2r=dcmplx(0.d0,0.d0) 
      h1fsu3l=dcmplx(0.d0,0.d0)
      h1fsu3r=dcmplx(0.d0,0.d0)
      h1fsu4l=dcmplx(0.d0,0.d0)
      h1fsu4r=dcmplx(0.d0,0.d0)
      h1fsu5l=dcmplx(0.d0,0.d0)
      h1fsu5r=dcmplx(0.d0,0.d0)
      h1fsu6l=dcmplx(0.d0,0.d0)
      h1fsu6r=dcmplx(0.d0,0.d0)
      h1fsu7l=dcmplx(0.d0,0.d0)
      h1fsu7r=dcmplx(0.d0,0.d0)
      h1fsu8l=dcmplx(0.d0,0.d0)
      h1fsu8r=dcmplx(0.d0,0.d0) 
      if(nfers.gt.0.and.nchiu.gt.0) then
***** changed Jan 2007, comment see st interference term
*       h1fsu=conjg(dsasdepro(s,kfers)) 
       do k=1,nchiu 
         do l=1,nfers
         h1fsu=conjg(dsasdepro(s,kfers(l))) 
         h1fsu1l=h1fsu1l+h1fsu*dsasdepro(u,kchiu(k))
     &          *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu1r=h1fsu1r+h1fsu*dsasdepro(u,kchiu(k))
     &          *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu2l=h1fsu2l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu2r=h1fsu2r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu3l=h1fsu3l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kfers(l))
     &          *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu3r=h1fsu3r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kfers(l))
     &          *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu4l=h1fsu4l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kp2)
     &          *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu4r=h1fsu4r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kp2)
     &          *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu5l=h1fsu5l+h1fsu*dsasdepro(u,kchiu(k))
     &       *mass(kfer)*mass(kchiu(k))*mass(kp2)*mass(kfers(l))
     &       *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &       *conjg(gr(kp1,kfers(l),kp2))
     &       *conjg(gr(khb,kfer,kfers(l)))
         h1fsu5r=h1fsu5r+h1fsu*dsasdepro(u,kchiu(k))
     &       *mass(kfer)*mass(kchiu(k))*mass(kp2)*mass(kfers(l))
     &       *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &       *conjg(gl(kp1,kfers(l),kp2))
     &       *conjg(gl(khb,kfer,kfers(l)))
         h1fsu6l=h1fsu6l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))
     &          *gl(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu6r=h1fsu6r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))
     &          *gr(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu7l=h1fsu7l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kp2)
     &          *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu7r=h1fsu7r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kp2)
     &          *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         h1fsu8l=h1fsu8l+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kfers(l))
     &          *gl(kp1,kfer,kchiu(k))*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kfers(l),kp2))
     &          *conjg(gr(khb,kfer,kfers(l)))
         h1fsu8r=h1fsu8r+h1fsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kfers(l))
     &          *gr(kp1,kfer,kchiu(k))*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kfers(l),kp2))
     &          *conjg(gl(khb,kfer,kfers(l)))
         enddo
       enddo
       h1fsu=dcmplx(1.d0,0.d0)
      endif

***** After the template file follows the form expression of the 
***** ('complex') amplitude squared: dsaschicasecas
***** which is (M_tM_t^\dagger + M_sM_s^\dagger + M_uM_u^\dagger
***** + 2M_tM_u^\dagger + 2M_tM_s^\dagger + 2M_uM_s^\dagger)
***** As M_tM_u^\dagger+M_uM_t^\dagger=2Re(M_tM_u^\dagger) 
***** we have to take the real part of the interference terms,
***** in order to obtain the true amplitude squared.
***** This is most easily done by taking the real part of the
***** whole expression dsaschicasecas.
***** This is done at the very end of the code.
        
c....Template file for dsaschicasec ends here






      dsaschicasecas =
     &  + msfer2*mhb2*h1fuu1l*cuu
     &  + msfer2*mhb2*h1fuu1r*cuu
     &  + msfer2*mhb2*h1fss*h1fss1l*css
     &  + msfer2*mhb2*h1fss*h1fss1r*css
     &  - 2.d0*msfer2*mhb2*h1fsu*h1fsu1l*csu
     &  - 2.d0*msfer2*mhb2*h1fsu*h1fsu1r*csu
     &  + msfer2*mfer2*h1fuu1l*cuu
     &  + msfer2*mfer2*h1fuu1r*cuu
     &  - msfer2*mfer2*h1fss*h1fss1l*css
     &  - msfer2*mfer2*h1fss*h1fss1r*css
     &  - msfer2*s*h1fuu1l*cuu
     &  - msfer2*s*h1fuu1r*cuu
     &  - msfer2*s*h1fss*h1fss1l*css
     &  - msfer2*s*h1fss*h1fss1r*css
     &  + 2.d0*msfer2*s*h1fsu*h1fsu1l*csu
     &  + 2.d0*msfer2*s*h1fsu*h1fsu1r*csu
     &  + msfer2*h1fuu6l*cuu
     &  + msfer2*h1fuu6r*cuu
     &  + 2.d0*msfer2*h1fuu7l*cuu
     &  + 2.d0*msfer2*h1fuu7r*cuu
     &  + msfer2*h1fuu8l*cuu
     &  + msfer2*h1fuu8r*cuu
     &  - msfer2*h1fss*h1fss6l*css
     &  - msfer2*h1fss*h1fss6r*css
     &  - msfer2*h1fss*h1fss8l*css
     &  - msfer2*h1fss*h1fss8r*css
     &  + 2.d0*msfer2*h1ftu*h1ftu3l*ctu
     &  + 2.d0*msfer2*h1ftu*h1ftu3r*ctu
     &  - 2.d0*msfer2*h1fts*h1fts3l*cts
     &  - 2.d0*msfer2*h1fts*h1fts3r*cts
     &  - 2.d0*msfer2*h1fsu*h1fsu6l*csu
     &  - 2.d0*msfer2*h1fsu*h1fsu6r*csu
     &  - 2.d0*msfer2*h1fsu*h1fsu7l*csu
     &  - 2.d0*msfer2*h1fsu*h1fsu7r*csu
     &  + 2.d0*msfer2*h1fsu*h1fsu8l*csu
     &  + 2.d0*msfer2*h1fsu*h1fsu8r*csu
     &  + mchi2*mhb2*h1fuu1l*cuu
     &  + mchi2*mhb2*h1fuu1r*cuu
     &  - mchi2*mhb2*h1fss*h1fss1l*css
     &  - mchi2*mhb2*h1fss*h1fss1r*css
     &  + 3.d0*mchi2*mfer2*h1fuu1l*cuu
     &  + 3.d0*mchi2*mfer2*h1fuu1r*cuu
     &  + mchi2*mfer2*h1fss*h1fss1l*css
     &  + mchi2*mfer2*h1fss*h1fss1r*css
     &  + 2.d0*mchi2*mfer2*h1fsu*h1fsu1l*csu
     &  + 2.d0*mchi2*mfer2*h1fsu*h1fsu1r*csu
     &  - 2.d0*mchi2*s*h1fuu1l*cuu
     &  - 2.d0*mchi2*s*h1fuu1r*cuu
     &  + 2.d0*mchi2*s*h1fsu*h1fsu1l*csu
     &  + 2.d0*mchi2*s*h1fsu*h1fsu1r*csu
     &  - mchi2*t*h1fuu1l*cuu
     &  - mchi2*t*h1fuu1r*cuu
     &  + mchi2*h1ftt*h1ftt1l*ctt
     &  + mchi2*h1ftt*h1ftt1r*ctt
     &  + mchi2*h1fuu2l*cuu
     &  + mchi2*h1fuu2r*cuu
     &  + mchi2*h1fuu3l*cuu
     &  + mchi2*h1fuu3r*cuu
     &  + mchi2*h1fuu4l*cuu
     &  + mchi2*h1fuu4r*cuu
     &  + 2.d0*mchi2*h1fuu6l*cuu
     &  + 2.d0*mchi2*h1fuu6r*cuu
     &  + 2.d0*mchi2*h1fuu7l*cuu
     &  + 2.d0*mchi2*h1fuu7r*cuu
     &  + 2.d0*mchi2*h1fuu8l*cuu
     &  + 2.d0*mchi2*h1fuu8r*cuu
     &  + mchi2*h1fss*h1fss3l*css
     &  + mchi2*h1fss*h1fss3r*css
     &  + mchi2*h1fss*h1fss6l*css
     &  + mchi2*h1fss*h1fss6r*css
     &  + mchi2*h1fss*h1fss8l*css
     &  + mchi2*h1fss*h1fss8r*css
     &  + 2.d0*mchi2*h1ftu*h1ftu1l*ctu
     &  + 2.d0*mchi2*h1ftu*h1ftu1r*ctu
     &  + 2.d0*mchi2*h1ftu*h1ftu2l*ctu
     &  + 2.d0*mchi2*h1ftu*h1ftu2r*ctu
     &  + 4.d0*mchi2*h1ftu*h1ftu3l*ctu
     &  + 4.d0*mchi2*h1ftu*h1ftu3r*ctu
     &  + 2.d0*mchi2*h1fts*h1fts1l*cts
     &  + 2.d0*mchi2*h1fts*h1fts1r*cts
     &  + 2.d0*mchi2*h1fts*h1fts3l*cts
     &  + 2.d0*mchi2*h1fts*h1fts3r*cts
     &  + 2.d0*mchi2*h1fsu*h1fsu2l*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu2r*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu3l*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu3r*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu6l*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu6r*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu7l*csu
     &  + 2.d0*mchi2*h1fsu*h1fsu7r*csu
     &  + 4.d0*mchi2*h1fsu*h1fsu8l*csu
     &  + 4.d0*mchi2*h1fsu*h1fsu8r*csu
     &  + mchi2**2*h1fuu1l*cuu
     &  + mchi2**2*h1fuu1r*cuu
     &  - mhb2*s*h1fuu1l*cuu
     &  - mhb2*s*h1fuu1r*cuu
     &  - mhb2*s*h1fss*h1fss1l*css
     &  - mhb2*s*h1fss*h1fss1r*css
     &  + 2.d0*mhb2*s*h1fsu*h1fsu1l*csu
     &  + 2.d0*mhb2*s*h1fsu*h1fsu1r*csu
     &  + mhb2*h1fuu2l*cuu
     &  + mhb2*h1fuu2r*cuu
     &  + mhb2*h1fuu4l*cuu
     &  + mhb2*h1fuu4r*cuu
     &  + 2.d0*mhb2*h1fuu7l*cuu
     &  + 2.d0*mhb2*h1fuu7r*cuu
     &  - mhb2*h1fss*h1fss2l*css
     &  - mhb2*h1fss*h1fss2r*css
     &  - mhb2*h1fss*h1fss4l*css
     &  - mhb2*h1fss*h1fss4r*css
     &  + 2.d0*mhb2*h1ftu*h1ftu2l*ctu
     &  + 2.d0*mhb2*h1ftu*h1ftu2r*ctu
     &  - 2.d0*mhb2*h1fts*h1fts2l*cts
     &  - 2.d0*mhb2*h1fts*h1fts2r*cts
     &  + 2.d0*mhb2*h1fsu*h1fsu2l*csu
     &  + 2.d0*mhb2*h1fsu*h1fsu2r*csu
     &  - 2.d0*mhb2*h1fsu*h1fsu4l*csu
     &  - 2.d0*mhb2*h1fsu*h1fsu4r*csu
     &  - 2.d0*mhb2*h1fsu*h1fsu7l*csu
     &  - 2.d0*mhb2*h1fsu*h1fsu7r*csu
     &  - 2.d0*mfer2*s*h1fuu1l*cuu
     &  - 2.d0*mfer2*s*h1fuu1r*cuu
     &  + 2.d0*mfer2*s*h1fsu*h1fsu1l*csu
     &  + 2.d0*mfer2*s*h1fsu*h1fsu1r*csu
     &  - mfer2*t*h1fuu1l*cuu
     &  - mfer2*t*h1fuu1r*cuu
     &  + mfer2*h1ftt*h1ftt1l*ctt
     &  + mfer2*h1ftt*h1ftt1r*ctt
     &  + 2.d0*mfer2*h1fuu2l*cuu
     &  + 2.d0*mfer2*h1fuu2r*cuu
     &  + mfer2*h1fuu3l*cuu
     &  + mfer2*h1fuu3r*cuu
     &  + 2.d0*mfer2*h1fuu4l*cuu
     &  + 2.d0*mfer2*h1fuu4r*cuu
     &  + mfer2*h1fuu6l*cuu
     &  + mfer2*h1fuu6r*cuu
     &  + 2.d0*mfer2*h1fuu7l*cuu
     &  + 2.d0*mfer2*h1fuu7r*cuu
     &  + mfer2*h1fuu8l*cuu
     &  + mfer2*h1fuu8r*cuu
     &  + mfer2*h1fss*h1fss2l*css
     &  + mfer2*h1fss*h1fss2r*css
     &  + mfer2*h1fss*h1fss3l*css
     &  + mfer2*h1fss*h1fss3r*css
     &  + mfer2*h1fss*h1fss4l*css
     &  + mfer2*h1fss*h1fss4r*css
     &  + 2.d0*mfer2*h1ftu*h1ftu1l*ctu
     &  + 2.d0*mfer2*h1ftu*h1ftu1r*ctu
     &  + 4.d0*mfer2*h1ftu*h1ftu2l*ctu
     &  + 4.d0*mfer2*h1ftu*h1ftu2r*ctu
     &  + 2.d0*mfer2*h1ftu*h1ftu3l*ctu
     &  + 2.d0*mfer2*h1ftu*h1ftu3r*ctu
     &  + 2.d0*mfer2*h1fts*h1fts1l*cts
     &  + 2.d0*mfer2*h1fts*h1fts1r*cts
     &  + 2.d0*mfer2*h1fts*h1fts2l*cts
     &  + 2.d0*mfer2*h1fts*h1fts2r*cts
     &  + 4.d0*mfer2*h1fsu*h1fsu2l*csu
     &  + 4.d0*mfer2*h1fsu*h1fsu2r*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu3l*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu3r*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu4l*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu4r*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu7l*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu7r*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu8l*csu
     &  + 2.d0*mfer2*h1fsu*h1fsu8r*csu
     &  + mfer2**2*h1fuu1l*cuu
     &  + mfer2**2*h1fuu1r*cuu
     &  + s*t*h1fuu1l*cuu
     &  + s*t*h1fuu1r*cuu
     &  + s*t*h1fss*h1fss1l*css
     &  + s*t*h1fss*h1fss1r*css
     &  - 2.d0*s*t*h1fsu*h1fsu1l*csu
     &  - 2.d0*s*t*h1fsu*h1fsu1r*csu
     &  - s*h1fuu2l*cuu
     &  - s*h1fuu2r*cuu
     &  - s*h1fuu4l*cuu
     &  - s*h1fuu4r*cuu
     &  - s*h1fuu6l*cuu
     &  - s*h1fuu6r*cuu
     &  - 2.d0*s*h1fuu7l*cuu
     &  - 2.d0*s*h1fuu7r*cuu
     &  - s*h1fuu8l*cuu
     &  - s*h1fuu8r*cuu
     &  + s*h1fss*h1fss2l*css
     &  + s*h1fss*h1fss2r*css
     &  + s*h1fss*h1fss4l*css
     &  + s*h1fss*h1fss4r*css
     &  + s*h1fss*h1fss6l*css
     &  + s*h1fss*h1fss6r*css
     &  + 2.d0*s*h1fss*h1fss7l*css
     &  + 2.d0*s*h1fss*h1fss7r*css
     &  + s*h1fss*h1fss8l*css
     &  + s*h1fss*h1fss8r*css
     &  - 2.d0*s*h1ftu*h1ftu2l*ctu
     &  - 2.d0*s*h1ftu*h1ftu2r*ctu
     &  - 2.d0*s*h1ftu*h1ftu3l*ctu
     &  - 2.d0*s*h1ftu*h1ftu3r*ctu
     &  + 2.d0*s*h1fts*h1fts2l*cts
     &  + 2.d0*s*h1fts*h1fts2r*cts
     &  + 2.d0*s*h1fts*h1fts3l*cts
     &  + 2.d0*s*h1fts*h1fts3r*cts
     &  - 2.d0*s*h1fsu*h1fsu2l*csu
     &  - 2.d0*s*h1fsu*h1fsu2r*csu
     &  + 2.d0*s*h1fsu*h1fsu4l*csu
     &  + 2.d0*s*h1fsu*h1fsu4r*csu
     &  + 2.d0*s*h1fsu*h1fsu6l*csu
     &  + 2.d0*s*h1fsu*h1fsu6r*csu
     &  - 2.d0*s*h1fsu*h1fsu8l*csu
     &  - 2.d0*s*h1fsu*h1fsu8r*csu
     &  + s**2*h1fuu1l*cuu
     &  + s**2*h1fuu1r*cuu
     &  + s**2*h1fss*h1fss1l*css
     &  + s**2*h1fss*h1fss1r*css
     &  - 2.d0*s**2*h1fsu*h1fsu1l*csu
     &  - 2.d0*s**2*h1fsu*h1fsu1r*csu
     &  - t*h1ftt*h1ftt1l*ctt
     &  - t*h1ftt*h1ftt1r*ctt
     &  - t*h1fuu2l*cuu
     &  - t*h1fuu2r*cuu
     &  - t*h1fuu3l*cuu
     &  - t*h1fuu3r*cuu
     &  - t*h1fuu4l*cuu
     &  - t*h1fuu4r*cuu
     &  - t*h1fuu6l*cuu
     &  - t*h1fuu6r*cuu
     &  - 2.d0*t*h1fuu7l*cuu
     &  - 2.d0*t*h1fuu7r*cuu
     &  - t*h1fuu8l*cuu
     &  - t*h1fuu8r*cuu
     &  - t*h1fss*h1fss3l*css
     &  - t*h1fss*h1fss3r*css
     &  - 2.d0*t*h1ftu*h1ftu1l*ctu
     &  - 2.d0*t*h1ftu*h1ftu1r*ctu
     &  - 2.d0*t*h1ftu*h1ftu2l*ctu
     &  - 2.d0*t*h1ftu*h1ftu2r*ctu
     &  - 2.d0*t*h1ftu*h1ftu3l*ctu
     &  - 2.d0*t*h1ftu*h1ftu3r*ctu
     &  - 2.d0*t*h1fts*h1fts1l*cts
     &  - 2.d0*t*h1fts*h1fts1r*cts
     &  - 2.d0*t*h1fsu*h1fsu2l*csu
     &  - 2.d0*t*h1fsu*h1fsu2r*csu
     &  - 2.d0*t*h1fsu*h1fsu3l*csu
     &  - 2.d0*t*h1fsu*h1fsu3r*csu
     &  - 2.d0*t*h1fsu*h1fsu8l*csu
     &  - 2.d0*t*h1fsu*h1fsu8r*csu
     &  + 2.d0*h1ftt*h1ftt2l*ctt
     &  + 2.d0*h1ftt*h1ftt2r*ctt
     &  + 2.d0*h1fuu5l*cuu
     &  + 2.d0*h1fuu5r*cuu
     &  + 2.d0*h1fss*h1fss5l*css
     &  + 2.d0*h1fss*h1fss5r*css
     &  + 4.d0*h1ftu*h1ftu4l*ctu
     &  + 4.d0*h1ftu*h1ftu4r*ctu
     &  + 4.d0*h1fts*h1fts4l*cts
     &  + 4.d0*h1fts*h1fts4r*cts
     &  + 4.d0*h1fsu*h1fsu5l*csu
     &  + 4.d0*h1fsu*h1fsu5r*csu

      if (dble(dsaschicasecas).lt.0.0d0) then     
        if(aszeroprint) then
        write(*,*) ' '
        write(*,*) 
     &    'DS: ERROR IN dsaschicasec with negative cross section:'
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
        write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
        write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
        write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
        write(*,*) 'DS: dsaschicasecas = ',dsaschicasecas
        write(*,*) 'DS: mass1=',mass1
        write(*,*) 'DS: mass2=',mass2
        write(*,*) 'DS: mass3=',mass3
        write(*,*) 'DS: mass4=',mass4
        write(*,*) 'DS: s=',s
        write(*,*) 'DS: t=',t
        write(*,*) 'DS: u=',u
        write(*,*) 'DS: p12=',p12
        write(*,*) 'DS: costheta=',costheta
        write(*,*) 'DS: (s-(mass(khb)+mass(kfer))**2)/s',
     &    (s-(mass(khb)+mass(kfer))**2)/s
        endif
        dsaschicasecas=dcmplx(0.0d0,0.0d0)
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsaschicasecas)*k34/(8.0d0*pi*gg1c*gg2c*s34*dsqrt(s))
      return
      end

c...This subroutine is automatically generated from form output by
c...parsing it through form2f (version 1.34, October 8, 2001, edsjo@physto.se)
c....Template file for dsaschicased begins here

**************************************************************
*** SUBROUTINE dsaschicased                                ***
*** computes dW_{ij}/dcostheta                             ***
***                                                        ***
*** anti-sfermion(i) + neutralino(j)/chargino^+(j)         ***
*** -> higgs-boson + anti-fermion                          ***
***                                                        ***
*** The anti-sfermion must be the first mentioned          ***
*** particle (kp1) and the neutralino/chargino             ***
*** the other (kp2) -- not the opposite.                   ***
*** For the final state the higgs boson must be mentioned  ***
*** first (i.e. kp3) and next the anti-fermion (kp4) --    ***
*** not the opposite.                                      ***
***                                                        ***
***                                                        ***
*** Author:Mia Schelke, schelke@physto.se                  ***
*** Date: 01-10-11                                         ***
*** Trivial color factors included: 02-03-21               ***
*** Added flavour changing charged exchange:               ***
*** for the 1st, 2nd, 3rd, 4th and 5th case                ***
*** added by Mia Schelke January 2007                      ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

***** Note that it is assumed that coupling constants that do 
***** not exist have already been set to zero!!!!!
***** Thus, many of the coefficients defined in this code 
***** simplify when the diagrams contain sneutrinos 
***** or neutrinos.

      subroutine dsaschicased(kp1,kp2,kp3,kp4,icase,par)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsidtag.h'
      integer kp1,kp2,kp3,kp4,icase
      real*8 par,s,t,u
      complex*16 dsasdepro
      complex*16  dsaschicasedas
      real*8      msfer2,mchi2,mhb2,mfer2
      complex*16  h1fbtt,h1fbtt1l,h1fbtt1r,h1fbtt2l,h1fbtt2r
      complex*16  h1fbuu1l,h1fbuu1r,h1fbuu2l,h1fbuu2r
      complex*16  h1fbuu3l,h1fbuu3r,h1fbuu4l,h1fbuu4r
      complex*16  h1fbuu5l,h1fbuu5r,h1fbuu6l,h1fbuu6r
      complex*16  h1fbuu7l,h1fbuu7r,h1fbuu8l,h1fbuu8r
      complex*16  h1fbss,h1fbss1l,h1fbss1r,h1fbss2l
      complex*16  h1fbss2r,h1fbss3l,h1fbss3r,h1fbss4l
      complex*16  h1fbss4r,h1fbss5l,h1fbss5r,h1fbss6l
      complex*16  h1fbss6r,h1fbss7l,h1fbss7r,h1fbss8l,h1fbss8r
      complex*16  h1fbtu,h1fbtu1l,h1fbtu1r,h1fbtu2l,h1fbtu2r
      complex*16  h1fbtu3l,h1fbtu3r,h1fbtu4l,h1fbtu4r
      complex*16  h1fbts,h1fbts1l,h1fbts1r,h1fbts2l,h1fbts2r
      complex*16  h1fbts3l,h1fbts3r,h1fbts4l,h1fbts4r
      complex*16  h1fbsu,h1fbsu1l,h1fbsu1r,h1fbsu2l,h1fbsu2r 
      complex*16  h1fbsu3l,h1fbsu3r,h1fbsu4l,h1fbsu4r 
      complex*16  h1fbsu5l,h1fbsu5r,h1fbsu6l,h1fbsu6r 
      complex*16  h1fbsu7l,h1fbsu7r,h1fbsu8l,h1fbsu8r
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
***** anti-up-type-sfermion(i) + chargino^+(j) 
***** -> H^0_1/H^0_2 + anti-down-type-fermion
      if(icase.eq.1) then
***** set particles in the intermediate states 
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
        goto 400
      endif        
*****
***** the second case
***** anti-up-type-sfermion(i) + chargino^+(j) 
***** -> H^0_3 + anti-down-type-fermion
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
***** in the u-channel (2 charginos):
        do k=1,2
          kchiu(k)=kcha(k)
          kchiu(k+2)=0
        enddo
        nchiu=2 
***** in the s-channel (1 fermion):
        kfers(1)=kfer
        nfers=1
        goto 400
      endif  
*****
***** the third case
***** anti-up-type-sfermion(i) + neutralino(j) 
***** -> H^- + anti-down-type-fermion
      if(icase.eq.3) then
***** set particles in the intermediate states 
***** in the t-channel (2 sfermions):
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
        goto 400
      endif  
*****
***** the fourth case
***** anti-up-type-sfermion(i) + chargino^+(j) 
***** -> H^+ + anti-up-type-fermion
      if(icase.eq.4) then
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
***** in the s-channel (1 or 3 fermions):
        nfers=ncferd
        do i=1,nfers
          kfers(i)=kcferd(i)
        enddo
        goto 400
      endif  
*****
***** the fifth case     
***** anti-down-type-sfermion(i) + chargino^+(j) 
***** -> H^+ + anti-down-type-fermion
      if(icase.eq.5) then
***** set particles in the intermediate states 
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
        goto 400
      endif         
      
      write(*,*) 'DS: dsaschicased called with wrong icase : ',icase   
      write(*,*) 'DS: initial or final states : '  
      write(*,*)   pname(kp1),pname(kp2),' -> ',
     &             pname(kp3),pname(kp4)   
      stop

 400  continue

c....specify the coefficients used in the form-code
c....
c....note that when writing down the Feynman rules for the diagrams 
c....we have started with the anti-fermion instead of ending there
c....(the two results are identical, which can be seen by 
c....transposing one of them)
c....in this way we got a result that looked very similar to the 
c....one in dsaschicasec
c....
c....first for the t-channel(sfermion exchange)
      h1fbtt=dcmplx(0.d0,0.d0)
      h1fbtt1l=dcmplx(0.d0,0.d0)
      h1fbtt1r=dcmplx(0.d0,0.d0)
      h1fbtt2l=dcmplx(0.d0,0.d0)
      h1fbtt2r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0) then
      h1fbtt=1.d0
      do k=1,nsfert
        do l=1,nsfert 
         h1fbtt1l=h1fbtt1l+
     &     dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,kp1,ksfert(k))*conjg(gl(khb,kp1,ksfert(l)))
     &     *gr(ksfert(k),kp2,kfer)*conjg(gr(ksfert(l),kp2,kfer))
         h1fbtt1r=h1fbtt1r+
     &     dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,kp1,ksfert(k))*conjg(gl(khb,kp1,ksfert(l)))
     &     *gl(ksfert(k),kp2,kfer)*conjg(gl(ksfert(l),kp2,kfer))
         h1fbtt2l=h1fbtt2l+mass(kfer)*mass(kp2)
     &     *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,kp1,ksfert(k))*conjg(gl(khb,kp1,ksfert(l)))
     &     *gl(ksfert(k),kp2,kfer)*conjg(gr(ksfert(l),kp2,kfer))
         h1fbtt2r=h1fbtt2r+mass(kfer)*mass(kp2)
     &     *dsasdepro(t,ksfert(k))*conjg(dsasdepro(t,ksfert(l)))
     &     *gl(khb,kp1,ksfert(k))*conjg(gl(khb,kp1,ksfert(l)))
     &     *gr(ksfert(k),kp2,kfer)*conjg(gl(ksfert(l),kp2,kfer))
        enddo
      enddo
      endif 


c....then the u-channel(neutralino/chargino exchange)
      h1fbuu1l=dcmplx(0.d0,0.d0)
      h1fbuu1r=dcmplx(0.d0,0.d0)
      h1fbuu2l=dcmplx(0.d0,0.d0)
      h1fbuu2r=dcmplx(0.d0,0.d0)
      h1fbuu3l=dcmplx(0.d0,0.d0)
      h1fbuu3r=dcmplx(0.d0,0.d0)
      h1fbuu4l=dcmplx(0.d0,0.d0)
      h1fbuu4r=dcmplx(0.d0,0.d0)
      h1fbuu5l=dcmplx(0.d0,0.d0)
      h1fbuu5r=dcmplx(0.d0,0.d0)
      h1fbuu6l=dcmplx(0.d0,0.d0)
      h1fbuu6r=dcmplx(0.d0,0.d0)
      h1fbuu7l=dcmplx(0.d0,0.d0)
      h1fbuu7r=dcmplx(0.d0,0.d0)
      h1fbuu8l=dcmplx(0.d0,0.d0)
      h1fbuu8r=dcmplx(0.d0,0.d0)
      if(nchiu.gt.0) then
      do k=1,nchiu
        do l=1,nchiu
         h1fbuu1l=h1fbuu1l
     &    +dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu1r=h1fbuu1r
     &    +dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu2l=h1fbuu2l+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu2r=h1fbuu2r+mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu3l=h1fbuu3l+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu3r=h1fbuu3r+mass(kchiu(k))*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu4l=h1fbuu4l+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu4r=h1fbuu4r+mass(kchiu(k))*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu5l=h1fbuu5l+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu5r=h1fbuu5r+mass(kfer)*mass(kchiu(k))
     &    *mass(kp2)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu6l=h1fbuu6l+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         h1fbuu6r=h1fbuu6r+mass(kfer)*mass(kchiu(k))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu7l=h1fbuu7l+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu7r=h1fbuu7r+mass(kfer)*mass(kp2)
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbuu8l=h1fbuu8l+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbuu8r=h1fbuu8r+mass(kfer)*mass(kchiu(l))
     &    *dsasdepro(u,kchiu(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
        enddo
      enddo
      endif       


c....then the s-channel(fermion-exchange)
      h1fbss=dcmplx(0.d0,0.d0) 
      h1fbss1l=dcmplx(0.d0,0.d0)
      h1fbss1r=dcmplx(0.d0,0.d0)
      h1fbss2l=dcmplx(0.d0,0.d0)
      h1fbss2r=dcmplx(0.d0,0.d0) 
      h1fbss3l=dcmplx(0.d0,0.d0)
      h1fbss3r=dcmplx(0.d0,0.d0)
      h1fbss4l=dcmplx(0.d0,0.d0)
      h1fbss4r=dcmplx(0.d0,0.d0) 
      h1fbss5l=dcmplx(0.d0,0.d0)
      h1fbss5r=dcmplx(0.d0,0.d0)
      h1fbss6l=dcmplx(0.d0,0.d0)
      h1fbss6r=dcmplx(0.d0,0.d0) 
      h1fbss7l=dcmplx(0.d0,0.d0)
      h1fbss7r=dcmplx(0.d0,0.d0)
      h1fbss8l=dcmplx(0.d0,0.d0)
      h1fbss8r=dcmplx(0.d0,0.d0)
      if(nfers.gt.0) then
***** the expression below has been changed Jan 2007
***** before we had only kfers now kfers(k), kfers(l)
      do k=1,nfers
        do l=1,nfers
         h1fbss=dsasdepro(s,kfers(k))*conjg(dsasdepro(s,kfers(l)))
         h1fbss1l=h1fbss1l+h1fbss*gr(khb,kfers(k),kfer)
     &          *gl(kp1,kp2,kfers(k))*conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer)) 
         h1fbss1r=h1fbss1r+h1fbss*gl(khb,kfers(k),kfer)
     &          *gr(kp1,kp2,kfers(k))*conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer)) 
         h1fbss2l=h1fbss2l+h1fbss*mass(kp2)*mass(kfers(l))
     &          *gr(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer)) 
         h1fbss2r=h1fbss2r+h1fbss*mass(kp2)*mass(kfers(l))
     &          *gl(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))  
         h1fbss3l=h1fbss3l+h1fbss*mass(kfers(k))*mass(kfers(l))
     &          *gr(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))  
         h1fbss3r=h1fbss3r+h1fbss*mass(kfers(k))*mass(kfers(l))
     &          *gl(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer)) 
         h1fbss4l=h1fbss4l+h1fbss*mass(kfers(k))*mass(kp2)
     &          *gr(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer)) 
         h1fbss4r=h1fbss4r+h1fbss*mass(kfers(k))*mass(kp2)
     &          *gl(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))  
         h1fbss5l=h1fbss5l+h1fbss
     &          *mass(kfer)*mass(kfers(k))*mass(kp2)*mass(kfers(l))
     &          *gl(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))  
         h1fbss5r=h1fbss5r+h1fbss
     &          *mass(kfer)*mass(kfers(k))*mass(kp2)*mass(kfers(l))
     &          *gr(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
         h1fbss6l=h1fbss6l+h1fbss*mass(kfer)*mass(kfers(k))
     &          *gl(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))
         h1fbss6r=h1fbss6r+h1fbss*mass(kfer)*mass(kfers(k))
     &          *gr(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer)) 
         h1fbss7l=h1fbss7l+h1fbss*mass(kfer)*mass(kp2)
     &          *gl(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer)) 
         h1fbss7r=h1fbss7r+h1fbss*mass(kfer)*mass(kp2)
     &          *gr(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer)) 
         h1fbss8l=h1fbss8l+h1fbss*mass(kfer)*mass(kfers(l))
     &          *gl(khb,kfers(k),kfer)*gr(kp1,kp2,kfers(k))
     &      *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer)) 
         h1fbss8r=h1fbss8r+h1fbss*mass(kfer)*mass(kfers(l))
     &          *gr(khb,kfers(k),kfer)*gl(kp1,kp2,kfers(k))
     &      *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
        enddo
      enddo
***** we must redefine the term h1fbss, because it will be
***** multiplied on every ss term in the form expression.
***** that was ok before when just one exchange particle,
***** but not now with 3 exchange particles (k,l)
***** where it has to be introduced on the terms here above
      h1fbss=dcmplx(1.d0,0.d0)
      endif


c....then t-channel amplitude multiplied by
c....hermitian conjugated u-channel amplitude
      h1fbtu=dcmplx(0.d0,0.d0)
      h1fbtu1l=dcmplx(0.d0,0.d0)
      h1fbtu1r=dcmplx(0.d0,0.d0)
      h1fbtu2l=dcmplx(0.d0,0.d0)
      h1fbtu2r=dcmplx(0.d0,0.d0)
      h1fbtu3l=dcmplx(0.d0,0.d0)
      h1fbtu3r=dcmplx(0.d0,0.d0)
      h1fbtu4l=dcmplx(0.d0,0.d0)
      h1fbtu4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nchiu.gt.0) then
      h1fbtu=1.d0
      do k=1,nsfert
        do l=1,nchiu
         h1fbtu1l=h1fbtu1l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kchiu(l))
     &    *gl(khb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbtu1r=h1fbtu1r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kchiu(l))
     &    *gl(khb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbtu2l=h1fbtu2l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kp2)
     &    *gl(khb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer))
         h1fbtu2r=h1fbtu2r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kp2)
     &    *gl(khb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbtu3l=h1fbtu3l+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kfer)
     &    *gl(khb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         h1fbtu3r=h1fbtu3r+dsasdepro(t,ksfert(k))
     &    *conjg(dsasdepro(u,kchiu(l)))*mass(kfer)
     &    *gl(khb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
         h1fbtu4l=h1fbtu4l+
     &    dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)*mass(kchiu(l))
     &    *gl(khb,kp1,ksfert(k))*gl(ksfert(k),kp2,kfer)
     &    *conjg(gr(khb,kchiu(l),kp2))*conjg(gr(kp1,kchiu(l),kfer)) 
         h1fbtu4r=h1fbtu4r+
     &    dsasdepro(t,ksfert(k))*conjg(dsasdepro(u,kchiu(l)))
     &    *mass(kfer)*mass(kp2)*mass(kchiu(l))
     &    *gl(khb,kp1,ksfert(k))*gr(ksfert(k),kp2,kfer)
     &    *conjg(gl(khb,kchiu(l),kp2))*conjg(gl(kp1,kchiu(l),kfer))
        enddo
      enddo
      endif


c....t-channel amplitude multiplied by
c....hermitian conjugated s-channel amplitude 
      h1fbts=dcmplx(0.d0,0.d0)
      h1fbts1l=dcmplx(0.d0,0.d0)
      h1fbts1r=dcmplx(0.d0,0.d0)
      h1fbts2l=dcmplx(0.d0,0.d0)
      h1fbts2r=dcmplx(0.d0,0.d0)
      h1fbts3l=dcmplx(0.d0,0.d0)
      h1fbts3r=dcmplx(0.d0,0.d0)
      h1fbts4l=dcmplx(0.d0,0.d0)
      h1fbts4r=dcmplx(0.d0,0.d0)
      if(nsfert.gt.0.and.nfers.gt.0) then
***** below changed Jan 2007 to have kfers(l)
***** must then include h1fbts in all terms from the start
***** instead of in the form result, therefore
***** set h1fbts=1 at the end
*       h1fbts=conjg(dsasdepro(s,kfers)) 
       do k=1,nsfert
         do l=1,nfers
         h1fbts=conjg(dsasdepro(s,kfers(l)))
         h1fbts1l=h1fbts1l+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfers(l))*gr(ksfert(k),kp2,kfer)
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))
         h1fbts1r=h1fbts1r+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfers(l))*gl(ksfert(k),kp2,kfer)
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
         h1fbts2l=h1fbts2l+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kp2)*gr(ksfert(k),kp2,kfer)
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))
         h1fbts2r=h1fbts2r+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kp2)*gl(ksfert(k),kp2,kfer)
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
         h1fbts3l=h1fbts3l+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfer)*gl(ksfert(k),kp2,kfer)
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))
         h1fbts3r=h1fbts3r+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfer)*gr(ksfert(k),kp2,kfer)
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
         h1fbts4l=h1fbts4l+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfer)*mass(kp2)*mass(kfers(l))
     &         *gl(ksfert(k),kp2,kfer)
     &     *conjg(gr(kp1,kp2,kfers(l)))*conjg(gr(khb,kfers(l),kfer))
         h1fbts4r=h1fbts4r+h1fbts
     &         *dsasdepro(t,ksfert(k))*gl(khb,kp1,ksfert(k))
     &         *mass(kfer)*mass(kp2)*mass(kfers(l))
     &         *gr(ksfert(k),kp2,kfer)
     &     *conjg(gl(kp1,kp2,kfers(l)))*conjg(gl(khb,kfers(l),kfer))
         enddo
       enddo
       h1fbts=dcmplx(1.d0,0.d0)
      endif  


c....u-channel amplitude multiplied by
c....hermitian conjugated s-channel amplitude
      h1fbsu=dcmplx(0.d0,0.d0) 
      h1fbsu1l=dcmplx(0.d0,0.d0)
      h1fbsu1r=dcmplx(0.d0,0.d0)
      h1fbsu2l=dcmplx(0.d0,0.d0)
      h1fbsu2r=dcmplx(0.d0,0.d0) 
      h1fbsu3l=dcmplx(0.d0,0.d0)
      h1fbsu3r=dcmplx(0.d0,0.d0)
      h1fbsu4l=dcmplx(0.d0,0.d0)
      h1fbsu4r=dcmplx(0.d0,0.d0)
      h1fbsu5l=dcmplx(0.d0,0.d0)
      h1fbsu5r=dcmplx(0.d0,0.d0)
      h1fbsu6l=dcmplx(0.d0,0.d0)
      h1fbsu6r=dcmplx(0.d0,0.d0)
      h1fbsu7l=dcmplx(0.d0,0.d0)
      h1fbsu7r=dcmplx(0.d0,0.d0)
      h1fbsu8l=dcmplx(0.d0,0.d0)
      h1fbsu8r=dcmplx(0.d0,0.d0) 
      if(nfers.gt.0.and.nchiu.gt.0) then
***** below changed Jan 2007 to have kfers(l)
***** must then include h1fbsu in all terms from the start
***** instead of in the form result, therefore
***** set h1fbsu=1 at the end
*       h1fbsu=conjg(dsasdepro(s,kfers)) 
       do k=1,nchiu 
         do l=1,nfers 
         h1fbsu=conjg(dsasdepro(s,kfers(l))) 
         h1fbsu1l=h1fbsu1l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu1r=h1fbsu1r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu2l=h1fbsu2l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu2r=h1fbsu2r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kp2)*mass(kfers(l))
     &          *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu3l=h1fbsu3l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kfers(l))
     &          *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu3r=h1fbsu3r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kfers(l))
     &          *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu4l=h1fbsu4l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kp2)
     &          *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu4r=h1fbsu4r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kchiu(k))*mass(kp2)
     &          *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu5l=h1fbsu5l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))*mass(kp2)*mass(kfers(l))
     &          *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu5r=h1fbsu5r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))*mass(kp2)*mass(kfers(l))
     &          *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu6l=h1fbsu6l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))
     &          *gl(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu6r=h1fbsu6r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kchiu(k))
     &          *gr(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu7l=h1fbsu7l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kp2)
     &          *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu7r=h1fbsu7r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kp2)
     &          *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         h1fbsu8l=h1fbsu8l+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kfers(l))
     &          *gl(kp1,kchiu(k),kfer)*gr(khb,kchiu(k),kp2)
     &          *conjg(gr(kp1,kp2,kfers(l)))
     &          *conjg(gr(khb,kfers(l),kfer))
         h1fbsu8r=h1fbsu8r+h1fbsu*dsasdepro(u,kchiu(k))
     &          *mass(kfer)*mass(kfers(l))
     &          *gr(kp1,kchiu(k),kfer)*gl(khb,kchiu(k),kp2)
     &          *conjg(gl(kp1,kp2,kfers(l)))
     &          *conjg(gl(khb,kfers(l),kfer))
         enddo
       enddo
       h1fbsu=dcmplx(1.d0,0.d0)
      endif

***** After the template file follows the form expression of the 
***** ('complex') amplitude squared: dsaschicasedas
***** which is (M_tM_t^\dagger + M_sM_s^\dagger + M_uM_u^\dagger
***** + 2M_tM_u^\dagger + 2M_tM_s^\dagger + 2M_uM_s^\dagger)
***** As M_tM_u^\dagger+M_uM_t^\dagger=2Re(M_tM_u^\dagger) 
***** we have to take the real part of the interference terms,
***** in order to obtain the true amplitude squared.
***** This is most easily done by taking the real part of the
***** whole expression dsaschicasedas.
***** This is done at the very end of the code.
        
c....Template file for dsaschicased ends here






      dsaschicasedas =
     &  + msfer2*mhb2*h1fbuu1l*cuu
     &  + msfer2*mhb2*h1fbuu1r*cuu
     &  + msfer2*mhb2*h1fbss*h1fbss1l*css
     &  + msfer2*mhb2*h1fbss*h1fbss1r*css
     &  - 2.d0*msfer2*mhb2*h1fbsu*h1fbsu1l*csu
     &  - 2.d0*msfer2*mhb2*h1fbsu*h1fbsu1r*csu
     &  + msfer2*mfer2*h1fbuu1l*cuu
     &  + msfer2*mfer2*h1fbuu1r*cuu
     &  - msfer2*mfer2*h1fbss*h1fbss1l*css
     &  - msfer2*mfer2*h1fbss*h1fbss1r*css
     &  - msfer2*s*h1fbuu1l*cuu
     &  - msfer2*s*h1fbuu1r*cuu
     &  - msfer2*s*h1fbss*h1fbss1l*css
     &  - msfer2*s*h1fbss*h1fbss1r*css
     &  + 2.d0*msfer2*s*h1fbsu*h1fbsu1l*csu
     &  + 2.d0*msfer2*s*h1fbsu*h1fbsu1r*csu
     &  + msfer2*h1fbuu6l*cuu
     &  + msfer2*h1fbuu6r*cuu
     &  + 2.d0*msfer2*h1fbuu7l*cuu
     &  + 2.d0*msfer2*h1fbuu7r*cuu
     &  + msfer2*h1fbuu8l*cuu
     &  + msfer2*h1fbuu8r*cuu
     &  - msfer2*h1fbss*h1fbss6l*css
     &  - msfer2*h1fbss*h1fbss6r*css
     &  - msfer2*h1fbss*h1fbss8l*css
     &  - msfer2*h1fbss*h1fbss8r*css
     &  + 2.d0*msfer2*h1fbtu*h1fbtu3l*ctu
     &  + 2.d0*msfer2*h1fbtu*h1fbtu3r*ctu
     &  - 2.d0*msfer2*h1fbts*h1fbts3l*cts
     &  - 2.d0*msfer2*h1fbts*h1fbts3r*cts
     &  - 2.d0*msfer2*h1fbsu*h1fbsu6l*csu
     &  - 2.d0*msfer2*h1fbsu*h1fbsu6r*csu
     &  - 2.d0*msfer2*h1fbsu*h1fbsu7l*csu
     &  - 2.d0*msfer2*h1fbsu*h1fbsu7r*csu
     &  + 2.d0*msfer2*h1fbsu*h1fbsu8l*csu
     &  + 2.d0*msfer2*h1fbsu*h1fbsu8r*csu
     &  + mchi2*mhb2*h1fbuu1l*cuu
     &  + mchi2*mhb2*h1fbuu1r*cuu
     &  - mchi2*mhb2*h1fbss*h1fbss1l*css
     &  - mchi2*mhb2*h1fbss*h1fbss1r*css
     &  + 3.d0*mchi2*mfer2*h1fbuu1l*cuu
     &  + 3.d0*mchi2*mfer2*h1fbuu1r*cuu
     &  + mchi2*mfer2*h1fbss*h1fbss1l*css
     &  + mchi2*mfer2*h1fbss*h1fbss1r*css
     &  + 2.d0*mchi2*mfer2*h1fbsu*h1fbsu1l*csu
     &  + 2.d0*mchi2*mfer2*h1fbsu*h1fbsu1r*csu
     &  - 2.d0*mchi2*s*h1fbuu1l*cuu
     &  - 2.d0*mchi2*s*h1fbuu1r*cuu
     &  + 2.d0*mchi2*s*h1fbsu*h1fbsu1l*csu
     &  + 2.d0*mchi2*s*h1fbsu*h1fbsu1r*csu
     &  - mchi2*t*h1fbuu1l*cuu
     &  - mchi2*t*h1fbuu1r*cuu
     &  + mchi2*h1fbtt*h1fbtt1l*ctt
     &  + mchi2*h1fbtt*h1fbtt1r*ctt
     &  + mchi2*h1fbuu2l*cuu
     &  + mchi2*h1fbuu2r*cuu
     &  + mchi2*h1fbuu3l*cuu
     &  + mchi2*h1fbuu3r*cuu
     &  + mchi2*h1fbuu4l*cuu
     &  + mchi2*h1fbuu4r*cuu
     &  + 2.d0*mchi2*h1fbuu6l*cuu
     &  + 2.d0*mchi2*h1fbuu6r*cuu
     &  + 2.d0*mchi2*h1fbuu7l*cuu
     &  + 2.d0*mchi2*h1fbuu7r*cuu
     &  + 2.d0*mchi2*h1fbuu8l*cuu
     &  + 2.d0*mchi2*h1fbuu8r*cuu
     &  + mchi2*h1fbss*h1fbss3l*css
     &  + mchi2*h1fbss*h1fbss3r*css
     &  + mchi2*h1fbss*h1fbss6l*css
     &  + mchi2*h1fbss*h1fbss6r*css
     &  + mchi2*h1fbss*h1fbss8l*css
     &  + mchi2*h1fbss*h1fbss8r*css
     &  + 2.d0*mchi2*h1fbtu*h1fbtu1l*ctu
     &  + 2.d0*mchi2*h1fbtu*h1fbtu1r*ctu
     &  + 2.d0*mchi2*h1fbtu*h1fbtu2l*ctu
     &  + 2.d0*mchi2*h1fbtu*h1fbtu2r*ctu
     &  + 4.d0*mchi2*h1fbtu*h1fbtu3l*ctu
     &  + 4.d0*mchi2*h1fbtu*h1fbtu3r*ctu
     &  + 2.d0*mchi2*h1fbts*h1fbts1l*cts
     &  + 2.d0*mchi2*h1fbts*h1fbts1r*cts
     &  + 2.d0*mchi2*h1fbts*h1fbts3l*cts
     &  + 2.d0*mchi2*h1fbts*h1fbts3r*cts
     &  + 2.d0*mchi2*h1fbsu*h1fbsu2l*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu2r*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu3l*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu3r*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu6l*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu6r*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu7l*csu
     &  + 2.d0*mchi2*h1fbsu*h1fbsu7r*csu
     &  + 4.d0*mchi2*h1fbsu*h1fbsu8l*csu
     &  + 4.d0*mchi2*h1fbsu*h1fbsu8r*csu
     &  + mchi2**2*h1fbuu1l*cuu
     &  + mchi2**2*h1fbuu1r*cuu
     &  - mhb2*s*h1fbuu1l*cuu
     &  - mhb2*s*h1fbuu1r*cuu
     &  - mhb2*s*h1fbss*h1fbss1l*css
     &  - mhb2*s*h1fbss*h1fbss1r*css
     &  + 2.d0*mhb2*s*h1fbsu*h1fbsu1l*csu
     &  + 2.d0*mhb2*s*h1fbsu*h1fbsu1r*csu
     &  + mhb2*h1fbuu2l*cuu
     &  + mhb2*h1fbuu2r*cuu
     &  + mhb2*h1fbuu4l*cuu
     &  + mhb2*h1fbuu4r*cuu
     &  + 2.d0*mhb2*h1fbuu7l*cuu
     &  + 2.d0*mhb2*h1fbuu7r*cuu
     &  - mhb2*h1fbss*h1fbss2l*css
     &  - mhb2*h1fbss*h1fbss2r*css
     &  - mhb2*h1fbss*h1fbss4l*css
     &  - mhb2*h1fbss*h1fbss4r*css
     &  + 2.d0*mhb2*h1fbtu*h1fbtu2l*ctu
     &  + 2.d0*mhb2*h1fbtu*h1fbtu2r*ctu
     &  - 2.d0*mhb2*h1fbts*h1fbts2l*cts
     &  - 2.d0*mhb2*h1fbts*h1fbts2r*cts
     &  + 2.d0*mhb2*h1fbsu*h1fbsu2l*csu
     &  + 2.d0*mhb2*h1fbsu*h1fbsu2r*csu
     &  - 2.d0*mhb2*h1fbsu*h1fbsu4l*csu
     &  - 2.d0*mhb2*h1fbsu*h1fbsu4r*csu
     &  - 2.d0*mhb2*h1fbsu*h1fbsu7l*csu
     &  - 2.d0*mhb2*h1fbsu*h1fbsu7r*csu
     &  - 2.d0*mfer2*s*h1fbuu1l*cuu
     &  - 2.d0*mfer2*s*h1fbuu1r*cuu
     &  + 2.d0*mfer2*s*h1fbsu*h1fbsu1l*csu
     &  + 2.d0*mfer2*s*h1fbsu*h1fbsu1r*csu
     &  - mfer2*t*h1fbuu1l*cuu
     &  - mfer2*t*h1fbuu1r*cuu
     &  + mfer2*h1fbtt*h1fbtt1l*ctt
     &  + mfer2*h1fbtt*h1fbtt1r*ctt
     &  + 2.d0*mfer2*h1fbuu2l*cuu
     &  + 2.d0*mfer2*h1fbuu2r*cuu
     &  + mfer2*h1fbuu3l*cuu
     &  + mfer2*h1fbuu3r*cuu
     &  + 2.d0*mfer2*h1fbuu4l*cuu
     &  + 2.d0*mfer2*h1fbuu4r*cuu
     &  + mfer2*h1fbuu6l*cuu
     &  + mfer2*h1fbuu6r*cuu
     &  + 2.d0*mfer2*h1fbuu7l*cuu
     &  + 2.d0*mfer2*h1fbuu7r*cuu
     &  + mfer2*h1fbuu8l*cuu
     &  + mfer2*h1fbuu8r*cuu
     &  + mfer2*h1fbss*h1fbss2l*css
     &  + mfer2*h1fbss*h1fbss2r*css
     &  + mfer2*h1fbss*h1fbss3l*css
     &  + mfer2*h1fbss*h1fbss3r*css
     &  + mfer2*h1fbss*h1fbss4l*css
     &  + mfer2*h1fbss*h1fbss4r*css
     &  + 2.d0*mfer2*h1fbtu*h1fbtu1l*ctu
     &  + 2.d0*mfer2*h1fbtu*h1fbtu1r*ctu
     &  + 4.d0*mfer2*h1fbtu*h1fbtu2l*ctu
     &  + 4.d0*mfer2*h1fbtu*h1fbtu2r*ctu
     &  + 2.d0*mfer2*h1fbtu*h1fbtu3l*ctu
     &  + 2.d0*mfer2*h1fbtu*h1fbtu3r*ctu
     &  + 2.d0*mfer2*h1fbts*h1fbts1l*cts
     &  + 2.d0*mfer2*h1fbts*h1fbts1r*cts
     &  + 2.d0*mfer2*h1fbts*h1fbts2l*cts
     &  + 2.d0*mfer2*h1fbts*h1fbts2r*cts
     &  + 4.d0*mfer2*h1fbsu*h1fbsu2l*csu
     &  + 4.d0*mfer2*h1fbsu*h1fbsu2r*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu3l*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu3r*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu4l*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu4r*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu7l*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu7r*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu8l*csu
     &  + 2.d0*mfer2*h1fbsu*h1fbsu8r*csu
     &  + mfer2**2*h1fbuu1l*cuu
     &  + mfer2**2*h1fbuu1r*cuu
     &  + s*t*h1fbuu1l*cuu
     &  + s*t*h1fbuu1r*cuu
     &  + s*t*h1fbss*h1fbss1l*css
     &  + s*t*h1fbss*h1fbss1r*css
     &  - 2.d0*s*t*h1fbsu*h1fbsu1l*csu
     &  - 2.d0*s*t*h1fbsu*h1fbsu1r*csu
     &  - s*h1fbuu2l*cuu
     &  - s*h1fbuu2r*cuu
     &  - s*h1fbuu4l*cuu
     &  - s*h1fbuu4r*cuu
     &  - s*h1fbuu6l*cuu
     &  - s*h1fbuu6r*cuu
     &  - 2.d0*s*h1fbuu7l*cuu
     &  - 2.d0*s*h1fbuu7r*cuu
     &  - s*h1fbuu8l*cuu
     &  - s*h1fbuu8r*cuu
     &  + s*h1fbss*h1fbss2l*css
     &  + s*h1fbss*h1fbss2r*css
     &  + s*h1fbss*h1fbss4l*css
     &  + s*h1fbss*h1fbss4r*css
     &  + s*h1fbss*h1fbss6l*css
     &  + s*h1fbss*h1fbss6r*css
     &  + 2.d0*s*h1fbss*h1fbss7l*css
     &  + 2.d0*s*h1fbss*h1fbss7r*css
     &  + s*h1fbss*h1fbss8l*css
     &  + s*h1fbss*h1fbss8r*css
     &  - 2.d0*s*h1fbtu*h1fbtu2l*ctu
     &  - 2.d0*s*h1fbtu*h1fbtu2r*ctu
     &  - 2.d0*s*h1fbtu*h1fbtu3l*ctu
     &  - 2.d0*s*h1fbtu*h1fbtu3r*ctu
     &  + 2.d0*s*h1fbts*h1fbts2l*cts
     &  + 2.d0*s*h1fbts*h1fbts2r*cts
     &  + 2.d0*s*h1fbts*h1fbts3l*cts
     &  + 2.d0*s*h1fbts*h1fbts3r*cts
     &  - 2.d0*s*h1fbsu*h1fbsu2l*csu
     &  - 2.d0*s*h1fbsu*h1fbsu2r*csu
     &  + 2.d0*s*h1fbsu*h1fbsu4l*csu
     &  + 2.d0*s*h1fbsu*h1fbsu4r*csu
     &  + 2.d0*s*h1fbsu*h1fbsu6l*csu
     &  + 2.d0*s*h1fbsu*h1fbsu6r*csu
     &  - 2.d0*s*h1fbsu*h1fbsu8l*csu
     &  - 2.d0*s*h1fbsu*h1fbsu8r*csu
     &  + s**2*h1fbuu1l*cuu
     &  + s**2*h1fbuu1r*cuu
     &  + s**2*h1fbss*h1fbss1l*css
     &  + s**2*h1fbss*h1fbss1r*css
     &  - 2.d0*s**2*h1fbsu*h1fbsu1l*csu
     &  - 2.d0*s**2*h1fbsu*h1fbsu1r*csu
     &  - t*h1fbtt*h1fbtt1l*ctt
     &  - t*h1fbtt*h1fbtt1r*ctt
     &  - t*h1fbuu2l*cuu
     &  - t*h1fbuu2r*cuu
     &  - t*h1fbuu3l*cuu
     &  - t*h1fbuu3r*cuu
     &  - t*h1fbuu4l*cuu
     &  - t*h1fbuu4r*cuu
     &  - t*h1fbuu6l*cuu
     &  - t*h1fbuu6r*cuu
     &  - 2.d0*t*h1fbuu7l*cuu
     &  - 2.d0*t*h1fbuu7r*cuu
     &  - t*h1fbuu8l*cuu
     &  - t*h1fbuu8r*cuu
     &  - t*h1fbss*h1fbss3l*css
     &  - t*h1fbss*h1fbss3r*css
     &  - 2.d0*t*h1fbtu*h1fbtu1l*ctu
     &  - 2.d0*t*h1fbtu*h1fbtu1r*ctu
     &  - 2.d0*t*h1fbtu*h1fbtu2l*ctu
     &  - 2.d0*t*h1fbtu*h1fbtu2r*ctu
     &  - 2.d0*t*h1fbtu*h1fbtu3l*ctu
     &  - 2.d0*t*h1fbtu*h1fbtu3r*ctu
     &  - 2.d0*t*h1fbts*h1fbts1l*cts
     &  - 2.d0*t*h1fbts*h1fbts1r*cts
     &  - 2.d0*t*h1fbsu*h1fbsu2l*csu
     &  - 2.d0*t*h1fbsu*h1fbsu2r*csu
     &  - 2.d0*t*h1fbsu*h1fbsu3l*csu
     &  - 2.d0*t*h1fbsu*h1fbsu3r*csu
     &  - 2.d0*t*h1fbsu*h1fbsu8l*csu
     &  - 2.d0*t*h1fbsu*h1fbsu8r*csu
     &  + 2.d0*h1fbtt*h1fbtt2l*ctt
     &  + 2.d0*h1fbtt*h1fbtt2r*ctt
     &  + 2.d0*h1fbuu5l*cuu
     &  + 2.d0*h1fbuu5r*cuu
     &  + 2.d0*h1fbss*h1fbss5l*css
     &  + 2.d0*h1fbss*h1fbss5r*css
     &  + 4.d0*h1fbtu*h1fbtu4l*ctu
     &  + 4.d0*h1fbtu*h1fbtu4r*ctu
     &  + 4.d0*h1fbts*h1fbts4l*cts
     &  + 4.d0*h1fbts*h1fbts4r*cts
     &  + 4.d0*h1fbsu*h1fbsu5l*csu
     &  + 4.d0*h1fbsu*h1fbsu5r*csu

      if (dble(dsaschicasedas).lt.0.0d0) then     
        if(aszeroprint) then
        write(*,*) ' '
        write(*,*) 
     &    'DS: ERROR IN dsaschicased with negative cross section:'
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: kp1=',kp1,'  ',pname(kp1)
        write(*,*) 'DS: kp2=',kp2,'  ',pname(kp2)
        write(*,*) 'DS: kp3=',kp3,'  ',pname(kp3)
        write(*,*) 'DS: kp4=',kp4,'  ',pname(kp4)
        write(*,*) 'DS: dsaschicasedas = ',dsaschicasedas
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
        dsaschicasedas=dcmplx(0.0d0,0.0d0)
      endif

c...Convert to dW_eff/dcos(theta)
      par=dble(dsaschicasedas)*k34/(8.0d0*pi*gg1c*gg2c*s34*dsqrt(s))
      return
      end

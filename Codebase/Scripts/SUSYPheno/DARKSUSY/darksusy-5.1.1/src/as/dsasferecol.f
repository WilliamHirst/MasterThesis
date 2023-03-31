**************************************************************
*** SUBROUTINE dsasferecol                                 ***
*** computes dW_{ij}/dcostheta                             ***
*** sfermion(i) + sfermion(j)                              ***
*** -> fermion(k1) + fermion(k2)                           ***
*** version to be used if i or j are both squarks          ***
***                                                        *** 
*** input askin variables: p12,costheta                    ***
*** kpk1 and kpk2 are the fermion code                     ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-08-10                                         ***
*** modified: Piero Ullio to include a new labelling of    ***
*** states, 08-05-30                                       ***
**************************************************************

      SUBROUTINE dsasferecol(kpi,kpj,kpk1,kpk2,icase,result)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kpi,kpj,kpk1,kpk2,icase
      integer iii,jjj
      real*8 result
      real*8 ampl2
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c....check if the final state is kinematically allowed
      if((Svar-(mass(kpk1)+mass(kpk2))**2)/Svar.lt.thstep) then
        result=0.d0
        return
      endif      
***** type for kpk1 and kpk2:
      itype(3)=ivtype(kpk1)
      itype(4)=ivtype(kpk2)
***** equal / non equal particle final state: 
      if(kpk1.eq.kpk2) then
        s34=2.d0
      else  
        s34=1.d0
      endif  
***** fermion type in final state:
      q3=2
      q4=1
***** masses in final state:
      mass3=mass(kpk1)
      mass4=mass(kpk2)
***** define the kinematic variables
      call dsaskinset2
      call dsaskinset3
***** initialize ASx and ASy to zero
      do jjj=1,6
      do iii=1,4
      ASxplc(jjj,iii)=(0.d0,0.d0)
      ASxprc(jjj,iii)=(0.d0,0.d0)
      enddo
      ASylc(jjj)=(0.d0,0.d0)
      ASyrc(jjj)=(0.d0,0.d0)
      enddo
***** then fill them in:
*****
*************** case 1 ****************************************
*****
*****
***** up-type + up-type, or down-type + down-type *************
      if(icase.eq.1) then
***** f.s. particle consistent? *******************************
        if(itype(3).ne.itype(1).or.itype(4).ne.itype(1)) then
          write(*,*) 'DS: dsasferecol called with wrong particle'
          write(*,*) 'DS: in the final or initial state :'
          write(*,*) pname(kpi),pname(kpj)
     &               ,'and ',pname(kpk1),pname(kpk2)
          stop        
        endif  
***** neutralinos in t-channel **+*****************************
        do iii=1,4
          call dsasscscTCffcol(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
***** neutralinos in u-channel **+*****************************
        do iii=1,4
          call dsasscscUCffcol(kpi,kpj,kpk1,kpk2,kn(iii))
        enddo
***** gluino in t-channel *************************************
        call dsasscscTCffcol(kpi,kpj,kpk1,kpk2,kgluin)
***** gluino in u-channel *************************************
        call dsasscscUCffcol(kpi,kpj,kpk1,kpk2,kgluin)
*****
*************** case 2 ****************************************
*****
***** up-type + down-type *************************************

***** i. s.: particle - particle in same doublet **************
      elseif(icase.eq.2) then
***** f.s. particle consistent? *******************************
        if(itype(1).eq.itype(3).and.itype(2).eq.itype(4)) then
***** neutralinos in t-channel **+*****************************
          do iii=1,4
            call dsasscscTCffcol(kpi,kpj,kpk1,kpk2,kn(iii))
          enddo  
***** gluino in t-channel *************************************
          call dsasscscTCffcol(kpi,kpj,kpk1,kpk2,kgluin)
***** charginos in u-channel **********************************
          do iii=1,2
            call dsasscscUCffcol(kpi,kpj,kpk1,kpk2,kcha(iii))
          enddo
        elseif(itype(1).eq.itype(4).and.itype(2).eq.itype(3)) then 
***** neutralinos in u-channel **+*****************************
          do iii=1,4
            call dsasscscUCffcol(kpi,kpj,kpk1,kpk2,kn(iii))
          enddo  
***** gluino in u-channel *************************************
          call dsasscscUCffcol(kpi,kpj,kpk1,kpk2,kgluin)
***** charginos in t-channel **********************************
          do iii=1,2
            call dsasscscTCffcol(kpi,kpj,kpk1,kpk2,kcha(iii))
          enddo
        else  
          write(*,*) 'DS: dsasferecol called with wrong particle'
          write(*,*) 'DS: in the final or initial state :'
          write(*,*) pname(kpi),pname(kpj)
     &               ,'and ',pname(kpk1),pname(kpk2)
          stop        
        endif  
*****
*************** error *****************************************
*****
      else  
        write(*,*) 'DS: dsasferecol called with wrong particles'
        write(*,*) 'DS: in the initial state :'
        write(*,*) pname(kpi),pname(kpj)
        stop        
      endif  
*****
*****
***** compute the amplitude squared 
      call dsasffcol(ampl2)
***** computes dW_{ij}/dcostheta
      result=ampl2*k34/(8.d0*pi*gg1*gg2*s34*dsqrt(Svar))
      return
      end

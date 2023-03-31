! This file is part of HiggsBounds
!  -KW
!************************************************************      
program HiggsBounds
!************************************************************
! HiggsBounds uses theoretical cross section and branching ratios
! from theories containing an arbitrary number of neutral
! Higgs particles and compares to the experimental limits 
! from LEP, Tevatron and LHC
!************************************************************
 use usefulbits, only : theo,res,deallocate_usefulbits,debug, &
                   &    inputmethod,HiggsBounds_info,ndat, &
                   &    file_id_debug1,file_id_debug2
 use input, only : setup_input,do_input
 use S95tables, only : setup_S95tables,deallocate_S95tables
 use theory_BRfunctions, only : setup_BRSM,deallocate_BRSM
 use channels, only : setup_channels,check_channels
 use output, only : setup_output,do_output
 use theo_manip, only : complete_theo

#if defined(NAGf90Fortran)
 use F90_UNIX_IO, only : flush
#endif
       
 implicit none 
 !-----------------------------------internal      
 integer :: ii,jj,kk
 logical :: messages
 !-------------------------------------------               

!#define DEBUGGING
#ifdef DEBUGGING
  debug=.True.
#else
  debug=.False.
#endif

#ifdef WEBVERSION
 inputmethod='website'  
#else
 inputmethod='datfile'    !(inputmethod='subrout' is also possible, but do not set that here)
#endif

 messages=debug.or.(inputmethod=='datfile')

 if(inputmethod.eq.'datfile')call HiggsBounds_info

 if(debug)then  
  open(file_id_debug2,file='debug_predratio.txt')
  open(file_id_debug1,file='debug_channels.txt')
 endif
                                    
 if(messages)write(*,*)'doing some preliminary tasks...'      ; call flush(6)
 call setup_input

 if(messages)write(*,*)'reading in Standard Model tables...'   ; call flush(6) 
 call setup_BRSM 

 if(messages)write(*,*)'reading in S95tables...'               ; call flush(6)
 call setup_S95tables      
                        
 if(messages)write(*,*)'sorting out processes to be checked...'; call flush(6)
 call setup_channels      

 if(messages)write(*,*)'preparing output arrays...'            ; call flush(6)
 call setup_output
      
 if(messages)write(*,*)'getting theoretical input...'          ; call flush(6)
 call do_input
       
 if(messages)write(*,*)'manipulating input...'                 ; call flush(6)
 call complete_theo
            
 if(messages)write(*,*)'compare each data point to the experimental bounds...' 
                                                                 call flush(6)            
      
!$OMP PARALLEL PRIVATE(ii,jj,kk)
!$OMP DO
 do jj=1,ndat
  !do kk=1,1000000 ! for testing
  ! ii=1
  !enddo
  !if(mod(jj,1000).eq.0) write(*,*)'jj=',jj; call flush(6)

  call check_channels(theo(jj),res(jj))       
 enddo
!$OMP END DO 
!$OMP END PARALLEL
      
 if(debug)then
  close(file_id_debug2)
  close(file_id_debug1)
 endif

 if(messages)write(*,*)'beginning output...'                   ; call flush(6)
 call do_output

 if(messages)write(*,*)'finishing off...'                      ; call flush(6)
 call deallocate_S95tables
 call deallocate_BRSM
 call deallocate_usefulbits
      
 if(messages)write(*,*)'finished'                              ; call flush(6)             
end program HiggsBounds


      

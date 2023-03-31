
 !************************************************************
 program create_clsb_binary_files
 !************************************************************
 !make libHB
 !gfortran create_clsb_binary_files.f90 -o create_clsb_binary_files.exe -L. -lHB
 use S95tables, only : setup_S95tables,S95_t2
 use S95tables_type2
 use S95tables_type3
 use theory_BRfunctions, only : setup_BRSM

 implicit none            
 integer :: clsbtable_selection_id
 integer :: z

 allocate(clsb_t3(ntable3))
 call setup_BRSM !needed for setup_S95tables
 call setup_S95tables
 call initializetables_type3_blank(clsb_t3)
 call initializetables3(clsb_t3)
 call fillt3needs_M2_gt_2M1(clsb_t3,S95_t2)

 do z=lbound(clsb_t3,dim=1),ubound(clsb_t3,dim=1)
  clsbtable_selection_id = clsb_t3(z)%id
  call readclsbfiles(clsbtable_selection_id)
 enddo

 call writeclsbfiles_binary

 end program create_clsb_binary_files

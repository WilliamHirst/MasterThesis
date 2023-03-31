**********************************************************************
*** subroutine dsgalproptable
*** input: header string to pass to galprop
*** 
*** author: e.a. baltz 2/21/2006
**********************************************************************

      subroutine dsgalprop_run_galprop(header)
      implicit none
      character*80 header
      
      include 'dsgalpropcom.h'

      call galprop(header)
      end

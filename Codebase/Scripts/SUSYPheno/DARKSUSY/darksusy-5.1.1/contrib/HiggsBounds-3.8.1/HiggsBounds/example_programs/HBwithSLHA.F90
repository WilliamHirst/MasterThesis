
!******************************************************
program HBwithSLHA
!
! Run with
!   ./HBwithSLHA npoints <stem>
! where npoints is the number of parameter points you would like to
! look at and each parameter point has a corresponding SLHA file
! e.g. the corresponding SLHA for the 5th point should be found at
!   <stem>.5
!
! Output
! The block HiggsBoundsResults will be added to each SLHA file.
! In addition, the HiggsBounds results will be collected together in 
! the file
!   <stem>-fromHB
!
!******************************************************     
#ifdef NAGf90Fortran
 use F90_UNIX_ENV, only : iargc,getarg
#endif
  implicit none
  integer :: nH,nHplus,HBresult,chan,ncombined
  double precision :: obsratio
  integer :: i,npoints
  integer,parameter :: fileid=78
  character(len=8) :: istring
  character(len=300) :: inputfilename,outputfilename
  character(len=300) :: stem
  character(LEN=300) :: temp
  integer :: number_args
#ifndef NAGf90Fortran
  integer :: iargc
#endif

  nH=3
  nHplus=1

  number_args = IARGC() 
  
  if( number_args .ne. 2)then
   stop "Incorrect number of arguments given to HBwithSLHA"
  endif

  ! Read arguments into text strings.
  i=1
  temp=""
  call GETARG(i,temp)
  read(temp,*) npoints

  i=i+1  
  temp=""
  call GETARG(i,temp)
  stem = ""
  stem = trim(temp)

  call initialize_HiggsBounds(nH,nHplus,'LandH')

  outputfilename=trim(adjustl(stem))//'-fromHB'

  open(fileid, file=trim(outputfilename))

  do i=1,npoints

   if(i.gt.99999999)stop'need to increase the size of istring in HBwithSLHA'
   write(istring,'(I8)')i

   inputfilename=trim(adjustl(stem))//'.'//trim(adjustl(istring))

   call HiggsBounds_input_SLHA(inputfilename)

   call run_HiggsBounds( HBresult, chan, obsratio, ncombined )

   !This will add the block HiggsBoundsResults to the SLHA file
   call HiggsBounds_SLHA_output

   !This will collect all the HiggsBounds results together into one file
   write(fileid,*)i,HBresult,chan,obsratio,ncombined

  enddo
  
  close(fileid)

  call finish_HiggsBounds

end program HBwithSLHA

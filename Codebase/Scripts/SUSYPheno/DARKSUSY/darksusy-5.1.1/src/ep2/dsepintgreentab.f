      real*8 function dsepintgreentab(DeltaV)
****************************************************************
***                                                          ***
*** tabulated version of dsepintgreen, see header of the     ***
*** latter for details. on first call the table is either    ***
*** computed or read from file provided by the user.         ***
*** This is handled by the flag dhow in common block /ephow/ ***
***       dhow = 2 - table is tabulated on first             ***
***                  call, and then interpolated             ***
***              3 - as 2, but also write the table to disk  ***
***                  at the first call                       ***
***              4 - read table from disk on first call, and ***
***                  use the subsequent calls. If the file   ***
***                  does not exist, it will be created      ***
***                  (as in 3). (use as default)             ***
***                                                          ***
*** To load another table reset variable epcl_load           ***
*** in dsepcom.h to .true. (done automatically with call     ***
*** to dsepset.f)                                            ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-02-03                                         ***
*** Modified: Joakim Edjso, modified table loading           ***
****************************************************************
      implicit none
      include 'dsepcom.h'
      include 'dshmcom.h' ! for haloid
      include 'dsdirver.h'
      integer ii
      real*8 deltav,result,dsepintgreen
      real*8 xdv(1000),ydv(1000),ydv2(1000),realndv
      real*8 xmax,xmin,xint,ymax
      character*200 epfile,scr
      integer i,k,kk,npoints
      common/bigItab2com/xdv,ydv,ydv2,realndv

      integer dhow
      common/ephow/dhow


 5    if (epcl_load) then    ! make a new table or load a new table from disk

c...Generate file name
        epfile=dsinstall
        call dscharadd(epfile,'share/DarkSUSY/epgretab-')
        call dscharadd(epfile,epc)
        call dscharadd(epfile,'-')
        call dscharadd(epfile,haloid)
        call dscharadd(epfile,'.dat')

        if (dhow.eq.2.or.dhow.eq.3) then  ! generate a new table
          write(*,*) 'starting dsepintgreen tabulation'  
          ymax=0.d0
          xmin=0.1d0
          xmax=10000.d0
 110      continue
          ii=0
          do i=0,100
            deltav=dexp(dlog(xmin)+(dlog(xmax)-dlog(xmin))/100.d0*i)
            result=dsepintgreen(DeltaV)
            if(i.eq.0.and.dabs(1.d0-result).gt.1.d-2) then
              xmin=xmin/10.d0
              goto 110
            endif
            ii=ii+1
            xdv(ii)=dlog(deltav)
            ydv(ii)=result
            if(ydv(ii).gt.ymax) ymax=ydv(ii)
            write(*,*) 'ii,deltav, result = ',ii,deltav,result
            if(ydv(ii).lt.1.d-8) goto 111
          enddo
 111      npoints=ii
          write(*,*) 'Number of points: ',npoints
 130      continue
          do k=1,npoints-1
            if(dabs(ydv(k)-ydv(k+1)).gt.0.1d0*max(ydv(k),ydv(k+1)).and.
     &        max(ydv(k),ydv(k+1)).gt.1.d-4*ymax) then
              xint=(xdv(k)+xdv(k+1))/2.d0
              deltav=dexp(xint)
              result=dsepintgreen(deltav)
              do kk=npoints,k+1,-1
                xdv(kk+1)=xdv(kk)
                ydv(kk+1)=ydv(kk)
              enddo
              xdv(k+1)=xint
              ydv(k+1)=result
              write(*,*) 'adding deltav, result = ',deltav,result
              npoints=npoints+1  
              if(npoints.le.1000) then
                goto 130 
              else
                write(*,*) 'in dsepintgreentab exceeded the maximum dim'
                write(*,*) 'allowed for vectors in the compspcom block'
                write(*,*) 'which is set equal to 1000'
                stop
              endif
            endif  
          enddo
          write(*,*) 'dsepintgreen tabulation is over'  
          if (dhow.eq.3) then ! write table to disk
            write(*,*) 'Writing ep table to file ',epfile
            open(unit=13,file=epfile,status='unknown',form='formatted')
            write(13,1001) dsversion
            write(13,1002) 
 1001       format('# Made with DarkSUSY version ',A)
 1002       format('#','..log(deltav).',1x,'...table(i)...')
            do i=1,npoints
              write(13,1000)  xdv(i),ydv(i) 
            enddo  
            close(13)
            write(*,*) 'Done.'
          endif

        elseif (dhow.eq.4) then  ! read table from file
          write(*,*) 'Reading ep table from file ',epfile
          ii=0
          open(unit=13,file=epfile,status='old',form='formatted',
     &       err=200)
          read(13,'(a)') scr  ! read header line
          read(13,'(a)') scr  ! read header line
 10       ii=ii+1
          if(ii.gt.1000) then
            write(*,*) 'in dsepintgreentab epfile has more than 1000'
            write(*,*) 'entries, increase size of var in bigItab2com'
            write(*,*) 'program stopped'
            stop
          endif  
          read(13,1000,end=11)  xdv(ii),ydv(ii) 
          goto 10
 1000     format(60(1x,e14.8))
 11       close(13)
          npoints=ii-1  ! for later spline setup
          write(*,*) 'Done.'
        else
          write(*,*) 'ERROR in dsepspecm: how = ',dhow,
     &     'requested, but this is out of range. Stopping.'
          stop
        endif

c...Now set up splines
        realndv=dble(npoints)
        call dsspline(xdv,ydv,int(realndv),1.d31,1.d31,ydv2)

        epcl_load=.false.  ! do not reload on next call

      endif

c...Now do the actual table lookup
      if(deltav.lt.dexp(xdv(1))) then
        if(dabs(1.d0-ydv(1)).lt.1.d-2) then 
          dsepintgreentab=1.d0
        else  
          dsepintgreentab=ydv(1)
        endif
      elseif(deltav.gt.dexp(xdv(int(realndv)))) then
        dsepintgreentab=0.d0
      else   
        call dssplint(xdv,ydv,ydv2,int(realndv),dlog(deltav),result)
        dsepintgreentab=result
      endif

      return

c----------------------------------------------------------------------

c...Now take care of missing file
c...We get here if dhow=4, but file does not exist
 200  close(13)
      write(*,*) 'The requested ep table file ',epfile
      write(*,*) 'does not exist. Will recreate if for you.'
      dhow=3
      goto 5


      end


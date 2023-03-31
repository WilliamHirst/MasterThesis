      real*8 function dsdbtd15(tp,howinp)

**********************************************************************
*** function dsdbtd15 is the containment time in 10^15 sec
***   input:
***     tp - kinetic energy per nucleon in gev
***     how - 1 calculate t_diff only for requested momentum
***           2 tabulate t_diff for first call and use table for
***             subsequent calls
***           3 as 2, but also write the table to disk as dbtd.dat
***           4 read table from disk on first call, and use that for
***             subsequent calls
***   output:
***     t_diff in units of 10^15 sec
*** calls dsdbtd15x for the actual calculation.
*** author: joakim edsjo (edsjo@physto.se)
*** uses piero ullios propagation routines.
*** date: dec 16, 1998
*** modified: 98-07-13 paolo gondolo
**********************************************************************

      implicit none
      include 'dshmcom.h'
      include 'dspbcom.h'
      include 'dsdirver.h'
      real*8 tp
      integer how,howinp
      logical tabuldb(2)
      integer npt,i
      real*8 tmin,tmax,ltmin,ltmax,dlt,t,ppl,dsdbtd15x,ltlow,lthigh,tt,
     &  dsdbtd15point
      parameter(npt=250)
      parameter(tmin=0.01d0,
     &          tmax=1000.0d0)
      real*8 carrdb(0:npt,2)
      character*128 scr,dbfile

      data tabuldb /.false.,.false./
      save tabuldb,carrdb


      how=howinp

c...Define filename to save/read data from if how=3 or 4
      dbfile=dsinstall
      call dscharadd(dbfile,'share/DarkSUSY/dbtd-')
      call dscharadd(dbfile,pbc)
      call dscharadd(dbfile,'-')
      call dscharadd(dbfile,haloid)
      call dscharadd(dbfile,'.dat')

c...tabulate if needed
 10   if ((how.eq.2.or.how.eq.3).and.(.not.tabuldb(hclumpy))) then
        write(*,*) 'dsdbtd15: tabulation started...'
        ltmin=log(tmin)
        ltmax=log(tmax)
        dlt=ltmax-ltmin
        if (how.eq.3) then  ! write to disk
          open(unit=27,file=dbfile,status='unknown',
     &      form='formatted')
          write(27,1001)
 1001     format('#','......tp......',1x,'...t_diff....')
        endif
        do i=0,npt
          t=exp(ltmin+dlt*dble(i)/dble(npt))
          if (pbrcy.gt.0.0d0) then
            carrdb(i,hclumpy)=dsdbtd15point(pbrho2int,t)+dsdbtd15x(t)
          else
            carrdb(i,hclumpy)=dsdbtd15x(t)
          endif
          if (how.eq.3) then   ! write to disk
            write(27,1002) t,carrdb(i,hclumpy)
 1002       format(2(1x,e14.8))
          endif
        enddo
        if (how.eq.3) close(27)
        tabuldb(hclumpy)=.true.
        write(*,*) '  done.'
      endif

c...read table from disk if requested
      if ((how.eq.4).and.(.not.tabuldb(hclumpy))) then
        ltmin=log(tmin)
        ltmax=log(tmax)
        dlt=ltmax-ltmin
        write(*,*) 'dsdbtd15: reading dbtd from file ',dbfile
        open(unit=27,file=dbfile,status='unknown',
     &    form='formatted')
        read(27,'(a)',end=100) scr
        do i=0,npt
          t=exp(ltmin+dlt*dble(i)/dble(npt))
          read(27,1002,end=100) tt,carrdb(i,hclumpy)
c...          check if file is ok
          if (abs((tt-t)/t).gt.0.01d0) then
            write(*,*) 'warning in dsdbtd15: the file dbtd.dat'
            write(*,*) '  does not have correct entries: t mismatch'
            write(*,*) 'will create a new data file instead.'
            close(27)
            how=3
            goto 10
          endif
        enddo
        close(27)
        tabuldb(hclumpy)=.true.
        write(*,*) '  done.'
      endif

c...use tabulated if requested
      if (how.ge.2.and.how.le.4) then
        if (tp.lt.tmin) write(*,*)
     &    'error in dsdbtd15: tp=',tp,' is less than tmin=',tmin
        if (tp.gt.tmax) write(*,*)
     &    'error in dsdbtd15: tp=',tp,' is greater than tmax=',tmax
        ltmin=log(tmin)
        ltmax=log(tmax)
        dlt=ltmax-ltmin
        i=int((log(tp)-ltmin)*npt/dlt)
        ltlow=ltmin+dlt*dble(i)/dble(npt)
        lthigh=ltmin+dlt*dble(i+1)/dble(npt)
        ppl=(log(tp)-ltlow)/(lthigh-ltlow)
        dsdbtd15=(1.0d0-ppl)*carrdb(min(i,npt),hclumpy) +
     &    ppl*carrdb(min(i+1,npt),hclumpy)
      endif

c...call dsdbtd15 directly
      if (how.eq.1) then
        if (pbrcy.gt.0.0d0) then
          dsdbtd15=dsdbtd15point(pbrho2int,tp)+dsdbtd15x(tp)
        else
          dsdbtd15=dsdbtd15x(tp)
        endif
      endif

      return

c...bad file read
 100  write(*,*) 'warning in dsdbtd15: the file',dbfile
      write(*,*) '  does either not exist or is incomplete.'
      write(*,*) 'will create a new data file instead.'
      close(27)
      how=3
      goto 10


      end




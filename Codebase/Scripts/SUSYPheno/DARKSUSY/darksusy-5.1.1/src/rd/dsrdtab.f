      subroutine dsrdtab(wrate,xmin)
c_______________________________________________________________________
c  tabulate the invariant annihilation rate as a function of p.
c  input:
c    wrate - invariant annihilation rate (real*8, external)
c    xmin - minimum mass/temperature needed (real*8)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses dsrdnormlz, dsrdlny, dsrdspline.
c  authors: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994 and
c           joakim edsjo (edsjo@physto.se) 06-march-98
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wrate,xmin,p
      external wrate
      integer i,j,k,nadd
      real*8 pmin,dp,dp1,dp2,x1,x2,dy1,dy2,y1,y2,
     & ymin,ymax,deltap,deltay,dsrdlny,pres,dpth,pthlow,
     & wspline,wlin,wrmax,wrmin,pwmax,pwmin,p2,
     & dpmin,dpthmin,pp1,pp2,w1,w2,w3,padd(nrmax),dsrdthclose,
     & pthr,dsrddpmin
      real*8 c(nrmax)

cc      write (*,*) 'Entering dsrdtab'
cc      write (*,*) 'nth=',nth

c---------------------- minimum acceptable cosine
c      cosmin=cos(3.14d0/180.d0*5.0d0)
c     cosmin=cos(3.14d0/180.d0*5.0d0) ! je change to increase efficiency
c... cosmin set in block data

c--------------------- minimum acceptable spacing in p
c      dpmin=1.d-4
c      dpmin=0.005d0 ! set in block data
      dpmin=dpminr*mco(1) ! smallest allowed distance between points
      dpthmin=dpthr*mco(1)
      dpth=dpthmin*0.5d0 ! set scale for sampling around thresholds

c      if (dpthr.ge.dpminr) then
c        write(*,*) 'warning in dsrdtab: dpthr greater than dpminr'
c        write(*,*)
c     &    '  this can cause unexpected behaviour at thresholds.'
c      endif

c========================= first tabulation and rescaling to unit square
cc      write (*,*) 'dsrdtab: start first tabulation'

      pmin=0.d0
      if (wrate(0.d0).eq.0.d0) pmin=1.d-10*mco(1)

      pmax=mco(1)*umax*sqrt((1.0d0+0.25d0*umax**2/xmin)/xmin)
c...Avoid too high momenta, these are not needed except at
c...very small x~2 which is not important anyway. Instead freeze
c...w(p) at w(pmax) when p>pmax and lower pmax.
      pmax=min(10.0d0*mco(1),pmax)


c      write(*,*) 'model = ',rdtag
c      write(*,*) '  pmax = ',pmax
      pthlow=pmax
c      deltap=pmax-pmin
c      dp=1.0d0/(dble(nmin)-1.0d0)
      deltap=mco(1)  ! sets scale for angle routine

      indx(1)=1
      pp(1)=pmin/deltap
      yy(1)=dsrdlny(pmin,wrate)

c...first sample low p region dense, je change 97-03-05
cc      write (*,*) 'dsrdtab: start low p region'
      dp=(pdivr*mco(1)-pmin)/deltap/(dble(nlow)-1.0d0)  ! je test
      do i=2,nlow
        indx(i)=i
        pp(i)=pp(1)+dble(i-1)*dp
        yy(i)=dsrdlny(pp(i)*deltap,wrate)
      enddo

c...then sample high p region sparse
cc      write (*,*) 'dsrdtab: start high p region'
      dp=(pmax-pdivr*mco(1))/deltap/(dble(nhigh))  ! je test
      do i=nlow+1,nlow+nhigh
        indx(i)=i
        pp(i)=pp(nlow)+dble(i-nlow)*dp
        yy(i)=dsrdlny(pp(i)*deltap,wrate)
      enddo
      nr=nlow+nhigh


c------------------- add one point close to p=0 to insure correct point
c----------------------------- insertion at low p, je addition 96-04-09
c
c      call dsrdaddpt(wrate,pmin+dpmin*deltap,deltap)

c------------------------ add points at resonances je addition 96-03-28

cc      write (*,*) 'dsrdtab: start resonances'
      i=0
 6    i=i+1
      if (i.le.nres) then
 7      pres=sqrt(rgev(i)**2/4.0d0-mco(1)**2)
        if (pres.lt.pmax) then
c...first check if important
          j=-npres
          p2=(rgev(i)+dpres*dble(j)*rwid(i))**2/4.0d0-mco(1)**2
          if (p2.gt.0.0d0) then
            pres=sqrt(p2)
            call dsrdaddpt(wrate,pres,deltap)
            if (rderr.ne.0) goto 999 ! return
            w1=exp(yy(nr))
          else
            w1=0.0d0
          endif

          j=npres
          p2=(rgev(i)+dpres*dble(j)*rwid(i))**2/4.0d0-mco(1)**2
          if (p2.gt.0.0d0) then
            pres=sqrt(p2)
            call dsrdaddpt(wrate,pres,deltap)
            if (rderr.ne.0) goto 999 ! return
            w3=exp(yy(nr))
          else
            w3=0.0d0
          endif

          j=0
          p2=(rgev(i)+dpres*dble(j)*rwid(i))**2/4.0d0-mco(1)**2
          if (p2.gt.0.0d0) then
            pres=sqrt(p2)
            call dsrdaddpt(wrate,pres,deltap)
            if (rderr.ne.0) goto 999 ! return
            w2=exp(yy(nr))
          else
            w2=1.0d0
          endif

          w1=0.5d0*(w1+w3)
          if (abs((w2-w1)/w1).lt.wdiffr) then ! unimportant resonance
            do k=i,nres-1
              rgev(k)=rgev(k+1)
              rwid(k)=rwid(k+1)
            enddo
            nres=nres-1
            if (i.le.nres) goto 7
          else
            do j=-npres+1,npres-1
              if (j.ne.0) then
                p2=(rgev(i)+dpres*dble(j)*rwid(i))**2/4.0d0-mco(1)**2
                if (p2.gt.0.0d0) then
                  pres=sqrt(p2)
                  call dsrdaddpt(wrate,pres,deltap)
                  if (rderr.ne.0) goto 999 ! return
                endif
              endif
            enddo
          endif
        endif
        goto 6
      endif

c------------------------ add points at thresholds je addition 96-03-28
cc      write (*,*) 'dsrdtab: start thresholds'
      i=0
 16   i=i+1
      if (i.le.nth) then
 17     if (pth(i).lt.pthlow) pthlow=pth(i)
        if (pth(i).lt.pmax) then
          call dsrdaddpt(wrate,pth(i)-dpth,deltap) ! je test
          if (rderr.ne.0) goto 999 ! return 
          w1=exp(yy(nr))
          do j=nthup,1,-1
            call dsrdaddpt(wrate,pth(i)+dble(j)**2*dpth,deltap)
            if (rderr.ne.0) goto 999 ! return 
            w2=exp(yy(nr))
            if (j.eq.nthup.and.(w2-w1)/w1.lt.wdifft) then ! unimportant thr.
               incth(i)=0
c..   ! still included in spline split, but not in integration split
c...delete threshold
c              do k=i,nth-1
c                pth(k)=pth(k+1)
c              enddo
c              nth=nth-1
c              if (i.le.nth) goto 17
            endif
          enddo
        endif
        goto 16
      endif

c------------------------------------- add points below lowest threshold
      if (nth.ge.1) then
        do i=1,3
          call dsrdaddpt(wrate,pthlow*dble(i)/4.0d0,deltap)
          if (rderr.ne.0) goto 999 ! return 
        enddo
      endif

c======================================================= cos(theta) test

      if (cthtest.eq.1) then ! perform cos(theta) test
cc         write (*,*) 'dsrdtab: start cos(theta) test'

c----------------------------------------------- determine ymin and ymax
 5      ymin=yy(1)
        ymax=yy(1)
        do i=2,nr
          if (yy(i).lt.ymin) ymin=yy(i)
          if (yy(i).gt.ymax) ymax=yy(i)
        enddo
c-----------------------------------------------------------------------
        deltay=ymax-ymin
        if (deltay.eq.0.0d0) deltay=1.0d0
        do i=1,nr
          yy(i)=yy(i)/deltay
        enddo
c--------------------- tabulation of cosines between successive segments
        c(1)=1.0d0
        do i=2,nr-1
          dp1=(pp(indx(i))-pp(indx(i-1)))
          dy1=(yy(indx(i))-yy(indx(i-1)))
          dp2=(pp(indx(i+1))-pp(indx(i)))
          dy2=(yy(indx(i+1))-yy(indx(i)))
          call dsrdnormlz(dp1,dy1,x1,y1)
          call dsrdnormlz(dp2,dy2,x2,y2)
          c(indx(i))=x1*x2+y1*y2
        enddo
        c(nr)=1.0d0

c============================================================== scanning

        i=2
 10     continue
c----------------------------- if needed, insert points at the two sides
          if ((c(indx(i)).lt.cosmin).and.
     &     (pp(indx(i+1))-pp(indx(i-1))).gt.dpmin) then
            if(nr+2.gt.nrmax) then
              write (rdluerr,*)
     &          'error in dsrdtab: array capacity exceeded'
              write (rdluerr,*) '  for model ',rdtag
              write(rdluerr,*) 'Omega calculation stopping.'
              nr=nr-2
              rderr=ibset(rderr,0)
              goto 999 ! return
            endif
c------------------------------------- point and cosine on the left side
            nr=nr+1
            pp(nr)=0.5d0*(pp(indx(i))+pp(indx(i-1)))
            yy(nr)=dsrdlny(pp(nr)*deltap,wrate)/deltay
            dp1=(pp(nr)-pp(indx(i-1)))
            dp2=(pp(indx(i))-pp(nr))
            dy1=(yy(nr)-yy(indx(i-1)))
            dy2=(yy(indx(i))-yy(nr))
            call dsrdnormlz(dp1,dy1,x1,y1)
            call dsrdnormlz(dp2,dy2,x2,y2)
            c(nr)=x1*x2+y1*y2
c------------------------------------ point and cosine on the right side
            nr=nr+1
            pp(nr)=0.5d0*(pp(indx(i+1))+pp(indx(i)))
            dp1=(pp(nr)-pp(indx(i)))
            dp2=(pp(indx(i+1))-pp(nr))
            yy(nr)=dsrdlny(pp(nr)*deltap,wrate)/deltay
            dy1=(yy(nr)-yy(indx(i)))
            dy2=(yy(indx(i+1))-yy(nr))
            call dsrdnormlz(dp1,dy1,x1,y1)
            call dsrdnormlz(dp2,dy2,x2,y2)
            c(nr)=x1*x2+y1*y2
c--------------------------------------- new cosine at the central point
            dp1=(pp(nr)-pp(indx(i)))
            dp2=(pp(indx(i))-pp(nr-1))
            dy1=(yy(nr)-yy(indx(i)))
            dy2=(yy(indx(i))-yy(nr-1))
            call dsrdnormlz(dp1,dy1,x1,y1)
            call dsrdnormlz(dp2,dy2,x2,y2)
            c(indx(i))=x1*x2+y1*y2
c---------------------------- new cosine to the left of left added point
c...je addition 96-03-29
            if (i.gt.2) then
              dp1=(pp(indx(i-1))-pp(indx(i-2)))
              dp2=(pp(nr-1)-pp(indx(i-1)))
              dy1=(yy(indx(i-1))-yy(indx(i-2)))
              dy2=(yy(nr-1)-yy(indx(i-1)))
              call dsrdnormlz(dp1,dy1,x1,y1)
              call dsrdnormlz(dp2,dy2,x2,y2)
              c(indx(i-1))=x1*x2+y1*y2
            endif
c---------------------------- new cosine to the right of right added point
c...je addition 96-03-29
            if (i.lt.(nr-3)) then
              dp1=(pp(indx(i+1))-pp(nr))
              dp2=(pp(indx(i+2))-pp(indx(i+1)))
              dy1=(yy(indx(i+1))-yy(nr))
              dy2=(yy(indx(i+2))-yy(indx(i+1)))
              call dsrdnormlz(dp1,dy1,x1,y1)
              call dsrdnormlz(dp2,dy2,x2,y2)
              c(indx(i+1))=x1*x2+y1*y2
            endif
c--------------------------------------------------------- shift indexes
            do j=nr,i+3,-1
              indx(j)=indx(j-2)
            enddo
            indx(i+2)=nr
            indx(i+1)=indx(i)
            indx(i)=nr-1
            i=max(i-1,2)
c------------------------------------------- else, move on to next point
          else
            i=i+1
            if (i.ge.nr) goto 20
          endif
          if (rderr.ne.0) goto 999 ! return
        goto 10
c---------------------------------------- scale back to user coordinates
 20     do i=1,nr
          pp(i)=pp(i)*deltap
          yy(i)=yy(i)*deltay
        enddo

      else ! no cos(theta) test
        do i=1,nr
          pp(i)=pp(i)*deltap
        enddo
      endif

c================================================= spline setup and test

cc      write (*,*) 'dsrdtab: start spline setup'
 25   call dsrdspline

c-------------------------------------- check cubic spline interpolation
c...addition by je 96-04-09

cc      write (*,*) 'dsrdtab: start checking spline interpolation'
      if (spltest.eq.1) then
        wrmax=1.0d0
        wrmin=1.0d0
        nadd=0

cc        write (*,*) 'dsrdtab: ..... first loop'
        do i=1,nr-1
          pp1=pp(indx(i))
          pp2=pp(indx(i+1))
          p=0.5d0*(pp1+pp2)
          if ((pp2-pp1).gt.dsrddpmin(p,dpmin)) then
            call dsrdwintpch(p,wspline,wlin)
            if (wspline/wlin.gt.(1.0d0+waccd).or.
     &        wspline/wlin.lt.(1.0d0-waccd)) then
              nadd=nadd+1
              padd(nadd)=p
            endif

c          if (wspline/wlin.gt.wrmax) then
c            wrmax=wspline/wlin
c            pwmax=p
c          endif
c          if (wspline/wlin.lt.wrmin) then
c            wrmin=wspline/wlin
c            pwmin=p
c          endif
          endif
        enddo


        if (nadd.ge.1) then
cc          write(*,*) 'now adding ',nadd,' points to nr = ',nr
          do i=1,nr
            pp(i)=pp(i)/deltap
          enddo
          do i=1,nadd
            call dsrdaddpt(wrate,padd(i),deltap)
            if (rderr.ne.0) goto 999 ! return 
          enddo
          do i=1,nr
            pp(i)=pp(i)*deltap
          enddo
          if (rderr.ne.0) goto 999 ! return
          goto 25
        endif

c      if (wrmax.gt.(1.0d0+waccd).or.
c     &  wrmin.lt.(1.0d0-waccd)) then
c        do i=1,nr
c          pp(i)=pp(i)/deltap
c        enddo
c        if (wrmax.gt.(1.0d0+waccd))
c     &    call dsrdaddpt(wrate,pwmax,deltap)
c        if (wrmin.lt.(1.0d0-waccd))
c     &    call dsrdaddpt(wrate,pwmin,deltap)
c        if (rderr.gt.0) goto 999 ! return
c        do i=1,nr
c          pp(i)=pp(i)*deltap
c        enddo
c        goto 25   ! if goto 5, take the three lines above away
c      endif

c...check how bad we are doing when delta_p < dpmin

        wrmax=1.0d0
        wrmin=1.0d0

cc        write (*,*) 'dsrdtab: ..... second loop'        
        do i=1,nr-1
          pp1=pp(indx(i))
          pp2=pp(indx(i+1))
          p=0.5d0*(pp1+pp2)
c...check if so close to a thr. that it doesn't matter.
          pthr=dsrdthclose(p)
          if (abs(p-pthr).gt.2.0d0*dpth) then ! will be used
            call dsrdwintpch(p,wspline,wlin)
            if (wspline/wlin.gt.wrmax) then
              wrmax=wspline/wlin
              pwmax=p
            endif
            if (wspline/wlin.lt.wrmin) then
              wrmin=wspline/wlin
              pwmin=p
            endif
          endif
        enddo

cc        write (*,*) 'dsrdtab: ..... set warnings'
        if (wrmax.gt.(1.0d0+5.0d0*waccd)) then
          write(*,*) 'warning in dsrdtab: wsplin/wlin = ',wrmax
          write(*,*) '  for p=',pwmax,' and model ',rdtag
          write(*,*) '  points not inserted since delta_p<dpmin = ',
     &      dpmin
        endif
        if (wrmax.gt.(1.0d0+5.0d0*waccd).and.
     &    wrmax.le.(1.0d0+10.0d0*waccd)) rdwar=ibset(rdwar,0)
        if (wrmax.gt.(1.0d0+10.0d0*waccd).and.
     &    wrmax.le.(1.0d0+15.0d0*waccd)) rdwar=ibset(rdwar,1)
        if (wrmax.gt.(1.0d0+15.0d0*waccd)) rdwar=ibset(rdwar,2)

        if (wrmin.lt.(1.0d0-5.0d0*waccd)) then
          write(*,*) 'warning in dsrdtab: wsplin/wlin = ',wrmin
          write(*,*) '  for p=',pwmin,' and model ',rdtag
          write(*,*) '  points not inserted since delta_p<dpmin = ',
     &      dpmin
        endif
        if (wrmin.lt.(1.0d0-5.0d0*waccd).and.
     &    wrmin.ge.(1.0d0-10.0d0*waccd)) rdwar=ibset(rdwar,0)
        if (wrmin.lt.(1.0d0-10.0d0*waccd).and.
     &    wrmin.ge.(1.0d0-15.0d0*waccd)) rdwar=ibset(rdwar,1)
        if (wrmin.lt.(1.0d0-15.0d0*waccd)) rdwar=ibset(rdwar,2)

      endif ! end of spline test

c--------------------------------------------- add endpoint and sort pth

cc      write (*,*) 'dsrdtab: add endpoint and start final sorting'
cc      write (*,*) 'nth=',nth
      pth(nth+1) = pmax*0.9999d0
      incth(nth+1)=1

      if (nth.gt.1) then
        do i=0,nth
          do j=nth,i,-1
cc            write (*,*) 'i,j',i,j
            if (pth(j).gt.pth(j+1)) then
              pres=pth(j)
              pth(j)=pth(j+1)
              pth(j+1)=pres
            endif
          enddo
        enddo
      endif

cc      write (*,*) 'dsrdtab: ending'


 999  continue

cc      write (*,*) 'Exiting dsrdtab'
      return
      end

      subroutine dsrdaddpt(wrate,pres,deltap)
c_______________________________________________________________________
c  add a point in rdrate table
c  input:
c    wrate - invariant annihilation rate (real, external)
c    pres - momentum of the point to add
c    deltap - scaling factor used in dsrdtab
c    pmax - maximum p used in dsrdtab (from common block)
c  common:
c    'dsrdcom.h' - included common blocks
c  used by dsrdtab
c  author: joakim edsjo (edsjo@physto.se)
c  modified: 01-01-31 paolo gondolo (paolo@mamma-mia.phys.cwru.edu) 
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 wrate
      external wrate
      integer i,j
      real*8 pres,deltap
      logical pexist
c-----------------------------------------------------------------------
      if (pres.le.0.0d0) return
c...check if p already exists [serious bug corrected pg 01-01-31]
c...Floating Point exception bug fixed 2008-09-09 JE
      pexist=.false.
      do i=1,nr
         if (pp(i).eq.0.d0.and.pres.eq.0.d0) then
            pexist=.true.
         else if
     &     (abs((pp(i)-pres/deltap)/max(pp(i),1.d-9)).lt.1.d-9) then
            pexist=.true.
         endif
      enddo
      if (pexist) return

      if (pres.gt.pmax) then
        if (rdprt.gt.0) write(*,*)
     &    'dsrdaddpt: pmax raised from ',pmax,' to ',pres
        pmax=pres
        if (nr+1.ge.nrmax) then
          write (rdluerr,*) 
     &      'error in dsrdaddpt: array capacity exceeded'
          write (rdluerr,*) '  for model ',rdtag
          write(rdluerr,*) 'Omega calculation stopping.'
          rderr=ibset(rderr,0)
          return
        endif
        nr=nr+1
        pp(nr)=pres/deltap
        yy(nr)=log(wrate(pres))
c...JE TEST. Delete following two lines later
c        write(*,*) 'pres=',pres,' wrate=',wrate(pres)
c        if (wrate(pres).lt.0.0d0) stop
        indx(nr)=nr
      else
        if (nr+1.ge.nrmax) then
          write (rdluerr,*) 
     &      'error in dsrdaddpt: array capacity exceeded'
          write (rdluerr,*) '  for model ',rdtag
          write(rdluerr,*) 'Omega calculation stopping.'
          rderr=ibset(rderr,0)
          return
        endif
        nr=nr+1
        pp(nr)=pres/deltap
        yy(nr)=log(wrate(pres))
c...JE TEST. Delete following two lines later
c        write(*,*) 'pres=',pres,' wrate=',wrate(pres)
c        if (wrate(pres).lt.0.0d0) stop
        do i=1,nr-2
          if (pp(indx(i)).lt.pp(nr).and.pp(indx(i+1)).gt.pp(nr)) then
            do j=nr-1,i+1,-1
              indx(j+1)=indx(j)
            enddo
            indx(i+1)=nr
          endif
        enddo
      endif

      return
      end



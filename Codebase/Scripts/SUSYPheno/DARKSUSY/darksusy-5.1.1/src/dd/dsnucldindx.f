      function dsnucldindx(a,z)
      implicit none
      integer dsnucldindx,a,z,ans,i,iseek,ilo,ihi,klo,khi,kk
      include 'dsnuclides.h'
      include 'dsio.h'

ccc      write (*,*) 'enter dsnucldindx.....',a,z
      if (inucld.eq.0) inucld=1

c...out of range

      if (z.lt.0.or.z.gt.nucldz(nnucld).or.
     &     a.lt.1.or.a.gt.nuclda(nnucld)) then

ccc         write (*,*) 'out of range.....'
         if (prtlevel.gt.0) write (*,*)
     &        'dsnucldindx: non-existent nuclide A=',a,' Z=',z
         ans=0
      
c...already at right place

      else if (nuclda(inucld).eq.a.and.nucldz(inucld).eq.z) then

ccc         write (*,*) 'already at right place.....'
         ans=inucld

c...bisection search

      else

ccc         write (*,*) 'start bisection search.....'
         iseek=1000*z+a
         klo=1
         khi=nnucld
 100     if (khi-klo.gt.1) then
            kk=(klo+khi)/2
            i=1000*nucldz(kk)+nuclda(kk)
            if (iseek.ge.i) then
               klo=kk
            else
               khi=kk
            endif
            goto 100
         endif
         ilo=1000*nucldz(klo)+nuclda(klo)
         ihi=1000*nucldz(khi)+nuclda(khi)
         if (iseek.eq.ilo) then
            ans=klo
         else if (iseek.eq.ihi) then
            ans=khi
         else
            ans=0
            if (prtlevel.gt.0) write (*,*)
     &           'dsnucldindx: non-existent nuclide A=',a,' Z=',z
         endif
      endif

 1000 continue
ccc      write (*,*) 'ans.....',ans
      inucld=ans
      dsnucldindx=inucld
ccc      write (*,*) 'nuclide is..... ',nucldsym(ans)
ccc      write (*,*) 'exit dsnucldindx.....'
      return
      end

      subroutine dshealpixpoints(n,ppmax,pp)
c
c     generate longitude and latitude of HEALPIX points
c
c     Input:
c     - n, integerm resolution
c     - ppmax, integer, maximum number of points >= 12*n**2
c
c     Output:
c     - pp, real(npts,2), longitude and latitude of each point in radians
c
      implicit none
      integer ppmax
      real*8 pp(ppmax,2)
      integer n,pmax,pmaxnorth,npolar,nequatorial,pminequatorial,
     &     p,q,ph,i,j,s
      real*8 nn,z,phi,pph,PIhalf,PI
      if (ppmax.lt.12*n*n) then
         write (*,*) 'dshealpixpoints: size ppmax of pp is < 12*n*n'
         stop
      endif
      PIhalf = 2.0 * datan(1.d0)
      PI=2.d0*PIhalf
      pmax = ppmax-1
      pmaxnorth = (pmax+1)/2+2*n-1
      npolar = n-1
      nequatorial = n+1
      pminequatorial = 2*n*(n-1)
      nn = real(n)
      do p=0,pmax
         if (p.gt.pmaxnorth) then
            q = pmax-p
         else
            q = p
         endif
         if (q.ge.pminequatorial) then
            ph = p-pminequatorial
            i = int(0.25d0*ph/nn)+n
            j = mod(ph,4*n)+1
            z = (4.d0-2.d0*i/nn)/3.d0
            s = mod(i-n+1,2)
            phi = PIhalf/nn*(j-0.5d0*s)
            !write (*,*) '+',n,nn,p,q,ph,i,j,z,s,phi/PI
         else
            pph = 0.5d0*(q+1)
            i = int(sqrt(pph-sqrt(real(int(pph)))))+1
            j = q+1-2*i*(i-1)
            if (p.gt.pmaxnorth) j=2*i-j+1
            z = 1-(i/nn)*(i/nn)/3.d0
            if (p.gt.pmaxnorth) z=-z
            phi = PIhalf/real(i)*(j-0.5d0)
            if (p.gt.pmaxnorth) i=4*n-i
            !write (*,*) '-',n,nn,p,q,pph,i,j,z,'na',phi/PI
         endif
         pp(p+1,1)=phi
         pp(p+1,2)=dasin(z)
      enddo
      return
      end

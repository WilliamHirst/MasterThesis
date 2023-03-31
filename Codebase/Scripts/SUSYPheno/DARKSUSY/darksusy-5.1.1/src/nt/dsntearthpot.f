


**********************************************************************
*** dsntearthpot  gives the gravitational potential inside and outside
*** of the earth as a function of the radius r (in meters).
*** input: radius in meters
*** output units: s^-1
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntearthpot(r)
      implicit none

      real*8 radius(42),pot(42),ppl,r,gn,dsntearthmass,
     &  dsntearthpotint
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1
      integer i,j

      logical first
      data first/.true./
      save first

      data radius/0.0d0,7140.0d0,128140.0d0,378140.0d0,628140.0d0,
     &  878140.0d0,1128140.0d0,1228140.0d0,1228150.0d0,1378140.0d0,
     &  1628140.0d0,1878140.0d0,2128140.0d0,2378140.0d0,2628140.0d0,
     &  2878140.0d0,3128140.0d0,3378140.0d0,3478140.0d0,3478150.0d0,
     &  3628140.0d0,3878140.0d0,4128140.0d0,4378140.0d0,4628140.0d0,
     &  4878140.0d0,5128140.0d0,5378140.0d0,5608140.0d0,5708140.0d0,
     &  5708150.0d0,5778140.0d0,5878140.0d0,5978140.0d0,5978150.0d0,
     &  6158140.0d0,6158150.0d0,6298140.0d0,6354140.0d0,6363140.0d0,
     &  6375140.0d0,6378140.0d0/

      save radius,pot

      if (first) then
        do i=1,42
          pot(i)=dsntearthpotint(radius(i))
        enddo
        first=.false.
      endif

      if (r.ge.6378.140d3) then
        dsntearthpot=-dsntearthmass(6378.140d3)*gn/r
        return
      endif

      if (r.le.0.0d0) then
        dsntearthpot=pot(1)
        return
      endif

      call dshunt(radius,42,r,j)

 20   ppl=(r-radius(j))/(radius(j+1)-radius(j))

      dsntearthpot=pot(j)*(1.0d0-ppl)+pot(j+1)*ppl

      return
      end

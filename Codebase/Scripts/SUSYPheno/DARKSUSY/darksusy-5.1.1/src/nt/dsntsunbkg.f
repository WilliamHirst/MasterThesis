***********************************************************************
*** dsntsunbkg calculats the differential background of muons cosmic
*** ray interactions in the sun's corona. the muon neutrino fluxes are from
*** g. ingelman and m. thunman, prd 54 (1996) 4385.
***   input:
***     emu - muon energy in gev
***     fltype = 1 - flux of muons
***              2 - contained event rate
***   output:
***     muon flux in units of gev^-1 km^-2(3) yr^-1
*** partly based on routines by l. bergstrom.
*** author: j. edsjo (edsjo@physto.se)
*** date: 1998-06-03
***********************************************************************

      real*8 function dsntsunbkg(emu,flt)
      implicit none

      real*8 emu

      real*8 e_mux
      real*8 dsff2,dsff3,a,b,eps,res
      external dsff2,dsff3
      integer flt
      integer fltype
      common/lbe_int2/e_mux,fltype
      fltype=flt
      eps=1.d-5
      a=emu
      e_mux=emu
      b=1.d7 !
      if (emu.ge.b) then
        dsntsunbkg=0.0d0
        return
      endif
c      call dgadap(a,b,dsff3,eps,res)
      call dgadap(log(a),log(b),dsff3,eps,res)
      res=res/1.0d15
c      call dsgauss1(dsff3,log(a),log(b),100,res,eps)
      if (flt.eq.1) then
        dsntsunbkg=res*3.15d17  ! convert to km^-2 yr^-1 gev^-1
      else
        dsntsunbkg=res*3.15d22  ! convert to km^-3 yr^-1 gev^-1
      endif
      return
      end

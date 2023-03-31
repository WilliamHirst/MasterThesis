***********************************************************************
*** dsntismbkg calculats the differential background of muons cosmic
*** ray interactions with the interstellar medium.
*** the muon neutrino fluxes are from
*** g. ingelman and m. thunman, hep-ph/9604286.
***   input:
***     emu - muon energy in gev
***     fltype = 1 - flux of muons
***              2 - contained event rate
***     rdelta = column density in units of nucleons / cm^2 kpc/cm
***   output:
***     muon flux in units of gev^-1 km^-2(3) yr^-1 sr^-1
*** partly based on routines by l. bergstrom.
*** author: j. edsjo (edsjo@physto.se)
*** date: 1998-09-20
***********************************************************************

      real*8 function dsntismbkg(emu,flt,rdelta)
      implicit none

      real*8 emu

      real*8 e_mux,rdelta,rd
      real*8 dsfff2,dsfff3,a,b,eps,res,dsntnuism
      external dsfff2,dsfff3,dsntnuism
      integer flt
      integer fltype
      common/lbe_int3/e_mux,rd,fltype
      fltype=flt
      rd=rdelta
      eps=1.d-5
      a=emu
      e_mux=emu
      b=1.d8 ! how high up it fluxes go.
      if (emu.ge.b) then
        dsntismbkg=0.0d0
        return
      endif
c      call dgadap(a,b,dsfff3,eps,res)
      call dgadap(log(a),log(b),dsfff3,eps,res)
      res=res/1.0d15
c      call dsgauss1(dsfff3,log(a),log(b),100,res,eps)
      if (flt.eq.1) then
        dsntismbkg=res*3.15d17  ! convert to km^-2 yr^-1 gev^-1
      else
        dsntismbkg=res*3.15d22  ! convert to km^-3 yr^-1 gev^-1
      endif
c      write(*,*) 'dsntnuism = ',dsntnuism(1000.0d0)
      return
      end

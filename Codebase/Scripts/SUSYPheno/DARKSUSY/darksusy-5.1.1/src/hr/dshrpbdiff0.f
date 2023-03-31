      real*8 function dshrpbdiff0(tp,solarmod,how)

**********************************************************************
*** function dshrpbdiff0 calculates the differential flux of
*** antiprotons for the antiproton kinetic energy tp as a result of
*** neutralino annihilation in the halo.
*** dshrpbdiff0 uses the unrescaled local density
*** see dshrpbardiff for rescaling the local density
***   input:
***     tp - antiproton kinetic energy in gev
***     solarmod - 0 no solar modulation
***                1 solar modulation a la perko
***     how - 1 calculate t_diff only for requested momentum
***           2 tabulate t_diff for first call and use table for
***             subsequent calls
***           3 as 2, but also write the table to disk as pbtd.dat
***             at the first call
***           4 read table from disk on first call, and use that for
***             subsequent calls
*** units: gev^-1 cm^-2 s^-1 sr^-1
*** author: joakim edsjo
*** date: 98-02-10
*** modified: joakim edsjo, edsjo@physto.se
*** modified: 98-07-13, 00-07-19 paolo gondolo
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      implicit none

      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dspbcom.h'
      include 'dsmpconst.h'

c----------------------------------------------------------------------
      real*8 tp,tpb,dqde,vtdiff,
     &  difftime,pp,epb,ppb,dshaloyield,phiap,dspbtd15,dspbtpb,
     &  dspbgalpropdiff
      integer istat,solarmod,how
      real*8 ep

c----------------------------------------------------------------------

c...Check if model is setup as it should
      if (.not.dshasetupcalled) then
        write(*,*) 'DS: Error in dshaloyield: dshasetup must be called',
     &    ' before any halo yield calculation'
        write(*,*) 'begins for every new model. Program stopping.'
        stop
      endif



c...solar modulation: calculate kinetic energy at heliosphere boundary
c...antiproton flux before modulation at kinetic energy tpb
      if (solarmod.eq.1) then  ! solar modulation
        ep=tp+m_p
        epb=ep+smod
        tpb=tp+smod
        pp=dsqrt(dabs(ep**2-m_p**2))
        ppb=dsqrt(dabs(epb**2-m_p**2))
c        tpb=dspbtpb(tp)
c        pp=dsqrt(2*m_p*tp+tp**2)
c        ppb=dsqrt(2*m_p*tpb+tpb**2)
c        epb=m_p+tpb
      else               ! no solar modulation
        tpb=tp
        pp=dsqrt(2*m_p*tp+tp**2)
        ppb=pp
        epb=m_p+tpb
      endif

c...  antiproton flux after modulation
      if (pbc.eq.'galprop') then
         dshrpbdiff0=dspbgalpropdiff(tpb)
         if (solarmod.eq.1) then
            dshrpbdiff0=dshrpbdiff0*pp**2/ppb**2
         endif
      else
         if (solarmod.eq.1) then
            phiap=dshaloyield(tpb,154,istat)*pp**2/ppb**2
         else
            phiap=dshaloyield(tp,154,istat)
         endif

c     JE Corr 03-01-21
         dqde = (rho0/hamwimp)**2*hasv*phiap/(4.0d0*pi)*0.5d0
         difftime = dspbtd15(tpb,how)
         vtdiff = (ppb/epb*2.99792458d10)*difftime*1.d15
         dshrpbdiff0 = vtdiff*dqde
      endif
      
      return
      end

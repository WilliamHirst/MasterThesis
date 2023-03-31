      real*8 function dsff2(x)
      implicit real*8 (a-h,o-z)
      real*8 e_mu,dsntnusun
      integer fltype
      common/lbe_int2/e_mu,fltype
      e_nu=x
      fluxnu=0.5d0*dsntnusun(e_nu)
      fluxnubar=fluxnu   ! approximately, not really true at higher e_nu
      r=fluxnu/fluxnubar
      sigma=(0.72*r+0.09)/(1.+r)
      a=(0.69+0.06*r)/(0.09+0.72*r)
c      alpha=2.0    ! g&s (8)  mod.to 2.5 by je?
c      epsilon=510. ! g&s      (650 gev?)
      alpha=2.6    ! ice values, je phd thesis
      epsilon=750. ! ice values, je phd thesis
c      alpha=2.25   ! rock values, je phd thesis
c      epsilon=522. ! rock values, je phd thesis
      avogadro=6.022
      if (fltype.eq.1) then ! flux of muons
        dsff2=avogadro*sigma/alpha*(fluxnu+fluxnubar)
c      dsff2=dsff2/(1.+e_mu/epsilon)*
c     1      (e_nu-e_mu+a/3.*(e_nu-e_mu*(e_mu/e_nu)**2)/e_nu**2)
        dsff2=dsff2/(1.+e_mu/epsilon)*
     1        (e_nu-e_mu+a/3.*(e_nu-e_mu*(e_mu/e_nu)**2)) !je corr 980603
        dsff2=dsff2*1.d-12 ! 23 from avo; -38 from sigma; 3 from alpha
      else ! contained events
        dsff2=avogadro*sigma*(fluxnu+fluxnubar)
        dsff2=dsff2*(1.0+a*(e_mu/e_nu)**2)
        dsff2=dsff2*1.d-15  ! 23 from avo; -38 from sigma
      endif
      return
      end

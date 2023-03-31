       real*8 function dsatm_mu(e_mu,c_th,flt)
c ***********************************************************************
c    gives muon flux from atmospheric neutrinos. uses dshonda.f and
c    dsgauss1.f.
c    based on the approximation in gaisser and stanev prd30 (1984) 985.
c    variables:
c     e_mu  muon energy in gev
c     c_th  cosine of zenith angle
c     fltype - 1 flux of muons in units of cm^-2 s^-1 sr^-1 gev^-1
c              2 cont. event rates in units of cm^-3 s^-1 sr^-1 gev^-1
c
c    output is dn/de_mu in muons per cm**2(3) per sec per sr per gev
c    l. bergstrom 1996-09-02
c    modified by j. edsjo (edsjo@physto.se)
c    date: jun-03-98
c
c ***********************************************************************
       implicit real*8 (a-h,o-z)
       real*8 e_mu,c_th,e_mux,c_thx
       real*8 dsff,dslnff,a,b,eps,res
       integer flt
       external dsff,dslnff
       integer nu_type,fltype
       common/lbe_int/c_thx,e_mux,nu_type,fltype
       fltype=flt
       eps=1.d-3
       a=e_mu
       e_mux=e_mu
       c_thx=c_th
       b=3900.d0 ! upper nuflux limit
       call dsgauss1(dslnff,log(a),log(b),res,eps,100)
       dsatm_mu=res
       return
       end

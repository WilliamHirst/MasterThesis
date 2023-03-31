
       real*8 function dsntdkka(aa,mx)
c--------------------------------------------------------------------------
c   solar fringe capture rate constant k_a^s for given mass number aa
c   for low-velocity population described by damour and krauss (1998)
c
c   lars bergstrom 1995-12-12
c--------------------------------------------------------------------------
       implicit none
       real*8 aa,mx,v_s2,m_a,betaplus,ahat_a,v_02,q_a,
     1       r_a2,eta_a,betaminus,dsntmoderf,sla,alphac,aofr_s,v_02c,
     2       alpha_dk
       include 'dshmcom.h'
       include 'dsmpconst.h'
       v_s2=(648.3d0)**2  ! (km/s)**2; see d-k eq (4.8)
       m_a=aa*(m_p+m_n)/2.0d0   ! nuclear mass
       alpha_dk=(30.d0)**2  ! see discussion after   (4.8)
       alphac=alpha_dk/(300000.)**2 ! use units where c=1
       betaplus=4.d0*mx*m_a/(mx+m_a)**2   !d-k eq (2.12)
       betaminus=4.d0*mx*m_a/(mx-m_a)**2  !d-k eq (2.12)
c       v_02=2./3.*vobs**2   ! vobs \sim 270 km/s gives v_0 \sim 220 km/s
       v_02=v_sun  ! solar velocity squared
       v_02c=v_02/(300000.)**2 ! use units where c=1
       r_a2=25.7d0*(0.3d0+0.91*m_a**0.3333d0)**2 !(2.26) in gev**-2
       q_a=3./2./m_a/r_a2                        !(2.26)
       ahat_a=mx*v_02c/2./q_a                     !(2.23)
c	   write(*,*) 'ahat_a: ',ahat_a
       eta_a=1.d0/sqrt(1.+ahat_a)*1.d0  ! assume v_s=v_0 in (2.23)
       aofr_s=sqrt(1.+ahat_a)*betaminus/v_02*(v_s2-alpha_dk/betaplus)
c       write(*,*) aofr_s,ahat_a,betaminus,v_02,v_s2,alpha_dk,betaplus
       aofr_s=sqrt(max(aofr_s,0.0d0)) ! (4.8) but put =0 for mx too heavy je
c	   write(*,*) 'aofrs: ',aofr_s
c	   write(*,*) 'eta_a: ',eta_a
       sla=dsntmoderf(aofr_s-eta_a)-dsntmoderf(aofr_s+eta_a)+
     1	   2.*dsntmoderf(eta_a)
       !   moderf is sqrt(pi)/2*erf(x)
       sla=sla*exp(-mx*alphac/2./q_a)*
     1       exp(-eta_a**2*ahat_a)/sqrt(3.141592)
       sla=sla/betaplus/(1.+ahat_a)/eta_a        !(2.25)
       dsntdkka=sla
c       write(*,*) aofr_s, dsntdkka
       return
       end

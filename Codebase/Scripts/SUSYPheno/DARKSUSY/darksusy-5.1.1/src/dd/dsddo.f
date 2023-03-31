      function dsddo(k,z,qsq)
c     k=0 Og, k=1 Ou, k=2 Od, ...... k=6 Ob
c     z=1 proton, z=2 neutron
c
c     This function is somehow related to the second moment of the
c     twist-2 parton matrix element. In the notation of Drees and
c     Nojiri, and dropping the momentum dependence and the index 2
c     indicating the second moment, one has Sigma+G=1, where Sigma
c     is the sum of the moments for each quark, Sigma=\sum q.
c     One has
c     <N|O^{(2)}_{q mu nu}|N> = (p_mu p_nu-m_N^2 g_{mu nu}/4)/m_N q,
c     and one needs \sum g_q q. From DN eq. (32),
c     with Op=O_{+}, Om=O_{-}, and a=alpha(m0^2)/alpha(Q0^2),
c     \sum g_q q = 
c       g_u [2 (15/31+Op a^{62/69})/5 + (Ou-Od)/2 + Sigma/10 a^{32/69}]
c     + g_d [2 (15/31+Op a^{62/69})/5 - (Ou-Od)/2 + Sigma/10 a^{32/69}]
c     + g_b [(15/31+Op a^{62/69})/5 - 2 Sigma/10 a^{32/69}]
c     + top quark (treated using Dn (17), (18), (33))
c
      implicit none
      integer k,z,i
      real*8 dsddo,qsq,q0sq,aa
      real*8 lam5sq,og0(2),o0(6,2)
      q0sq = 5.d0**2
      ! MTB1
      lam5sq = 0.146d0**2
      o0(1,1) = 0.312d0
      o0(2,1) = 0.197d0
      o0(3,1) = 0.d0
      o0(4,1) = 0.d0
      o0(5,1) = 0.d0
      o0(6,1) = 0.d0
      o0(1,2) = 0.164d0
      o0(2,2) = 0.345d0
      o0(3,2) = 0.d0
      o0(4,2) = 0.d0
      o0(5,2) = 0.d0
      o0(6,2) = 0.d0
      og0(1)=1.d0
      og0(2)=1.d0
      do i=1,6
         og0(1)=og0(1)-o0(i,1)
         og0(2)=og0(2)-o0(i,2)
      enddo
      if (qsq.gt.q0sq) then
         aa = log(q0sq/lam5sq)/log(qsq/lam5sq)
         if (k.eq.0) then
            dsddo = 16.d0/31.d0 
     &           - (16.d0/31.d0-og0(z))*exp(62.d0/29.d0*log(aa))
         else if (k.ge.1.and.k.le.6) then
            dsddo = 3.d0/31.d0 + 
     &           (16.d0/31.d0-og0(z))/5.d0*exp(62.d0/29.d0*log(aa)) +
     &           (o0(k,z)-(1.d0-og0(z))/5.d0)*exp(32.d0/69.d0*log(aa))
         else
            write (*,*) 'dsddo: invalid k=',k
            dsddo = 1.d10
         endif
      else
         if (k.eq.0) then
            dsddo = og0(z)
         else if (k.ge.1.and.k.le.6) then
            dsddo = o0(k,z)
         else
            write (*,*) 'dsddo: invalid k=',k
            dsddo = 1.d10
         endif         
      endif
      return
      end

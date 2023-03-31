      function dsqindx(kp1,kp2,kp3,kp4)
c_______________________________________________________________________
c  auxiliary function to dsvertx for quartic couplings 
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c=======================================================================
      implicit none
      integer dsqindx,kp1,kp2,kp3,kp4,x,i,j,k,l,max,imax,
     &     q,s,ii,jj,p(50)
c     with 31 particles involved in quartic couplings,
c     there will be 31*32/2=496  496*497/2=123256 4-particle vertices
c          knue,ke,knumu,kmu,knutau,ktau,ku,kd,kc,ks,kt,kb,
c          kgamma,kw,kz,kgluon,kh1,kh2,kh3,khc,
c          ksnue,kse1,kse2,ksnumu,ksmu1,ksmu2,ksnutau,kstau1,kstau2,
c          ksu1,ksu2,ksd1,ksd2,ksc1,ksc2,kss1,kss2,kst1,kst2,ksb1,ksb2,
c          kn1,kn2,kn3,kn4,kcha1,kcha2,kgluin,kgold0,kgoldc
      data p/
     &     0,0,0,0,0,0,0,0,0,0,0,0,
     &     1,2,3,4,5,6,7,8,
     &     9,10,11,12,13,14,15,16,17,
     &     18,19,20,21,22,23,24,25,26,27,28,29,
     &     0,0,0,0,0,0,0,30,31/
      i=p(kp1)
      j=p(kp2)
      k=p(kp3)
      l=p(kp4)
      if (i.eq.0.or.j.eq.0.or.k.eq.0.or.l.eq.0) then
         q = 0
      else
c     order indices
         s = 1
         max=i
         imax=1
         if (j.gt.max) then
            max=j
            imax=2
         endif
         if (k.gt.max) then
            max=k
            imax=3
         endif
         if (l.gt.max) then
            max=l
            imax=4
         endif
         if (imax.eq.2) then
            j=k
            k=l
            l=i
            i=max
            s=-s
         elseif (imax.eq.3) then
            k=i
            i=max
         elseif (imax.eq.4) then
            l=k
            k=j
            j=i
            i=max
            s=-s
         endif
         if (j.lt.l) then
            x=j
            j=l
            l=x
         endif
         if (i.eq.j.and.k.lt.l) then
            x=k
            k=l
            l=x
            s=-s
         endif
c     index in single array
         ii=k+i*(i-1)/2
         jj=l+j*(j-1)/2
         q=jj+ii*(ii-1)/2
         q=s*q
      endif
      dsqindx=q
      return
      end

      function dsrdlny(p,wrate)
c_______________________________________________________________________
c  logarithm of the invariant rate.
c  input:
c    p - initial cm momentum (real)
c    wrate - invariant annihilation rate (real, external)
c  called by dsrdtab.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      real*8 dsrdlny,p,wrate,tmp,ans
      external wrate
      tmp=wrate(p)
      ans=-1.d10
      if (tmp.gt.0.d0) ans=log(tmp)
      dsrdlny=ans
      end

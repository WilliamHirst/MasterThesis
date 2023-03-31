c====================================================================
c
c   this subroutine gives the real and imaginary parts of the 
c   amplitude of the process of neutralino annihilation into 
c   one photon and one z boson in the limit of vanishing relative 
c   velocity of the neutralino pair
c
c   p. ullio & l. bergstrom, phys. rev. d 57 (1998) 1962
c
c   the present version assumes equal sfermion mass eigenstates in
c   fermion - sfermion loop diagrams
c
c   imres: imaginary part
c   imfbxz: contribution to imres from diagram 1a - 1c divided by 
c     imres
c   imftxz: contribution to imres from diagrams 1e - 1h  divided by 
c     imres
c   imgbxz: contribution to imres from diagram 3a & 3b  divided by 
c     imres
c   reres: real part
c   refbxz: contribution to reres from diagram 1a - 1d  divided by 
c     reres
c   reftxz: contribution to reres from diagrams 1e - 1h  divided by 
c     reres
c   rehbxz: contribution to reres from diagram 2a - 2d  divided by 
c     reres
c   rehtxz: contribution to reres from diagrams 2e - 2h  divided by 
c     reres
c   regbxz: contribution from diagram 3a - 3f & 4a -4f divided by 
c     reres
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanzg(reres,imres)
      implicit none
      real*8 imresz,imfbxz,imftxz,imgbxz,reresz,refbxz
     &  ,reftxz,rehbxz,rehtxz,regbxz,reres,imres
      common/dsanzgcom/imresz,imfbxz,imftxz,imgbxz
     &  ,reresz,refbxz,reftxz,rehbxz,rehtxz,regbxz
      call dsanzgpar(imresz,imfbxz,imftxz,imgbxz,reresz,refbxz,reftxz,
     &  rehbxz,rehtxz,regbxz)
      imres=imresz
      reres=reresz
      return
      end










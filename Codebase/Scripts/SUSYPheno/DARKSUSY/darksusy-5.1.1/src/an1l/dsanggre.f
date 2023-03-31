c====================================================================
c
c   this subroutine gives the real part of the amplitude of the 
c   process of neutralino annihilation into two photons in the limit 
c   of vanishing relative velocity of the neutralino pair
c
c   l. bergstrom & p. ullio, nucl. phys. b 504 (1997) 27
c
c   reres: real part
c   refbxg: contribution from diagram 1a & 1b  divided by reres
c   reftxg: contribution from diagrams 1c & 1d  divided by reres
c   rehbxg: contribution from diagram 2a & 2b  divided by reres
c   rehtxg: contribution from diagrams 2c & 2d  divided by reres
c   regbxg: contribution from diagram 3a - 3c & 4a -4b divided by 
c     reres
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________


      subroutine dsanggre(reres)
      implicit none
      real*8 reres,reresg,refbxg,reftxg,rehbxg,rehtxg
     &  ,regbxg,imresg,imfbxg,imftxg,imgbxg
      common/anggcom/imresg,imfbxg,imftxg,imgbxg
     &  ,reresg,refbxg,reftxg,rehbxg,rehtxg,regbxg
      call dsanggrepar(reresg,refbxg,reftxg,rehbxg,rehtxg,regbxg)
      reres=reresg
      return
      end

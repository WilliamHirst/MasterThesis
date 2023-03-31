c====================================================================
c
c   this subroutine gives the imaginary part of the amplitude of the 
c   process of neutralino annihilation into two photons in the limit 
c   of vanishing relative velocity of the neutralino pair
c
c   l. bergstrom & p. ullio, nucl. phys. b 504 (1997) 27
c
c   imres: imaginary part
c   imfbxg: contribution from diagram 1a  divided by imres
c   imftxg: contribution from diagrams 1c & 1d  divided by imres
c   imgbxg: contribution from diagram 3a  divided by imres
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dsanggim(imres)
      implicit none
      real*8 imres,imresg,imfbxg,imftxg,imgbxg
     &  ,reresg,refbxg,reftxg,rehbxg,rehtxg,regbxg
      common/anggcom/imresg,imfbxg,imftxg,imgbxg
     &  ,reresg,refbxg,reftxg,rehbxg,rehtxg,regbxg
      call dsanggimpar(imresg,imfbxg,imftxg,imgbxg)
      imres=imresg
      return
      end

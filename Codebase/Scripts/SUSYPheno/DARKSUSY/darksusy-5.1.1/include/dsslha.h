*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsslha.h                              ***
***      these are common blocks for SLHA reading and writing        ***
c----------------------------------------------------------------------c
      double complex slhadata(nslhadata)
      common /dsslhacom/slhadata
      save /dsslhacom/

      integer prl ! print level
      common /dsslhaopt/prl
      save /dsslhaopt/
***                                                                 ***
************************ end of dsslha.h ******************************


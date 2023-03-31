* GetPara.F
* get the derived SUSY parameters from the internal arrays
* this file is part of FeynHiggs
* last modified 30 May 12 th

#include "externals.h"
#include "types.h"
#include "debug.h"


	subroutine FHGetPara(error,
     &    fv_, MASf_, UASf_,
     &    MCha_, UCha_, VCha_, MNeu_, ZNeu_,
     &    Deltab_, MGl_, MHtree_, SAtree_)
	implicit none
	integer error
	integer fv_
	RealType MASf_(6,4)
	RealType MCha_(2), MNeu_(4), MGl_
	ComplexType UASf_(6*6,4), UCha_(2*2), VCha_(2*2)
	ComplexType ZNeu_(4*4), Deltab_
	RealType MHtree_(4), SAtree_

#include "FH.h"
#define __SUBROUTINE__ "FHGetPara"

	integer i, t, ty

	error = 0
	CheckSf(error)
	if( error .gt. 0 ) return

	fv_ = fv

	do t = 1, 4
	  ty = t + (bBR - 4)*isB(t)
	  do i = 1, 6
	    MASf_(i,t) = MASf(i,ty)
	  enddo
	  do i = 1, 6*6
	    UASf_(i,t) = UASf_flat(i,ty)
	  enddo
	enddo

	do i = 1, 2
	  MCha_(i) = MCha(i)
	enddo

	do i = 1, 2*2
	  UCha_(i) = UCha_flat(i)
	enddo

	do i = 1, 2*2
	  VCha_(i) = VCha_flat(i)
	enddo

	do i = 1, 4
	  MNeu_(i) = MNeu(i)
	enddo

	do i = 1, 4*4
	  ZNeu_(i) = ZNeu_flat(i)
	enddo

	Deltab_ = Deltaf(bTR,3)
	MGl_ = MGl

	MHtree_(1) = MHtree(1)
	MHtree_(2) = MHtree(2)
	MHtree_(3) = MHtree(3)
	MHtree_(4) = MHtree(4)
	SAtree_ = SA
	end


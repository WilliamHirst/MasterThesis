* HEigensystem.F
* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 27 Sep 07 th
* adapted to darksusy 04 June 08 pg

************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

	subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	include 'dsdiag.h'

	integer p, q, j
	double precision red, off, thresh
	double precision delta, t, invc, s
	double complex x, y, Apq
	double precision ev(2,16)

	integer sweep
	common /nsweeps/ sweep

	double precision sq
	double complex c
	sq(c) = DBLE(c*DCONJG(c))

	if( n .gt. 16 ) then
	  print *, "HEigensystem: Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = DBLE(A(p,p))
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(A(p,q))
	    enddo
	  enddo
	  if( off .lt. 2D0**(-103) ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = sq(A(p,q))
	      if( sweep .gt. 4 .and. off .lt.
     &              2D0**(-103)*max(ev(2,p)**2, ev(2,q)**2) ) then
	        A(p,q) = 0
	      else
	        if( off .gt. thresh ) then
	          t = .5D0*(ev(2,p) - ev(2,q))
	          t = 1/(t + sign(sqrt(t**2 + off), t))

	          delta = t*off
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          s = t/invc
	          t = delta/(invc + 1)

	          Apq = A(p,q)

	          do j = 1, p - 1
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + s*(DCONJG(Apq)*y - t*x)
	            A(j,q) = y - s*(Apq*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = A(p,j)
	            y = A(j,q)
	            A(p,j) = x + s*(Apq*DCONJG(y) - t*x)
	            A(j,q) = y - s*(Apq*DCONJG(x) + t*y)
	          enddo

	          do j = q + 1, n
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + s*(Apq*y - t*x)
	            A(q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo

	          A(p,q) = 0

	          do j = 1, n
	            x = U(p,j)
	            y = U(q,j)
	            U(p,j) = x + s*(Apq*y - t*x)
	            U(q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "HEigensystem: Bad convergence in HEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end


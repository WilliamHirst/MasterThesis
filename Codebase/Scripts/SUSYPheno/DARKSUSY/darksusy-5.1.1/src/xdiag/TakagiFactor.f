* TakagiFactor.F
* computes the Takagi factorization of a complex symmetric matrix
* code adapted from the "Handbook" routines
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 27 Sep 07 th
* adapted to darksusy 04 June 08 pg

************************************************************************
** TakagiFactor factorizes a complex symmetric n-by-n matrix
** Input: n, A = n-by-n matrix, complex symmetric
** (only the upper triangle of A needs to be filled).
** Output: d = vector of diagonal values, U = transformation matrix
** these fulfill diag(d) = U^* A U^+ with U unitary.

	subroutine TakagiFactor(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	include 'dsdiag.h'

	integer p, q, j
	double precision red, off, thresh
	double precision sqp, sqq, fsq, t, invc, s
	double complex f, x, y
	double complex ev(2,16)

	integer sweep
	common /nsweeps/ sweep

	double precision sq
	double complex c
	sq(c) = DBLE(c*DCONJG(c))

	if( n .gt. 16 ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
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
	      sqp = sq(ev(2,p))
	      sqq = sq(ev(2,q))
	      if( sweep .gt. 4 .and.
     &            off .lt. 2D0**(-103)*max(sqp, sqq) ) then
	        A(p,q) = 0
	      else
	        if( off .gt. thresh ) then
	          f = sign(1D0, sqp - sqq)*
     &              (ev(2,q)*DCONJG(A(p,q)) +
     &                DCONJG(ev(2,p))*A(p,q))
	          fsq = DBLE(f)**2 + DIMAG(f)**2
	          if( fsq .eq. 0 ) then
	            f = 1
	            fsq = 1
	          endif
	          t = .5D0*abs(sqp - sqq)
	          t = 1/(t + sqrt(t**2 + fsq))

	          ev(1,p) = ev(1,p) + t*A(p,q)*DCONJG(f)
	          ev(2,p) = A(p,p) + ev(1,p)
	          ev(1,q) = ev(1,q) - t*A(p,q)*f
	          ev(2,q) = A(q,q) + ev(1,q)

	          invc = sqrt(fsq*t**2 + 1)
	          s = t/invc
	          t = t*fsq/(invc + 1)

	          do j = 1, p - 1
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + s*(DCONJG(f)*y - t*x)
	            A(j,q) = y - s*(f*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = A(p,j)
	            y = A(j,q)
	            A(p,j) = x + s*(DCONJG(f)*y - t*x)
	            A(j,q) = y - s*(f*x + t*y)
	          enddo

	          do j = q + 1, n
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + s*(DCONJG(f)*y - t*x)
	            A(q,j) = y - s*(f*x + t*y)
	          enddo

	          A(p,q) = 0

	          do j = 1, n
	            x = U(p,j)
	            y = U(q,j)
	            U(p,j) = x + s*(f*y - t*x)
	            U(q,j) = y - s*(DCONJG(f)*x + t*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    A(p,p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in TakagiFactor"

1	continue

* make the diagonal elements nonnegative

	do p = 1, n
	  d(p) = abs(A(p,p))
	  if( d(p) .gt. 2D0**(-52) .and. d(p) .ne. DBLE(A(p,p)) ) then
	    f = sqrt(A(p,p)/d(p))
	    do q = 1, n
	      U(p,q) = U(p,q)*f
	    enddo
	  endif
	enddo

	if( sort .eq. 0 ) return

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


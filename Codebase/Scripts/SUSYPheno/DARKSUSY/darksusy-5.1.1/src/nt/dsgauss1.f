      subroutine dsgauss1(f,a,b,result,eps,lambda)

      implicit real*8 (a-h,o-z)
      external f
      dimension w(12),x(12)
      data w
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      data x
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      result=0.d0
      eps=0.d0
      au=a
      c=(b-a)/(lambda)
      do 1 i = 1,lambda
      fi=i
      ao=a+fi*c
      c1=0.5*(au+ao)
      c2=c1-au
      s8=0.
      s16=0.
      do 2 j = 1,4
      u=x(j)*c2
    2 s8=s8+w(j)*(f(c1+u)+f(c1-u))
      do 3 j = 5,12
      u=x(j)*c2
    3 s16=s16+w(j)*(f(c1+u)+f(c1-u))
      result=result+c2*s16
      eps=eps+abs(c2*(s16-s8))
    1 au=ao
      return
      end

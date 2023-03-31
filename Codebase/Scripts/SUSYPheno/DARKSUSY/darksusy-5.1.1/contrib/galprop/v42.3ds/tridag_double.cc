
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * tridag_double.cc *                            galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/*****************************************************************************\
*                 *** I.Moskalenko (MPE,Garching) *** version of 04.08.99 *** *
* PURPOSE:                                                                    *
*   tridag(a,b,c,r,u,n) - solves for vector u[] the tridiagonal linear set:   *
*                 | b1 c1  0 ...                     | |  u1  | |  r1  |      *
*                 | a2 b2 c2 ...                     | |  u2  | |  r2  |      *
*                 |          ...                     |*|  ... |=|  ... |      *
*                 |          ... a(n-1) b(n-1) c(n-1)| |u(n-1)| |r(n-1)|      *
*                 |          ...   0    a(n)   b(n)  | | u(n) | | r(n) |      *
* INPUT parameters:                                                           *
*   (double) a[],b[],c[] are the input vectors (j=1..n),                      *
*   (int)    n - the length of the vectors.                                   *
* OUTPUT:                                                                     *
*   returns 0 when solution exists, in this case                              *
*      (double) u[] is a vector of solutions (j=1..n);                        *
*   returns 1 for a bad defined matrix;                                       * 
*   returns 2 when solution doesn't exist.                                    *
* NOTES:                                                                      *
* - The algorithm is adopted from the "Numerical Recipes" by W.H.Press et al. *
\*****************************************************************************/
// 19991124  replaced dynamic gam with static 

#include <stdio.h>

int tridag(double a[], double b[], double c[], double r[], double u[], int n) 
{
   int j,key;
   double bet;
#define NMAX 1000
   static double gam[NMAX];
   if(n>NMAX) {printf("\ntridag ... n>NMAX\n\n");return 3;}
   if(b[0] == 0) {  printf("\ntridag ... rewrite equations\n\n"); return (1); }

   bet = b[0];
   u[0]= r[0]/bet;
   for(j=1; j<n; j++) 
   {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if(bet == 0.) { printf("\ntridag ... tridag failed\n\n"); return (2); }
      u[j] = (r[j]-a[j]*u[j-1])/bet;
   }
   for(j=n-2; j>=0; u[j]=u[j]-gam[j+1]*u[j+1],j--);
   return (0);
}
/*
int tridag(double*,double*,double*,double*,double*,int);
main() {
   int i,n = 3;
   double a[3],b[3],c[3],r[3],u[3];

   for(i=0; i<n; i++) {
      a[i] = 0.;
      b[i] = 1.;
      c[i] = 1.;
      r[i] = 1.;
   }
   tridag(a,b,c,r,u,n);
   for(i=0,printf("\n"); i<n; printf(" %e ",u[i]),i++);
   printf("\n");
   return (0);
}
*/


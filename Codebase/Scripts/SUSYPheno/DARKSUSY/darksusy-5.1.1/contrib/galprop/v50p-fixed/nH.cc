
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nH.cc *                                       galprop package * 2/20/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// The formalism is described in Moskalenko I.V. et al. 2002, ApJ 565, 280
//
//1 The routine nH2 calculates H2 number density in mol/cm3
//1 nH2(R,Z) = epsilon0(R) *X /kpc2cm *exp(-ln2 (Z-Z0)^2/Zh^2), 
//1 where epsilon0(R) (K km s^-1 kpc^-1) -CO volume emissivity, and
//1 Z0(R), Zh(R) are taken from [B88]/Table 3/Cols 4,7,10;
//1 X = nH2/nCO =1.9e20 is the conversion factor taken from [SM96].
//1 [B88] Bronfman L. et al. 1988, ApJ 324, 248
//1 [SM96] Strong & Mattox 1996, A&A 308, L21
//
//2 The routine nHI calculates HI number density in cm^-3
//2 nHI(R,Z) = Y(R)/nGB *fZ(Z), 
//2 where Y(R) (cm^-3) - nHI number density taken from [GB76],
//2 nGB = 0.33 cm^-3 their average density of the disk at 4-8 kpc,
//2 fZ(Z) - Z-distribution: 
//2 0<R<8 [DL90], R>10 kpc [C86], with linear interpolation in 8-10 kpc.
//2 The total distribution is thus renormalized to the average surface density
//2 =6.2e20 cm^-2 at 4-8 kpc [DL90].
//2 References:
//2 [C86]  Cox, Kruegel, Mezger 1986, A&A 155, 380
//2 [DL90] Dickey & Lockman 1990, Ann. Rev. Astron. Astrophys. 28, 215
//2 [GB76] Gordon & Burton 1976, ApJ 208, 346/Table 1
//
//3 The routine nHII calculates HII number density in cm^-3
//3 Cordes J.M., et al. 1991, Nature 354, 121 / Equation (6) and Table 1
//3 Annular Gaussian preferred model for narrow component
//
//4 nH*_av routines make averaging over smaller steps.
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

using namespace std;//AWS20050624
#include<cmath>
#include"constants.h"


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH2(double Rkpc, double Zkpc)
{
   int i;
   double nH2_ = 0.0, fR,fZ0,fZh;                                              // [B88]/Table 3
   double R[18] ={ 0.00, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
                         6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,10.25},
          Y[18] ={ 0.00,  1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1
		          9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,  0.0},// (col.4)
	  Z0[18]={0.039,0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
                        -.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,-.020},
	  Zh[18]={0.077,0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
                        0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,0.147};
   double H2toCO = 1.9e20;              // [SM96]

   if(Rkpc > R[17]) return nH2_;
   for (i=0; i<17; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

   fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(Rkpc - R[i]); 
   fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   nH2_ =  fR * exp( -log(2.)*pow( (Zkpc-fZ0)/fZh, 2 ) )  *H2toCO/kpc2cm;
   return nH2_< 0. ? 0.: nH2_;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHI(double Rkpc, double Zkpc)
{
   int i;                                                             // Table 1 [GB76]
   double R[30] ={ 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,  // kpc, col.1
                   6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,10.0,10.5,11.0,
                  11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0},
          Y[30] ={ .10, .13, .14, .16, .19, .25, .30, .33, .32, .31,  // nHI, cm^-3
		   .30, .37, .38, .36, .32, .29, .38, .40, .25, .23,  // (col.3)
                   .32, .36, .32, .25, .16, .10, .09, .08, .06, .00};
   double fR, fZ,fZ1=0.,fZ2=0., R1,R2=R[29], Y1,Y2=Y[29];
   double nGB =0.33, nDL =0.57;       // cm^-3, disk density @ 4-8 kpc; [GB76], [DL90]
   double A1=0.395,     z1=0.212/2.,  // cm^-3, kpc; Z-distribution parameters from [DL90]
          A2=0.107,     z2=0.530/2.,
          B =0.064,     zh=0.403;

   for (i=0; i<29; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

   R1 = (R[i]+R[i+1])/2;   Y1 = Y[i];
   if(Rkpc < R1)
   {  
      if(i> 0)    { R2 = (R[i-1]+R[i])/2;   Y2 = Y[i-1]; }
      else        { R2 = R[0];              Y2 = Y[0];   }
   }
   else  if(i<28) { R2 = (R[i+1]+R[i+2])/2; Y2 = Y[i+1]; }

   fR = Y1 +(Y2 -Y1)/(R2 -R1)*(Rkpc -R1);                             // interpolation in R

   R2 = (R[28] +R[29]) /2;
   if(Rkpc > R2) fR = Y[28]*exp(-(Rkpc-R2)/3);                        // extrapolation in R

// calculation of Z-dependence
   if(Rkpc <10.)                                                      // [DL90]
      fZ1 =A1*exp(-log(2.)*pow(Zkpc/z1,2))+A2*exp(-log(2.)*pow(Zkpc/z2,2))+B*exp(-fabs(Zkpc)/zh);
   if(Rkpc > 8.) fZ2=nDL*exp(-pow(Zkpc /(0.0523*exp(0.11*Rkpc)), 2)); // [C86] IMOS20010220

   if(Rkpc <= 8.) fZ = fZ1;
   else
   {   if(Rkpc >=10.) fZ = fZ2;
       else fZ = fZ1 +(fZ2 -fZ1)/2.*(Rkpc -8.);                       // interp. [DL90] & [C86]
   }
   return fZ *fR/nGB;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHII(double Rkpc,double Zkpc)
{
   double fne1=0.025, H1=1.00, A1=20.0;
   double fne2=0.200, H2=0.15, A2= 2.0;
   double R2=4.0;
   double ne1 = fne1 * exp(-fabs(Zkpc)/H1) * exp (-pow( Rkpc    /A1, 2));
   double ne2 = fne2 * exp(-fabs(Zkpc)/H2) * exp (-pow((Rkpc-R2)/A2, 2));
   return ne1+ne2;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH2_av(double x,double y,double z,double dz,double dzz)
{  
   double R=sqrt(x*x + y*y);
   double nH2_av_=0.0;
   int nuse=0;

   for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nH2_av_+=nH2(R,zz);
      nuse++;
   }
   return nH2_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHI_av(double x,double y,double z,double dz,double dzz)
{
  
   double R=sqrt(x*x + y*y);
   double nHI_av_=0.0;
   int nuse=0;

   for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nHI_av_+=nHI(R,zz);
      nuse++;
   }
   return nHI_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHII_av(double x,double y,double z,double dz,double dzz)
{  
   double R=sqrt(x*x + y*y);
   double nHII_av_=0.0;
   int nuse=0;

   for(double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nHII_av_+=nHII(R,zz);
      nuse++;
   }
   return nHII_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/*
#include<stdio.h>
main()
{
   for(double R=0.; R<20.; R+=0.05) printf("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
     R,2*nH2(R,0.),2*nH2(R,0.1),2*nH2(R,0.2),nHI(R,0.),nHI(R,0.1),nHI(R,0.2),nHII(R,0.),nHII(R,0.1),nHII(R,0.2));
}
*/

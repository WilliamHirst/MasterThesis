
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_distribution.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include<math.h>

#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"

// cosmic ray source distribution: x,y,z in kpc

double source_distribution (double x,double y,double z)
{
   double r=sqrt(x*x + y*y);
   double alpha,beta,rmax,ro,rc,zscale,Value;
   double source_distribution_ =0.0;

   if (galdef.source_model==1)  //arbitrary parametrization
   {         
      alpha= galdef.source_parameters[1];
      beta = galdef.source_parameters[2]; 
      rmax = galdef.source_parameters[3]; 
      ro=8.5;
      zscale=0.200;  //nominal value

      source_distribution_ =pow(r/ro,alpha) * exp(-beta*(r-ro)/ro) * exp(-fabs(z)/zscale);
      if (r >= rmax) source_distribution_ = 0.0;
   }

   if (galdef.source_model==2)  //SNR:  Case and Bhattacharya 3rd Compton p437
   {
      alpha=1.69;
      beta=3.33;
      ro=8.5;
      zscale=0.200;  // nominal value;
      source_distribution_ =pow(r/ro,alpha) * exp(-beta*(r-ro)/ro) * exp(-fabs(z)/zscale);
   }

   if (galdef.source_model==3)  //pulsars: Taylor, Manchester and Lyne 1993 ApJSupp 88,259
   {
      rc=3.5;
      zscale=0.200; //nominal value, not for pulsars
      source_distribution_ =cosh(ro/rc)/cosh(r/rc)  * exp(-fabs(z)/zscale);
   }    


   if (galdef.source_model==5)  //Strong and Mattox gamma-ray distribution E> 100 MeV
   {
      zscale=0.200;
                     Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 7.4;
      source_distribution_ = Value * exp(-fabs(z)/zscale);
   }

   if (galdef.source_model==6)  //Strong and Mattox gamma-ray distribution E> 100 MeV R<15 kpc, zero beyond
   {
      zscale=0.200;
                     Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 0.0;
      source_distribution_ =  Value  * exp(-fabs(z)/zscale);
   }         


   if(galdef.n_spatial_dimensions==3)
   {
      for(int i_cr_source=0; i_cr_source<galdef.n_cr_sources; i_cr_source++)
      {
         double r2= pow(galdef.cr_source_x[i_cr_source]-x,2)
                   +pow(galdef.cr_source_y[i_cr_source]-y,2)
                   +pow(galdef.cr_source_z[i_cr_source]-z,2);

         double s=galdef.cr_source_L[i_cr_source]*exp(-r2/(2.*pow(galdef.cr_source_w[i_cr_source],2)));
         source_distribution_+=s;
// cout<<" source_distribution: x y z cr source r2 L s:"<<x<<" "<<y<<" "<<z<<" "<<i_cr_source<<" "<<r2<<" "<<s<<endl;
      }
   }   //  galdef.n_spatial_dimensions==3
// cout<<"source distribution R z model source "<< r <<" "<<z<<" "<<galdef.source_model<<" "<<source_distribution_<<endl;

   return source_distribution_;
}

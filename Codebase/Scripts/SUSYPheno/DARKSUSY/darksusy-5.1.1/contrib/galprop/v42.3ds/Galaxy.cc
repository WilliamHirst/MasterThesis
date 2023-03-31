
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galaxy.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"Galaxy.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::init( double r_min_, double r_max_, double dr_, 
                   double z_min_, double z_max_, double dz_) 
{
   cout<<" Galaxy  : initializing 2D\n";
   n_spatial_dimensions=2;       //2D    

   r_min=r_min_;
   r_max=r_max_;
   dr   =   dr_;
   z_min=z_min_;
   z_max=z_max_;
   dz   =   dz_;

   n_rgrid=(int)((r_max-r_min)/dr + 1.5);
   n_zgrid=(int)((z_max-z_min)/dz + 1.5);

   r=new double[n_rgrid];
   z=new double[n_zgrid];
     
   int ir,iz;    
   for(ir=0; ir<n_rgrid; ir++) r[ir]=r_min+ir*dr;
   for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz;

   n_HI    .init(n_rgrid, n_zgrid,1);
   n_H2    .init(n_rgrid, n_zgrid,1);
   n_HII   .init(n_rgrid, n_zgrid,1);
   B_field .init(n_rgrid, n_zgrid,1);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::init( double x_min_, double x_max_, double dx_, 
                   double y_min_, double y_max_, double dy_, 
                   double z_min_, double z_max_, double dz_) 
{  
   cout<<" Galaxy  : initializing 3D\n";
   n_spatial_dimensions=3;       //3D    

   x_min=x_min_;
   x_max=x_max_;
   dx   =   dx_;
   y_min=y_min_;
   y_max=y_max_;
   dy   =   dy_;
   z_min=z_min_;
   z_max=z_max_;
   dz   =   dz_;

   n_xgrid=(int)((x_max-x_min)/dx + 1.5);
   n_ygrid=(int)((y_max-y_min)/dy + 1.5);
   n_zgrid=(int)((z_max-z_min)/dz + 1.5);

   x=new double[n_xgrid];
   y=new double[n_ygrid];
   z=new double[n_zgrid];

   int ix,iy,iz;    
   for(ix=0; ix<n_xgrid; ix++) x[ix]=x_min+ix*dx; 
   for(iy=0; iy<n_ygrid; iy++) y[iy]=y_min+iy*dy; 
   for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz; 

   n_HI    .init(n_xgrid, n_ygrid, n_zgrid,1);
   n_H2    .init(n_xgrid, n_ygrid, n_zgrid,1);
   n_HII   .init(n_xgrid, n_ygrid, n_zgrid,1);	
   B_field .init(n_xgrid, n_ygrid, n_zgrid,1);

   SNR_cell_time  .init(n_xgrid, n_ygrid, n_zgrid,1);
   SNR_cell_phase .init(n_xgrid, n_ygrid, n_zgrid,1);

   SNR_electron_dg.init(n_xgrid, n_ygrid, n_zgrid,1); //AWS20010410
   SNR_nuc_dg     .init(n_xgrid, n_ygrid, n_zgrid,1); //AWS20010410
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::print()
{
   cout<<"Galaxy: n_spatial_dimensions="<<n_spatial_dimensions<<endl;
   cout<<"        n_ISRF_components   ="<<n_ISRF_components   <<endl;

   if(n_spatial_dimensions==2)
   {
      cout<<"r_min="<<r_min<<endl;
      cout<<"r_max="<<r_max<<endl;
      cout<<"dr   ="<<dr   <<endl;
      cout<<"z_min="<<z_min<<endl;
      cout<<"z_max="<<z_max<<endl;
      cout<<"dz   ="<<dz   <<endl;

      cout<<"n_rgrid="<<n_rgrid<<endl;
      cout<<"n_zgrid="<<n_zgrid<<endl;

      cout<<"r grid:"<<endl;
      for(int ir=0; ir<n_rgrid; cout<<r[ir++]<<" "); cout<<endl;
      cout<<"z grid:"<<endl;
      for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" "); cout<<endl;

      cout<<"n_HI[0][0].s[0]="<<n_HI.d2[0][0].s[0]<<endl;
   }//(n_spatial_dimensions==2

   if(n_spatial_dimensions==3)
   {
      cout<<"x_min="<<x_min<<endl;
      cout<<"x_max="<<x_max<<endl;
      cout<<"dx   ="<<dx   <<endl;
      cout<<"y_min="<<y_min<<endl;
      cout<<"y_max="<<y_max<<endl;
      cout<<"dy   ="<<dy   <<endl;
      cout<<"z_min="<<z_min<<endl;
      cout<<"z_max="<<z_max<<endl;
      cout<<"dz   ="<<dz   <<endl;

      cout<<"n_xgrid="<<n_xgrid<<endl;
      cout<<"n_ygrid="<<n_ygrid<<endl;
      cout<<"n_zgrid="<<n_zgrid<<endl;

      cout<<"x grid:"<<endl;
      for(int ix=0; ix<n_xgrid; cout<<x[ix++]<<" ");  cout<<endl;
      cout<<"y grid:"<<endl;
      for(int iy=0; iy<n_ygrid; cout<<y[iy++]<<" ");  cout<<endl;
      cout<<"z grid:"<<endl;
      for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" ");  cout<<endl;

      cout<<"n_E_gammagrid="<<n_E_gammagrid<<endl;
      cout<<"E_gamma grid:"<<endl;
      for(int iEgamma=0; iEgamma<n_E_gammagrid; cout<<E_gamma[iEgamma++]<<" ");  cout<<endl;

      cout<<"n_HI[0][0][0].s[0]="<<n_HI.d3[0][0][0].s[0]<<endl;
      cout<<"n_H2[0][0][0].s[0]="<<n_H2.d3[0][0][0].s[0]<<endl;
   } //(n_spatial_dimensions==3

}


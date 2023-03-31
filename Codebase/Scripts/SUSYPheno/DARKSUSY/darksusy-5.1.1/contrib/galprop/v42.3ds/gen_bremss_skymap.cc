
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_skymap.cc *                        galprop package * 4/29/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate bremsstrahlung skymaps 

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

int gen_bremss_skymap()
{
   cout<<" >>>>gen_bremss_skymap"<<endl;

   int stat=0;
   read_HIR();
   read_COR();

   Distribution emiss_av; // emissivity * gas density averaged for each ring
   Distribution emiss_HII;// emissivity * gas density averaged for each ring //IMOS20020429
   Distribution n_H_av;   //              gas density averaged for each ring
//                                       V axis 3 is ring number
   emiss_av. init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid);
   emiss_HII.init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid); //IMOS20020429
   n_H_av   .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,1);

   double Ro= 8.5; // Galactocentric distance of Sun, kpc
          Ro= 8.3; // to avoid discontinuity in maps
   double dd    =0.01; //  integration step in kpc
   double zzmax =1.0 ; //     maximum z for integration (i.e. limit of gas)
   double dtr=acos(-1.)/180.;
   int ir,ix,iy,iz;
   int i_Ring, n_Ring;

   n_Ring=galaxy.HIR.n_zgrid;
   for(int i_long=0; i_long<galaxy.n_long; i_long++)
   {
      for(int i_lat =0; i_lat<galaxy.n_lat; i_lat++)
      {
         double l=galaxy.long_min + i_long*galaxy.d_long;
         double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
         if(galdef.verbose>=1) cout<<"  gen_bremss_skymap l b ="<<l<<" "<<b<<endl;
         double sinb=sin(b*dtr);
         double cosb=cos(b*dtr);
         double sinl=sin(l*dtr);
         double cosl=cos(l*dtr);
         double d=0;
         int complete=0;
         while(complete==0)
         {
            d += dd;
            double zz=d*sinb;
            double RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl); // Galactocentric distance of point
            double costheta=(Ro-d*cosb*cosl)/RR;
            if(costheta> 1.0) costheta= 1.0;
            if(costheta<-1.0) costheta=-1.0;

            if(gcr[0].n_spatial_dimensions==2)
            {
               if(RR>galaxy.r_max)    complete=1;

               ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.0001);
               if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
               if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
               if(fabs(zz) > zzmax                  ) complete=1;

               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.0001);
               if(iz<0               ) { complete=1; iz=0;                } 
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }

               i_Ring=HIR_iRing(RR);//AWS20010316
//cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
            }
            if(gcr[0].n_spatial_dimensions==3)
            {
               double xx=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro
               double yy=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y since Z=X^Y system
               if(galdef.use_symmetry==1)
               {
                  xx=fabs(xx);
                  yy=fabs(yy);
                  zz=fabs(zz);
               }
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.0001);
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.0001);
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.0001);

               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
               if (zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
               if (fabs(zz) > zzmax                  ) complete=1;

               i_Ring=HIR_iRing(RR);
//cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3

            if(gcr[0].n_spatial_dimensions==2) n_H_av.d3[i_long][i_lat][i_Ring].s[0] +=
	       dd*kpc2cm *(galaxy.n_HI. d2[ir][iz]    .s[0] +2*galaxy.n_H2. d2[ir][iz]    .s[0]); //IMOS20020429

            if(gcr[0].n_spatial_dimensions==3) n_H_av.d3[i_long][i_lat][i_Ring].s[0] +=
               dd*kpc2cm *(galaxy.n_HI. d3[ix][iy][iz].s[0] +2*galaxy.n_H2. d3[ix][iy][iz].s[0]); //IMOS20020429

            for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) //IMOS20020429 whole loop
            {
               float delta,delta1;
               if(gcr[0].n_spatial_dimensions==2)
	       { 
                  delta =dd*kpc2cm *galaxy.bremss_emiss.        d2[ir][iz].s[iEgamma]
                     *(galaxy.n_HI. d2[ir][iz].s[0] +2*galaxy.n_H2. d2[ir][iz].s[0]);
                  delta1=dd*kpc2cm *galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]
		      *galaxy.n_HII.d2[ir][iz].s[0];                                          
               }
               if(gcr[0].n_spatial_dimensions==3) 
	       { 
                  delta =dd*kpc2cm *galaxy.bremss_emiss.        d3[ix][iy][iz].s[iEgamma]
                     *(galaxy.n_HI. d3[ix][iy][iz].s[0] +2*galaxy.n_H2. d3[ix][iy][iz].s[0]);
                  delta1=dd*kpc2cm *galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma]
                      *galaxy.n_HII.d3[ix][iy][iz].s[0];
               }
               emiss_av .d3[i_long][i_lat][i_Ring].s[iEgamma]+=delta;
               emiss_HII.d3[i_long][i_lat][i_Ring].s[iEgamma]+=delta1;
//cout<<"l b RR zz i_Ring Egamma emiss emiss_av "<<l<<" "<<b<<" "<<RR<<" "<<zz<<" "<<i_Ring
//<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma]<<" "
//<<emiss_av.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;     
            } // iEgamma
         } // complete==0
 
         for (i_Ring=0; i_Ring<n_Ring; i_Ring++)
         {
            if(n_H_av.d3[i_long][i_lat][i_Ring].s[0]>1.e17)
            {
               double w=       (galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+ 
                  2*galdef.X_CO*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0])
                                   /n_H_av.d3[i_long][i_lat][i_Ring].s[0];

               for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
               {
                  galaxy.bremss_skymap.d2[i_long][i_lat].        s[iEgamma] +=
                             emiss_av .d3[i_long][i_lat][i_Ring].s[iEgamma] *w
		            +emiss_HII.d3[i_long][i_lat][i_Ring].s[iEgamma];     // IMOS20020429

                  if(galdef.verbose==-100) // selectable debug
//if(n_H_av.d3[i_long][i_lat][i_Ring].s[0]<1.e14)
                     cout<<"l b i_Ring w HIR n_H_av skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w<<" "
                         <<galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_H_av.d3[i_long][i_lat][i_Ring].s[0]
                         <<" "<<galaxy.bremss_skymap.d2[i_long][i_lat].s[iEgamma]<<endl;
               }//iEgamma
            }//if
         }//i_Ring
// cout<<"i_long i_lat "<<i_long<<" "<<i_lat<<endl;

      }//lat
      cout<<"i_long  "<<i_long<<endl;
   }//long
   emiss_av .delete_array();
   emiss_HII.delete_array();//IMOS20020429
   n_H_av   .delete_array();

   if(galdef.verbose==-101)       // selectable debug
   {
      cout<<" bremsstrahlung skymap "<<endl;
      galaxy.bremss_skymap.print();
   } // galdef.verbose>=2
   cout<<" <<<< gen_bremss_skymap"<<endl;
   return stat;
}

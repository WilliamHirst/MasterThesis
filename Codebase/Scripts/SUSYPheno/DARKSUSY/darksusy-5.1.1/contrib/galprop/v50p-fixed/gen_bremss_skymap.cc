
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_skymap.cc *                     galprop package * 4/29/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

// generate pi0-decay skymaps 
// for HIR.10058 8 rings Ro=8.5 kpc 

/*
defined in gen_pi0_decay_skymap
int HIR_iRing(double RR)
{
   int i_Ring;
   i_Ring=0; // no ring RR<1.5 but assign to ring 0
   if (RR> 1.5) i_Ring=0;
   if (RR> 3.5) i_Ring=1;
   if (RR> 5.5) i_Ring=2;
   if (RR> 7.5) i_Ring=3;
   if (RR> 9.5) i_Ring=4;
   if (RR>11.5) i_Ring=5;
   if (RR>13.5) i_Ring=6;
   if (RR>15.5) i_Ring=7;
   return i_Ring;
}
 */

/* for HIR.10051 6 rings Ro=10 kpc
int HIR_iRing(double RR)
{
   int i_Ring;
   i_Ring=0;
   if (RR> 4.0) i_Ring=1;
   if (RR> 8.0) i_Ring=2;
   if (RR> 8.5) i_Ring=3; // fix, should be 10 but HI maps have Ro=10
   if (RR>12.0) i_Ring=4;
   if (RR>15.0) i_Ring=5;
   return i_Ring;
}
*/

int Galprop::gen_bremss_skymap()
{
   cout<<" >>>> gen_bremss_skymap"<<endl;
   int stat=0;

   read_HIR();
   read_COR();

   Distribution emiss_av; // emissivity * gas density averaged for each ring
   Distribution emiss_HII;// emissivity * HII density averaged for each ring //IMOS20020429
   Distribution emiss_HI ;// emissivity * HI  density averaged for each ring            //AWS20041214
   Distribution emiss_H2 ;// emissivity * H2  density averaged for each ring            //AWS20041214

   Distribution n_H_av;   //              gas density averaged for each ring
   Distribution n_HI_av;  //         NHI  gas density averaged for each ring            //AWS20041214
   Distribution n_H2_av;  //         NH2  gas density averaged for each ring            //AWS20041214

//                                       V axis 3 is ring number
   emiss_av. init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid);
   emiss_HII.init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid); //IMOS20020429
   emiss_HI .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid);  //AWS20041214
   emiss_H2 .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,galaxy.n_E_gammagrid);  //AWS20041214


   n_H_av   .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,1);
   n_HI_av  .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,1);                     //AWS20041214
   n_H2_av  .init(galaxy.n_long,galaxy.n_lat,galaxy.HIR.n_zgrid,1);                     //AWS20041214

   double Ro= 8.5; // Galactocentric distance of Sun, kpc
          Ro= 8.3; // to avoid discontinuity in maps
   double dd    =0.01; //     integration step in kpc 
   double zzmax =1.0 ; //     maximum z for integration (i.e. limit of gas)
   double dtr=acos(-1.)/180.;
   int ir,ix,iy,iz;
   int i_Ring,n_Ring;

   double w,w_HI,w_H2;                                                                 //AWS20041215

   n_Ring=galaxy.HIR.n_zgrid;
   for(int i_long=0; i_long<galaxy.n_long; i_long++)
   {
      for(int i_lat =0; i_lat <galaxy.n_lat; i_lat++)
      {
         double l=galaxy.long_min + i_long*galaxy.d_long;
         double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
         if(galdef.verbose>=1) cout<<"  gen_pi0_skymap l b ="<<l<<" "<<b<<endl;
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
               if(RR>galaxy.r_max) complete=1;       

               ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
               if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
               if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
               if(fabs(zz) > zzmax )     complete=1;

               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
               if(iz<0               ) { complete=1; iz=0; } 
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }

               i_Ring=HIR_iRing(RR);//AWS20010316
// cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==2

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
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }

               if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
               if(fabs(zz) > zzmax                  ) complete=1;

               i_Ring=HIR_iRing(RR);
//  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            }//particle.n_spatial_dimensions==3

            if(gcr[0].n_spatial_dimensions==2)
	      {
              n_H_av. d3[i_long][i_lat][i_Ring].s[0] += 
                   (galaxy.n_HI. d2[ir][iz]    .s[0] +2*galaxy.n_H2. d2[ir][iz]    .s[0])*dd*kpc2cm; //IMOS20020429

              n_HI_av.d3[i_long][i_lat][i_Ring].s[0] += 
                   (galaxy.n_HI. d2[ir][iz]    .s[0]                                    )*dd*kpc2cm; //AWS20041214

              n_H2_av.d3[i_long][i_lat][i_Ring].s[0] += 
               (                                      2*galaxy.n_H2. d2[ir][iz]    .s[0])*dd*kpc2cm; //AWS20041214
	      }

            if(gcr[0].n_spatial_dimensions==3)
	    {
              n_H_av.d3[i_long][i_lat][i_Ring].s[0] += 
               (galaxy.n_HI. d3[ix][iy][iz].s[0] +2*galaxy.n_H2. d3[ix][iy][iz].s[0])*dd*kpc2cm; //IMOS20020429


              n_HI_av.d3[i_long][i_lat][i_Ring].s[0] += 
               (galaxy.n_HI. d3[ix][iy][iz].s[0]                                    )*dd*kpc2cm; //AWS20041214

              n_H2_av.d3[i_long][i_lat][i_Ring].s[0] += 
               (                                  2*galaxy.n_H2. d3[ix][iy][iz].s[0])*dd*kpc2cm; //AWS20041214


	    }

//cout<<"dd d RR theta zz i_comp pi0_aniso_factor "<<dd<<" "<<d<<" "<<RR<<" "
//<<theta<<" "<<zz<<" "<<i_comp<<" "<<pi0_aniso_factor<<endl;

            for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) //IMOS20020429 whole loop
            {
               float delta,delta1;
               float delta_HI,delta_H2;                                                          //AWS20041214
            
               if(gcr[0].n_spatial_dimensions==2) 
	       { 
                  delta =  dd*kpc2cm *galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]
                     *(galaxy.n_HI. d2[ir][iz].s[0] +2*galaxy.n_H2. d2[ir][iz].s[0]);
                  delta1=  dd*kpc2cm *galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]
                      *galaxy.n_HII.d2[ir][iz].s[0];

                  delta_HI =dd*kpc2cm *galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]              //AWS20041214
                     *(galaxy.n_HI. d2[ir][iz].s[0]                                );            //AWS20041214
                  delta_H2 =dd*kpc2cm *galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]              //AWS20041214
                     *(                              2*galaxy.n_H2. d2[ir][iz].s[0]);            //AWS20041214

               }
               if(gcr[0].n_spatial_dimensions==3) 
	       { 
                  delta =  dd*kpc2cm *galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma] 
	             *(galaxy.n_HI. d3[ix][iy][iz].s[0] +2*galaxy.n_H2. d3[ix][iy][iz].s[0]);
                  delta1=  dd*kpc2cm *galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma] 
	              *galaxy.n_HII.d3[ix][iy][iz].s[0];
	         
                  delta_HI =dd*kpc2cm *galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma] 
	             *(galaxy.n_HI. d3[ix][iy][iz].s[0]                                    );    //AWS20041214
                  delta_H2 =dd*kpc2cm *galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma] 
	             *(                                  2*galaxy.n_H2. d3[ix][iy][iz].s[0]);    //AWS20041214

               }

               emiss_av. d3[i_long][i_lat][i_Ring].s[iEgamma]+= delta;
               emiss_HII.d3[i_long][i_lat][i_Ring].s[iEgamma]+= delta1;
               emiss_HI .d3[i_long][i_lat][i_Ring].s[iEgamma]+= delta_HI;                        //AWS20041214
               emiss_H2 .d3[i_long][i_lat][i_Ring].s[iEgamma]+= delta_H2;                        //AWS20041214

               if(galdef.verbose==-100) // selectable debug
                  cout<<"l b d RR zz i_Ring Egamma n_H_av n_HI_av  n_H2_av  emiss_av complete "<<l<<" "<<b<<" "<<d<<" "<<RR<<" "
                      <<zz<<" "<<i_Ring<<" "<<galaxy.E_gamma[iEgamma]<<" "
                      << n_H_av. d3[i_long][i_lat][i_Ring].s[0]<<" "
                      << n_HI_av.d3[i_long][i_lat][i_Ring].s[0]<<" "
                      << n_H2_av.d3[i_long][i_lat][i_Ring].s[0]<<" "
                      <<" "<<emiss_av.d3[i_long][i_lat][i_Ring].s[iEgamma]<<" "<<complete<<endl;     
            } // iEgamma
         } // complete==0
 
         for (i_Ring=0; i_Ring<n_Ring; i_Ring++)
         {
            if(n_H_av.d3[i_long][i_lat][i_Ring].s[0]>1.e17)
            {
                      w        =(galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]
       	  +2*galdef.X_CO[i_Ring]*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0])   //AWS20040227
                                    /n_H_av.d3[i_long][i_lat][i_Ring].s[0];

                      w_HI     = 0.;
	       	if(n_HI_av.d3[i_long][i_lat][i_Ring].s[0]>1.e17 )
                      w_HI     = galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]
                                   /n_HI_av.d3[i_long][i_lat][i_Ring].s[0];

	              w_H2     = 0.;
	       	if(n_H2_av.d3[i_long][i_lat][i_Ring].s[0]>1.e17 )
                      w_H2     =(                                                
       	  +2*galdef.X_CO[i_Ring]*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0])   //AWS20041214
                                   /n_H2_av.d3[i_long][i_lat][i_Ring].s[0];


               for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
               {
                  galaxy.bremss_skymap.d2[i_long][i_lat].        s[iEgamma] +=
                                emiss_av. d3[i_long][i_lat][i_Ring].s[iEgamma] *w
		               +emiss_HII.d3[i_long][i_lat][i_Ring].s[iEgamma];     //  IMOS20020429


                  galaxy.bremss_HIR_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] +=    //AWS20041214
                                emiss_HI     .d3[i_long][i_lat][i_Ring].s[iEgamma] *w_HI //AWS20041214
		               +emiss_HII    .d3[i_long][i_lat][i_Ring].s[iEgamma];      //AWS20041214


                  galaxy.bremss_H2R_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] +=    //AWS20041214
                                emiss_H2     .d3[i_long][i_lat][i_Ring].s[iEgamma] *w_H2;//AWS20041214
		          


                  if(galdef.verbose==-100) // selectable debug
		    {
                      //if(n_H_av.d3[i_long][i_lat][i_Ring].s[0]<1.e14)
                     cout<<"l b i_Ring w HIR n_H_av skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w<<" "
                         <<galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_H_av.d3[i_long][i_lat][i_Ring].s[0]
                         <<" "<<galaxy.bremss_skymap.d2[i_long][i_lat].s[iEgamma]<<endl;

                     cout<<"l b i_Ring w_HI HIR n_HI_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_HI<<" "
                         <<galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_HI_av.d3[i_long][i_lat][i_Ring].s[0]
                         <<" "<<galaxy.bremss_HIR_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;

                     cout<<"l b i_Ring w_H2 COR n_H2_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_H2<<" "
                         <<galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_H2_av.d3[i_long][i_lat][i_Ring].s[0]
                         <<" "<<galaxy.bremss_H2R_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;

		    }
               }//iEgamma
            }//if

         }//i_Ring
//cout<<"i_long i_lat "<<i_long<<" "<<i_lat<<endl;
      }//i_lat
      cout<<"i_long  "<<i_long<<endl;
   }//i_long
 
   emiss_av. delete_array();
   emiss_HII.delete_array();//IMOS20020429
   n_H_av.   delete_array();

   if(galdef.verbose==-101)       // selectable debug
   {
      cout<<" bremss skymap "<<endl;
      galaxy.bremss_skymap.print();
      cout<<" bremss_HIR skymap "<<endl;
      galaxy.bremss_HIR_skymap.print();
      cout<<" bremss_H2R skymap "<<endl;
      galaxy.bremss_H2R_skymap.print();

   }//galdef.verbose>=2

   cout<<" <<<< gen_bremss_skymap"<<endl;
   return stat;
}




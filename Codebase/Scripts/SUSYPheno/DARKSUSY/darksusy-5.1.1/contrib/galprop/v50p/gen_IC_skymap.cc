
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_IC_skymap.cc *                            galprop package * 4/20/2006
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate inverse Compton skymaps 
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

int Galprop::gen_IC_skymap()
{
   cout<<" >>>> gen_IC_skymap"<<endl;

   int stat=0;
   double Ro= 8.5; // Galactocentric distance of Sun, kpc
          Ro= 8.3; // to avoid discontinuity in maps
   double dd0    =0.1; // max integration step in kpc at b=90 deg.
   double ddmax = 0.5; // max integration step allowed
   double dtr=acos(-1.)/180.;
   int i,ir,ix,iy,iz;

//IMOS20060420 full calculation of anisoIC
   double RG=15.;                                  // kpc,radius of Galactic disk
   double RS=8.5;                                  // kpc,distance of the sun from Galactic center
   double factor, DENS;
   double *ISRF_over_nu, *Etarget;
   int ielectrons;
   Distribution electrons;

   if(galdef.IC_anisotropic==1) //IMOS20060420 full calculation of anisoIC
     {
// identify the electrons/positrons IMOS20030217
       if(galdef.n_spatial_dimensions==2) electrons.init(gcr[0].n_rgrid,                 gcr[0].n_zgrid, gcr[0].n_pgrid);
       if(galdef.n_spatial_dimensions==3) electrons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
       electrons = 0.;
       for(i=0, ielectrons=-1; i<n_species; i++)  
	 if(100==100*abs(gcr[i].Z)+gcr[i].A)
	   {
	     ielectrons=i;
	     electrons+=gcr[ielectrons].cr_density;
	     cout<<"  CR "<<gcr[ielectrons].name<<" found as species #"<<ielectrons<<endl;
	   }
       if(ielectrons==-1) { cout<<"CR electrons/positrons not found!"<<endl; electrons.delete_array(); return 1; }
       
       cout<<"gen_IC_skymap: generating ISRF_w and ISRF_over_nu"<<endl;
       Etarget     =new double[galaxy.ISRF[0].n_pgrid];
       ISRF_over_nu=new double[galaxy.ISRF[0].n_pgrid];
       
       for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
	 Etarget[inu]=h_planck * galaxy.nu_ISRF[inu] * erg_to_eV;  // target photon energy in eV
       
       factor=(eV_to_erg /h_planck) *log(galaxy.nu_ISRF[1] /galaxy.nu_ISRF[0]) *log(galdef.Ekin_factor);
     }
   
   for(int i_long=0; i_long<galaxy.n_long; i_long++)
   {
      for(int i_lat =0; i_lat<galaxy.n_lat; i_lat++)
      {
         double l=galaxy.long_min + i_long*galaxy.d_long;
         double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
         cout<<"  gen_IC_skymap: (l, b) = "<<l<<" "<<b<<endl;
         double sinb=sin(b*dtr);
         double cosb=cos(b*dtr);
         double sinl=sin(l*dtr);
         double cosl=cos(l*dtr);
         double dd=dd0/(fabs(sinb)+1.e-6);
         if(dd>ddmax) dd=ddmax;
         double d=0;
         int complete=0;

         while(complete==0)
         {
            d += dd;
            double Z=d*sinb;
            double R=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
            double costheta=(Ro-d*cosb*cosl)/R;
            if(costheta> 1.0) costheta= 1.0;
            if(costheta<-1.0) costheta=-1.0;
            double theta=acos(costheta);

            if(gcr[0].n_spatial_dimensions==2)
            {
               if(R>galaxy.r_max)                   complete=1;
               if(Z<galaxy.z_min || Z>galaxy.z_max) complete=1;
       
               ir=(int)((R-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
               iz=(int)((Z-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }

// cout<<"d R Z ir iz "<<d<<" "<<R<<" "<<Z<<" "<<ir<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==2

            if(gcr[0].n_spatial_dimensions==3)
            {
               double X=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro
               double Y=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y since Z=X^Y system
               if(galdef.use_symmetry==1)
               {
                  X=fabs(X);
                  Y=fabs(Y);
                  Z=fabs(Z);
               }
               ix=(int)((X-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((Y-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((Z-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }

//  cout<<"d R X Y Z ix iy iz "<<d<<" "<<R<<" "<<X<<" "<<Y<<" "<<Z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3

            for(int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)//IMOS20060420
            {
//cout<<"dd d R theta Z i_comp IC_aniso_factor "<<dd<<" "<<d<<" "<<R<<" "
//<<theta<<" "<<Z<<" "<<i_comp<<" "<<IC_aniso_factor<<endl;

	      if(galdef.IC_anisotropic==1) //IMOS20060420 full calculation of anisoIC
		{
		  if(gcr[0].n_spatial_dimensions==2)
		    for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
		      ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d2[ir][iz]    .s[inu]/galaxy.nu_ISRF[inu];

		  if(gcr[0].n_spatial_dimensions==3)
		    for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
		      ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d3[ix][iy][iz].s[inu]/galaxy.nu_ISRF[inu];
		}

	      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)//IMOS20060420
		{
		  double IC_aniso_emiss=0., IC_aniso_factor=1.;

		  if(galdef.IC_anisotropic==1 && i_comp<2) //IMOS20060420 full calculation of anisoIC// i_comp=2 -CMB
		    {
		      for(int ip=0; ip<gcr[ielectrons].n_pgrid; ip++)
			{
			  double sum=0.0;
			  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
			    {
			      if(galaxy.E_gamma[iEgamma] > 4.*Etarget[inu]*1.e-6 *pow(gcr[ielectrons].gamma[ip],2)  // IMOS20030217
				 /(1.+4.*Etarget[inu]*1.e-6/Mele*gcr[ielectrons].gamma[ip])) continue;

			      aic_cc(   0,1,Etarget[inu]*1.e-6/Mele,galaxy.E_gamma[iEgamma]/Mele,gcr[ielectrons].gamma[ip],RG,R,theta,fabs(Z),RS, DENS);

			      sum+= ISRF_over_nu[inu] 
				*aic_cc(1,1,Etarget[inu]*1.e-6/Mele,galaxy.E_gamma[iEgamma]/Mele,gcr[ielectrons].gamma[ip],RG,R,theta,fabs(Z),RS, DENS);
			    }

//cout<<"ir iz i_comp Ekin E_gamma sum cr_density "<<ir<<" "<<iz<<" "<<i_comp<<" "<<gcr[ielectrons].Ekin[ip]
//<<" "<<galaxy.E_gamma[iEgamma]<<" "<<sum<<" "<<gcr[ielectrons].cr_density.d2[ir][iz].s[ip]<<endl;

			  sum*=factor; // to avoid loss of accuracy

			  if(gcr[0].n_spatial_dimensions==2) 
			    IC_aniso_emiss+=sum *electrons.d2[ir][iz]    .s[ip]*gcr[ielectrons].Ekin[ip];

			  if(gcr[0].n_spatial_dimensions==3) 
			    IC_aniso_emiss+=sum *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip];
			}//ip

		      IC_aniso_factor=1.;

		      if(gcr[0].n_spatial_dimensions==2 && galaxy.IC_iso_emiss[i_comp].d2[ir][iz]    .s[iEgamma]!=0.) 
			IC_aniso_factor=IC_aniso_emiss/galaxy.IC_iso_emiss[i_comp].d2[ir][iz]    .s[iEgamma];

		      if(gcr[0].n_spatial_dimensions==3 && galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma]!=0.) 
			IC_aniso_factor=IC_aniso_emiss/galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma];
		    }

		  if(galdef.IC_anisotropic==2) IC_aniso_factor=IC_anisotropy_factor(galaxy.E_gamma[iEgamma],R,theta,fabs(Z),i_comp);//IMOS20060420
		  if(galdef.IC_anisotropic==3) IC_aniso_factor=1.0;            // anisotropic=isotropic skymap  IMOS20060420
		  if(galdef.IC_anisotropic) cout<<"ir iz i_comp E_gamma IC_aniso_factor-> "<<ir<<" "<<iz<<" "<<i_comp<<" "<<galaxy.E_gamma[iEgamma]<<" "<<IC_aniso_factor<<endl;

                  float delta;
                  if(gcr[0].n_spatial_dimensions==2) delta =dd*kpc2cm *galaxy.IC_iso_emiss[i_comp].d2[ir][iz]    .s[iEgamma];
                  if(gcr[0].n_spatial_dimensions==3) delta =dd*kpc2cm *galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma];

                  galaxy.IC_iso_skymap  [i_comp].d2[i_long][i_lat].s[iEgamma] +=delta;
                  if(galdef.IC_anisotropic) galaxy.IC_aniso_skymap[i_comp].d2[i_long][i_lat].s[iEgamma] +=delta*IC_aniso_factor;
		  
//cout<<"ir iz i_comp  E_gamma IC_iso_emiss "<<ir<<" "<<iz<<" "<<i_comp<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.IC_iso_emiss[i_comp].d2[ir][iz].s[iEgamma]<<endl;     
               }//iEgamma
           }//ISRF_components
         }//complete==0
      }//lat
   }//long

   // apply user-defined factors to components          AWS20050301
   for(int i_long=0; i_long<galaxy.n_long; i_long++)
     for(int i_lat =0; i_lat<galaxy.n_lat; i_lat++)
       for(int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)
	 for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	   {
	     galaxy.IC_iso_skymap  [i_comp].d2[i_long][i_lat].s[iEgamma] *= galdef.ISRF_factors[i_comp];
	     if(galdef.IC_anisotropic) galaxy.IC_aniso_skymap[i_comp].d2[i_long][i_lat].s[iEgamma] *= galdef.ISRF_factors[i_comp];
	   }

   if(galdef.verbose>=2)
   {
      for(int icomp=0; icomp<galaxy.n_ISRF_components; icomp++)
      {
         cout<<" inverse Compton isotropic skymap for ISRF component #"<<icomp<<endl;
         galaxy.IC_iso_skymap[icomp]  .print();
         if(galdef.IC_anisotropic) 
	   {
	     cout<<" inverse Compton anisotropic skymap for ISRF component #"<<icomp<<endl;
	     galaxy.IC_aniso_skymap[icomp].print();
	   }
      }
   } // galdef.verbose>=2
   cout<<" <<<< gen_IC_skymap"<<endl;
   if(galdef.IC_anisotropic==1) electrons.delete_array();  // IMOS20060420
   return stat;
}

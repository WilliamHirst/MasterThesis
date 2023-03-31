
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_DM_source.cc *                             galprop package * 9/09/2005 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine gen_DM_source calculates the source functions of the products of the
// dark matter (DM) particle annihilation [cm^-2 s^-1 sr^-1 MeV^-1].
// The routine can be used to calculate source function of positrons, electrons,
// and antiprotons.
// Use gen_DM_emiss to define gamma-ray emissivity (cm^-3 sr^-1 MeV^-1).
// The user must use the parameters DM_double0-9 and DM_int0-9 (galdef-file) to 
// specify the Galactic DM profile, branching, decay channels, and spectra (see 
// the template below). The DM profile is defined in the DM_profile routine.
// The profile is then averaged over the grid step (dR,dz) or (dx,dy,dz) with 
// a smaller step: normally 1/10 of the grid size.          IMOS20050912
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

extern "C" void rho_darksusy__(double*,double*,double*,double*);

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int gen_DM_source(Particle &particle)
{
   cout<<"gen_DM_source"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   double DMwidth,DMbranching,   // annihilation product distribution
     DMsecondary_spectrum,       // spectrum of secondaries from DM annihilation
     DME0,                       // delta function energy used for Green's function
     DMmass  =galdef.DM_double2, // DM particle mass
     DMcs_v  =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 
     dzz=0.01;                   // kpc, gas averaging step
   int stat=0;

 // define the spectra of annihilation products: positrons, electrons, antiprotons

   if(strcmp(particle.name,"DM_positrons")==0)
     {
       DMwidth     =galdef.DM_double3;
       DMbranching =galdef.DM_double4;
     }

   if(strcmp(particle.name,"DM_electrons")==0)
     {
       DMwidth     =galdef.DM_double5;
       DMbranching =galdef.DM_double6;
     }

   if(strcmp(particle.name,"DM_antiprotons")==0)
     {
       DMwidth     =galdef.DM_double7;
       DMbranching =galdef.DM_double8;
     }

   if(galdef.DM_int1==9) // Green's functions for all particles at the same Ekin
     {
       DMwidth     =galdef.Ekin_factor;
       DMbranching =1.;
       DME0=galdef.DM_double3;
       DMsecondary_spectrum=1./DME0/(1.-1./DMwidth);
     }

// assign the source function (2D)

   if(galaxy.n_spatial_dimensions==2)
     {
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
	       for(int ip=0; ip<particle.n_pgrid; ip++)
		 {
		   if(galdef.DM_int1==9) // Green's function (delta function)
		     {
		       if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
		       particle.secondary_source_function.d2[ir][iz].s[ip]+=pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz),2)*DMsecondary_spectrum*DMbranching/4./Pi*C;
		       continue;
		     }
		   if(particle.Etot[ip]<=DMmass)
		     {
		       DMsecondary_spectrum=exp(-pow((DMmass-particle.Etot[ip])/DMwidth,2));
		       particle.secondary_source_function.d2[ir][iz].s[ip]+= DMcs_v*
			 pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
			 /4./Pi*C*DMbranching*DMsecondary_spectrum/DMmass;
		     }
		 } // ip
	     }  //  iz
	 }  //  ir
     }  //  particle.n_spatial_dimensions==2
   
// assign the source function (3D)

   if(galaxy.n_spatial_dimensions==3)
     {
       for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int ip=0; ip<particle.n_pgrid; ip++)
		     {
		       if(galdef.DM_int1==9) // Green's function (delta function)
			 {
			   if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
			   particle.secondary_source_function.d3[ix][iy][iz].s[ip]+=pow(DM_profile_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz),2)*DMsecondary_spectrum*DMbranching/4./Pi*C;
			   continue;
			 }
		       if(particle.Etot[ip]<=DMmass) 
			 DMsecondary_spectrum=exp(-pow((DMmass-particle.Etot[ip])/DMwidth,2));
		       particle.secondary_source_function.d3[ix][iy][iz].s[ip]+= DMcs_v*
			 pow(DM_profile_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
			 /4./Pi*C*DMbranching*DMsecondary_spectrum/DMmass;
		     } //ip
		 }  //  iz
	     }  //  iy
	 }  //  ix
     }  //  particle.n_spatial_dimensions==3
 
 // test printout

   if(galdef.verbose>=2)
     {
       cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
       particle.secondary_source_function.print();
     }
   cout<<" <<<< gen_DM_source"<<endl;
   return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int gen_DM_emiss()
{
   cout<<"gen_DM_emiss"<<endl;
   double 
     DMmass      =galdef.DM_double2, // DM particle mass
     DMcs_v      =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 
     DMbranching =0.1,
     dzz=0.01;                       // kpc, gas averaging step
   int stat=0;

   galaxy.DM_emiss=0.;

// define the spectra of annihilation products: gammas
   
   if(galdef.n_spatial_dimensions==2)
     {
       cout<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions<<endl;
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
               for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		 {
		   if(galaxy.E_gamma[iEgamma]>DMmass) 
		     {
		       galaxy.DM_emiss.d2[ir][iz].s[iEgamma]=0;
		       continue;
		     }
		   galaxy.DM_emiss.d2[ir][iz].s[iEgamma]= DMcs_v *DMbranching
		     *pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
		     /galaxy.E_gamma[iEgamma];
		 }
	     }
	 }
     }
   if(galdef.n_spatial_dimensions==3)
     {
       cout<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions<<endl;
       for(int ix=0; ix<gcr[0].n_rgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_rgrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		     {
		       if(galaxy.E_gamma[iEgamma]>DMmass) 
			 {
			   galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=0;
			   continue;
			 }
		       galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=  DMcs_v *DMbranching
			 *pow(DM_profile_av(galaxy.x[ix], galaxy.y[ix], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
			 /galaxy.E_gamma[iEgamma];
		     }
		 }
	     }
	 }
     }
   cout<<" <<<< gen_DM_emiss"<<endl;
   return(stat);
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double DM_profile(double Xkpc, double Ykpc, double Zkpc)
{
  double R=sqrt(Xkpc*Xkpc+Ykpc*Ykpc+Zkpc*Zkpc),
    Rsun =8.5,                     //kpc, galactocentric distance of the solar system 
    Rc         =galdef.DM_double0, //core radius
    rho0       =galdef.DM_double1; //local DM mass density
  int profile_key =galdef.DM_int0; //profile type
  
  switch(profile_key)
    {
    case 0:   //NFW profile
      return(rho0*Rc/R*pow(1.+R/Rc,-2));
      
    case 1:   //isothermal profile
      return(rho0*(pow(Rc,2)+pow(Rsun,2))/(pow(Rc,2)+pow(R,2)));
      
    case 2:   //Evans profile
      return(rho0*pow(pow(Rc,2)+pow(Rsun,2),2)/(3.*pow(Rc,2)+pow(Rsun,2))
	     *(3.*pow(Rc,2)+pow(R,2))/pow(pow(Rc,2)+pow(R,2),2));
      
    case 3:   //alternative profile
      return(rho0*pow(Rc+Rsun,2)/pow(Rc+R,2));
      
    case 9:   //DarkSUSY profile (use only if the DarkSUSY and GALPROP combined) 
      rho_darksusy__(&Xkpc,&Ykpc,&Zkpc,&rho0);
      if(rho0<0.)
	{
	  cout<<"gen_DM_source: rho_darksusy() function is not defined"<<endl;
	  exit(0);
	}
      return(rho0);

    default:
      return(rho0);
    }
}
  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

 double DM_profile_av(double r,double z,double dr,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int Znuse=0, Rnuse=0;

     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       {
	 DM_profile_av_+=DM_profile(r,0,zz);
	 Znuse++;
       }
     
     for (double rr=r-dr/2.; rr<=r+dr/2.; rr+=dr/10.)
       {
	 if (rr<0.) continue;
	 DM_profile_av_+=DM_profile(rr,0,z);
	 Rnuse++;
       }
     return (DM_profile_av_/(Znuse+Rnuse));
   }
 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
 
 double DM_profile_av(double x,double y,double z,double dx,double dy,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int Znuse=0, Xnuse=0, Ynuse=0;
     
     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       {
	 DM_profile_av_+=DM_profile(x,y,zz);
	 Znuse++;
       }
     
     for (double xx=x-dx/2.; xx<=x+dx/2.; xx+=dx/10.)
       {
	 DM_profile_av_+=DM_profile(xx,y,z);
	 Xnuse++;
       }
     
     for (double yy=y-dy/2.; yy<=y+dy/2.; yy+=dy/10.)
       {
	 DM_profile_av_+=DM_profile(x,yy,z);
	 Ynuse++;
       }
     return DM_profile_av_/Znuse/Xnuse/Ynuse;
   }
 

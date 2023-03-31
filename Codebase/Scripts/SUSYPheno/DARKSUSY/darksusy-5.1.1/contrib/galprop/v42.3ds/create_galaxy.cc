
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_galaxy.cc *                            galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"


int create_galaxy()
{
   cout<<" >>>> create_galaxy"<<endl;

   int stat=0;
   double dzz=0.01; // average with this resolution in kpc

   if(galdef.n_spatial_dimensions==2)
      galaxy.init(galdef.r_min,  galdef.r_max, galdef.dr,  
                  galdef.z_min,  galdef.z_max, galdef.dz);
                   
   if(galdef.n_spatial_dimensions==3)
   {
      if(galdef.use_symmetry==1) galdef.x_min=galdef.y_min=galdef.z_min=0.; // IMOS20020419
      galaxy.init(galdef.x_min,  galdef.x_max, galdef.dx,
                  galdef.y_min,  galdef.y_max, galdef.dy,
                  galdef.z_min,  galdef.z_max, galdef.dz);
   }                  
// GAS and B-FIELD DISTRIBUTION

   galaxy.n_HI =1.e-6; //  non-zero values to avoid problem in energy loss logarithm
   galaxy.n_H2 =1.e-6;
   galaxy.n_HII=1.e-6;

   if(galdef.n_spatial_dimensions==2)
   {                                          // use average at y=0
      for(int ir=0; ir<galaxy.n_rgrid; ir++)
      {
         for(int iz=0; iz<galaxy.n_zgrid; iz++)
         {
            galaxy.n_HI. d2[ir][iz].s[0]=
               nHI_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
            galaxy.n_H2. d2[ir][iz].s[0]=
               nH2_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
            galaxy.n_HII.d2[ir][iz].s[0]=
               nHII_av(galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz); 
            galaxy.B_field.d2[ir][iz].s[0]=
               B_field_model(galaxy.r[ir],0.0,galaxy.z[iz],galdef.B_field_model);
         }
      }
   }
   if(galdef.n_spatial_dimensions==3)
   {
      for(int ix=0; ix<galaxy.n_xgrid; ix++)
      {
         for(int iy=0; iy<galaxy.n_ygrid; iy++)
         {
            for(int iz=0; iz<galaxy.n_zgrid; iz++)
            {
               galaxy.n_HI. d3[ix][iy][iz].s[0]=
                  nHI_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
               galaxy.n_H2. d3[ix][iy][iz].s[0]=
                  nH2_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
               galaxy.n_HII.d3[ix][iy][iz].s[0]=
                  nHII_av(galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);  
               galaxy.B_field.d3[ix][iy][iz].s[0]=
                  B_field_model(galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galdef.B_field_model);
            }
         }
      }
   }

   if(galdef.verbose>=1)
   {
      cout<<"galaxy.n_HI:   \n";galaxy.n_HI.print();
      cout<<"galaxy.n_H2:   \n";galaxy.n_H2.print();
      cout<<"galaxy.n_HII:  \n";galaxy.n_HII.print();
      cout<<"galaxy.B_field:\n";galaxy.B_field.print();
   }

// READING THE INTERSTELLAR RADIATION FIELD

   stat = read_isrf();
   if(gen_isrf_energy_density() !=0) return 1;

// STOCHASTIC SNR

   if(galdef.SNR_events==1) create_SNR();

// SKYMAP PARAMETERS

   galaxy.d_long=    galdef.d_long;
   galaxy.d_lat   =  galdef.d_lat;
   galaxy.long_min=  galdef.long_min;
   galaxy.long_max=  galdef.long_max;
   galaxy. lat_min=  galdef. lat_min;
   galaxy. lat_max=  galdef. lat_max;
   galaxy.n_long=(int)((galaxy.long_max-galaxy.long_min)/galaxy.d_long + 1.00001);
   galaxy.n_lat =(int)((galaxy. lat_max-galaxy. lat_min)/galaxy.d_lat  + 1.00001);

   cout<<"    gamma-ray, synchrotron skymap arrays: n_long,nlat="<<galaxy.n_long<<" "<<galaxy.n_lat<<endl;

// GAMMA RAYS

   galaxy.E_gamma_min=galdef.E_gamma_min;
   galaxy.E_gamma_max=galdef.E_gamma_max;
   galaxy.E_gamma_factor=galdef.E_gamma_factor;
   galaxy.n_E_gammagrid=(int)(log(galaxy.E_gamma_max/galaxy.E_gamma_min)/log(galaxy.E_gamma_factor)+1.001);
   galaxy.E_gamma=new double[galaxy.n_E_gammagrid];
   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
      galaxy.E_gamma[iEgamma]=exp(log(galaxy.E_gamma_min)+iEgamma*log(galaxy.E_gamma_factor));

   if(galdef.gamma_rays==1)
   {
      cout<<"   creating gamma-ray emissivity arrays"<<endl;
      galaxy.   IC_iso_emiss=new Distribution[galaxy.n_ISRF_components];

      if(galdef.n_spatial_dimensions==2)
      {
         galaxy.   bremss_emiss.        init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         galaxy.   bremss_ionized_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         for(int icomp=0; icomp<galaxy.n_ISRF_components; icomp++)
            galaxy.   IC_iso_emiss[icomp].init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         galaxy.pi0_decay_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);

	 galaxy.DM_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid); // DM: IMOS20050912
      }

      if(galdef.n_spatial_dimensions==3)
      {
         galaxy.   bremss_emiss.        init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         galaxy.   bremss_ionized_emiss.init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         for(int icomp=0; icomp<galaxy.n_ISRF_components; icomp++)
            galaxy.   IC_iso_emiss[icomp].init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
         galaxy.pi0_decay_emiss.init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);

	 galaxy.DM_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid); // DM: IMOS20050912
      }

// GAMMA-RAY SKYMAPS
      cout<<"   creating gamma-ray skymap arrays"<<endl;

      galaxy.bremss_skymap          .init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid);
      galaxy.bremss_ionized_skymap  .init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid);

      galaxy.IC_iso_skymap=new Distribution[galaxy.n_ISRF_components];
      for(int icomp=0; icomp<galaxy.n_ISRF_components; icomp++)
         galaxy.IC_iso_skymap[icomp].init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid);

      galaxy.IC_aniso_skymap=new Distribution[galaxy.n_ISRF_components];
      for(int icomp=0; icomp<galaxy.n_ISRF_components; icomp++)
         galaxy.IC_aniso_skymap[icomp].init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid);

      galaxy.pi0_decay_skymap.init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid);

      galaxy.DM_skymap.init(galaxy.n_long,galaxy.n_lat,galaxy.n_E_gammagrid); // DM: IMOS20050912 
   }

// SYNCHROTRON EMISSION

   if(galdef.synchrotron==1)
   {
      galaxy.nu_synch_min=   galdef.nu_synch_min;
      galaxy.nu_synch_max=   galdef.nu_synch_max;
      galaxy.nu_synch_factor=galdef.nu_synch_factor;
      galaxy.n_nu_synchgrid=(int)(log(galaxy.nu_synch_max/galaxy.nu_synch_min)/log(galaxy.nu_synch_factor)+1.001);

      galaxy.nu_synch=new double[galaxy.n_nu_synchgrid];
      cout<<"synchrotron frequency grid: ";
      for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
      {
         galaxy.nu_synch[inusynch]=exp(log(galaxy.nu_synch_min)+inusynch*log(galaxy.nu_synch_factor));
         cout<<galaxy.nu_synch[inusynch]<<" ";
      }
      cout<<endl;
      cout<<"   creating synchrotron emissivity array"<<endl;

      if(galdef.n_spatial_dimensions==2)
         galaxy.synchrotron_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,               galaxy.n_nu_synchgrid);

      if(galdef.n_spatial_dimensions==3)
         galaxy.synchrotron_emiss.init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,galaxy.n_nu_synchgrid);

      cout<<"   creating synchrotron skymap     array"<<endl;

      galaxy.synchrotron_skymap.init(galaxy.n_long,galaxy.n_lat,                 galaxy.n_nu_synchgrid);

   }

   cout<<"   creating ionization rate array"<<endl;

   if(galdef.n_spatial_dimensions==2)
      galaxy.ionization_rate.init(galaxy.n_rgrid,galaxy.n_zgrid,                1);

   if(galdef.n_spatial_dimensions==3)
      galaxy.ionization_rate.init(galaxy.n_xgrid,galaxy.n_ygrid, galaxy.n_zgrid,1);

   cout<<"      galaxy.print"<<endl;
   galaxy.print();
   cout<<" <<<< create_galaxy"<<endl;
   return stat;
}

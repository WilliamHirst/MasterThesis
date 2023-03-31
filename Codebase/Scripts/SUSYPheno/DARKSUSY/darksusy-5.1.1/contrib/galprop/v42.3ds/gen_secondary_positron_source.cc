
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_positron_source.cc *            galprop package * 6/02/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine to calculate the secondary positron and electron source function. 
//
// CR density gcr.cr_density is in c/4pi * n(E) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// The routine PP_MESON written in FORTRAN-77 is designed to calculate
// sec. positron (or sec. electron) production spectrum vs. energy (barn/GeV). 
// Positron/electron energy and total nucleus momentum (GeV) are used as input
// parameters as well as beam and target nuclei atomic numbers.
//
// The positron/electron source function [cm^-2 s^-1 sr^-1 MeV^-1] as used in
// galprop is defined as following:
//                 ___      ___  
//          beta c \        \    /             d sigma_ij(Etot,p')
// q(Etot)= ------ /__  n_i /__  \ dp' n_j(p') ------------------ ,
//           4pi  i=H,He     j   /                    dEtot  
// 
// where n_i is the gas density, d sigma_ij(Etot,p')/d Etot is
// the production cross section, n_j(p') is the CR species density, 
// and p' is the total momentum of a nucleus.
// Substitution of dp' with d(log Ekin') gives:
//                ___                            ___  
//           c A  \        /                     \              d sigma_ij(Etot,Ekin')
// q(Etot) = ---  /__ n_i  \ d(log Ekin')  Ekin' /__  n_j(Ekin') --------------------
//           4pi i=H,He    /                      j                      d Etot  
//                             ___     ___       ___
//           c A               \       \         \              d sigma_ij(Etot,Ekin')
//         = --- /\(log Ekin') /__ n_i /__ Ekin' /__ n_j(Ekin') --------------------- ,
//           4pi              i=H,He   Ekin       j                      d Etot             
// 
// where /\=Delta, and we used dp'=1/beta A Ekin' d(log Ekin').
// Since positrons/electrons are suggested massless Etot=p.
// 
// To transfer to units cm^2/MeV we need a factor= 1.0e-24 *1.0e-3.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int gen_secondary_positron_source(Particle &particle)
{
   cout<<"gen_secondary_positron_source"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   int key1=-9999;
   if(strcmp(particle.name,"secondary_positrons")==0) key1=+3;  // e+ from pi+ decay including kaons
   if(strcmp(particle.name,"secondary_electrons")==0) key1=-3;  // e- from pi- decay including kaons
   if(key1==-9999) { cout<<"invalid particle "<<particle.name<<endl;  return 2; }

   int stat=0, iprotons=-1, iHelium =-1,  NA1, NA2;
   float cs_p_HI, cs_p_He, cs_He_HI, cs_He_He;
   Distribution protons;                 // IMOS20000606.1

// identify the CR protons               // IMOS20000606.2
   if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   protons=0.;
   for(int i=0; i<n_species; i++)  
      if(101==100*gcr[i].Z+gcr[i].A)
      {
         iprotons=i;
 	 protons+=gcr[iprotons].cr_density;
         cout<<"  CR protons found as species #"<<iprotons<<endl;
      }
   if(iprotons==-1) { cout<<"CR protons not found!"<<endl;  return 1; }
 
// identify the CR Helium
   for(int i=0; i<n_species; i++)  if(204==100*gcr[i].Z+gcr[i].A)  iHelium =i;
   if(iHelium ==-1) { cout<<"CR Helium  not found!"<<endl;  return 1; }
   cout<<"  CR Helium  found as species #"<<iHelium <<endl;

   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {
         NA1=1;   NA2=1;  //  beam+target: p+HI
         cs_p_HI =pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1 );

         NA1=1;   NA2=4;  //  beam+target: p+He
         cs_p_He =pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1 );

         NA1=4;   NA2=1;  //  beam+target: He+HI
         cs_He_HI=pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1 );

         NA1=4;   NA2=4;  //  beam+target: He+He
         cs_He_He=pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1 );

         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               {
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                     (galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
                  *( (cs_p_HI +cs_p_He *galdef.He_H_ratio) 
                     *protons.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip]                  // IMOS20000606.3
//                     *gcr[iprotons].cr_density.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip]
                    +(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                     *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A );
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2
 
         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=
                        (galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                     *( (cs_p_HI  +cs_p_He *galdef.He_H_ratio) 
                        *protons.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip]           // IMOS20000606.4
//                        *gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip]
                       +(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                        *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A );
                  }  //  iz
               }  //  iy
            }  //  ix
         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  iEgamma
 
   double factor=1.e-24 *1.e-3 *c *log(galdef.Ekin_factor); // transformation to cm2/MeV and constant factors
   particle.secondary_source_function *= factor;

   protons.delete_array();                  // IMOS20000606.5

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
   }
   cout<<" <<<< gen_secondary_positron_source"<<endl;
   return stat;
}

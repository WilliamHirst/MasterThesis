
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_synch_emiss.cc *                          galprop package * 5/02/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate synchrotron emissivity
/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1  * cm^-3 MeV^-1]

emissivity (erg cm^-3 sr^-1 Hz^-1 s^-1)=
(c/4pi)*integral[sigma{Eelectron,nu}  ) n(E)E dlog(E)]

emmissivity fron synchrotron_cc in erg s^-1 Hz-1
factor=  log(Ekin_factor)
*/

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

int gen_synch_emiss()
{
   cout<<"gen_synchs_emiss"<<endl;
   cout<<"generating synchrotron    emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, i, ielectrons, ipositrons;

// identify the primary electrons 
   ielectrons=-1;
   for(i=0; i<n_species; i++)  if(strcmp(gcr[i].name,"primary_electrons")==0)  ielectrons=i;
   if(ielectrons==-1) { cout<<"primary electrons not found!"<<endl; return 1; }
   cout<<"  primary electrons found as species #"<<ielectrons<<endl;

   Distribution electrons;                                            // IMOS20020429 >>>
   if(galdef.n_spatial_dimensions==2) 
      electrons.init(gcr[ielectrons].n_rgrid, gcr[ielectrons].n_zgrid, gcr[ielectrons].n_pgrid);
   if(galdef.n_spatial_dimensions==3) 
      electrons.init(gcr[ielectrons].n_xgrid, gcr[ielectrons].n_ygrid, gcr[ielectrons].n_zgrid, gcr[ielectrons].n_pgrid);
   electrons += gcr[ielectrons].cr_density;

// identify the secondary electrons
   if(galdef.secondary_electrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"secondary_electrons")==0) ipositrons=i;
      if(ipositrons!=-1)
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  secondary electrons found as species #"<<ipositrons<<endl;
      }
   }
// identify the secondary positrons
   if(galdef.secondary_positrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"secondary_positrons")==0) ipositrons=i;
      if(ipositrons!=-1) 
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  secondary positrons found as species #"<<ipositrons<<endl;
      }
   }                                                                  // IMOS20020429 <<<

// identify the DM positrons  IMOS20050912
   if(galdef.DM_positrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"DM_positrons")==0) ipositrons=i;
      if(ipositrons!=-1) 
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  DM positrons found as species #"<<ipositrons<<endl;
      }
   }

// identify the DM electrons  IMOS20050912
   if(galdef.DM_electrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"DM_electrons")==0) ipositrons=i;
      if(ipositrons!=-1)
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  DM electrons found as species #"<<ipositrons<<endl;
      }
   }

   for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
   {
      for(int ip=0; ip<gcr[ielectrons].n_pgrid; ip++)
      {

         if(gcr[ielectrons].n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[ielectrons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
               {                                                                                    
		  double synch_per_electron =synchrotron_cc(gcr[ielectrons].gamma[ip], // Tesla ->Gauss
                     galaxy.nu_synch[inusynch], galaxy.B_field.d2[ir][iz].s[0]*10000.);// IMOS20020429

                  galaxy.synchrotron_emiss.d2[ir][iz].s[inusynch] +=synch_per_electron 
	             *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip];             // IMOS20020429
//cout<<"ir iz  E_gamma bremss_emiss "<<ir<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]<<endl;  
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2
 
         if(gcr[ielectrons].n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[ielectrons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[ielectrons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
                  {
                     double synch_per_electron =synchrotron_cc(gcr[ielectrons].gamma[ip],     // Tesla ->Gauss
                        galaxy.nu_synch[inusynch], galaxy.B_field.d3[ix][iy][iz].s[0]*10000.);// IMOS20020429

                     galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch] += synch_per_electron
	                *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip];             // IMOS20020429
                  }  //  iz
               }  //  iy
            }  //  ix
         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  inusynch

   double factor= log(galdef.Ekin_factor);
   galaxy.synchrotron_emiss *= factor;

   if(galdef.verbose>=2)
   {
      cout<<"   synchrotron    emissivity "<<endl;
      galaxy.synchrotron_emiss.print();
   }  //  galdef.verbose
   electrons.delete_array();  // IMOS20020429
   cout<<" <<<< gen_synch_emiss"<<endl;
   return stat;
}

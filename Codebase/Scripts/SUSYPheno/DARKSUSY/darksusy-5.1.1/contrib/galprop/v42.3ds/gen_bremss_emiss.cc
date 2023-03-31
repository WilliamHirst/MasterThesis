
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_emiss.cc *                         galprop package * 4/29/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate bremsstrahlung emissivity
/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1  * cm^-3 MeV^-1]
emissivity (cm^-3 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,Eelectron}  ) n(E)E dlog(E)]
cross section from bremss_spec in barns MeV^-1
factor= 1.0e-24* log(Ekin_factor)
*/

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

int gen_bremss_emiss()
{
   cout<<"gen_bremss_emiss"<<endl;
   cout<<"generating bremsstrahlung emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, i, IZ1, Ne1, ielectrons, ipositrons;
   float cs_HII, cs_HI, cs_He;

// identify the primary electrons 
   ielectrons=-1;
   for(i=0; i<n_species; i++) if(strcmp(gcr[i].name,"primary_electrons")==0) ielectrons=i;
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

// identify the DM positrons   IMOS20050912
   if(galdef.DM_positrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"DM_positrons")==0) ipositrons=i;
      if(ipositrons!=-1) 
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  DM positrons found as species #"<<ipositrons<<endl;
      }
   }

// identify the DM electrons   IMOS20050912
   if(galdef.DM_electrons)
   {
      for(ipositrons=-1, i=0; i<n_species; i++) if(strcmp(gcr[i].name,"DM_electrons")==0) ipositrons=i;
      if(ipositrons!=-1)
      { 
         electrons += gcr[ipositrons].cr_density;
         cout<<"  DM electrons found as species #"<<ipositrons<<endl;
      }
   }

   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
   {
      for(int ip=0;ip<gcr[ielectrons].n_pgrid;ip++)
      {
         IZ1=1; Ne1=0; //              ionized hydrogen
         cs_HII=bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);
	 IZ1=1; Ne1=1; //              1-electron atom 
         cs_HI =bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);
         IZ1=2; Ne1=2; //              2-electron atom 
         cs_He =bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);

         if(gcr[ielectrons].n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[ielectrons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
               {
                  galaxy.bremss_emiss.        d2[ir][iz].s[iEgamma] +=
                     (cs_HI +cs_He*galdef.He_H_ratio) *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
                  galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma] +=
                      cs_HII                          *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
//cout<<"ir iz  E_gamma bremss_emiss "<<ir<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]<<endl;
               }//iz
            }//ir
         }//particle.n_spatial_dimensions==2
 
         if(gcr[ielectrons].n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[ielectrons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[ielectrons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
                  {
                     galaxy.bremss_emiss.        d3[ix][iy][iz].s[iEgamma] +=
                        (cs_HI +cs_He*galdef.He_H_ratio) *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
                     galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma] +=
                         cs_HII                          *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
//cout<<"ix iy  iz  E_gamma bremss_emiss "<<ix<<" "<<iy<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma]<<endl; 
                  }//iz
               }//iy
            }//ix
         }//particle.n_spatial_dimensions==3
      }//ip
   }//iEgamma

   double factor=1.0e-24* log(galdef.Ekin_factor);
   galaxy.bremss_emiss        *= factor;
   galaxy.bremss_ionized_emiss*= factor; //IMOS20020429

   if(galdef.verbose>=2)
   {
      cout<<"   bremsstrahlung emissivity "<<endl;
      galaxy.bremss_emiss.print();
   }//galdef.verbose
   electrons.delete_array();  // IMOS20020429
   cout<<" <<<< gen_bremss_emiss"<<endl;
   return stat;
}

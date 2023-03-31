
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_pi0_decay_emiss.cc *                      galprop package * 4/23/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1  * cm^-3 MeV^-1]

emissivity (cm^-3 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,p_beam }  ) n(E)E dlog(E)]

pp_meson has Egamma in GeV, beam momentum in GeV
The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
BUT UNITS OF density/momentum = flux/(KE/nucleon)..... CHECK CAREFULLY, ALSO factor A
cross section from pp_meson in barns GeV^-1
factor= 1.0e-24* *1.0e-3 log(Ekin_factor)
*/

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

int gen_pi0_decay_emiss()   // generate pi0-decay emissivity
{
   cout<<"gen_pi0_decay_emiss"<<endl;
   cout<<"generating pi0-decay emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, NA1, NA2;

// identify the CR protons
   int iprotons=-1;
   for(int i=0; i<n_species; i++) if(gcr[i].A==1 && gcr[i].Z==1) iprotons=i;
   if(iprotons==-1) { cout<<"CR protons not found!"<<endl; return 1; }
   cout<<"  CR protons found as species #"<<iprotons<<endl;
 
// identify the CR Helium
   int iHelium =-1;
   for(int i=0; i<n_species; i++) if(gcr[i].A==4 && gcr[i].Z==2) iHelium =i;
   if(iHelium ==-1) { cout<<"CR Helium  not found!"<<endl; return 1; }
   cout<<"  CR Helium  found as species #"<<iHelium <<endl;

   int key1=0; // pi0-decay
   galaxy.pi0_decay_emiss=0.;
   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)  //pp_meson_cc( Esec, Pp1, NA1, NA2, key1 );
      {
         NA1=1;    NA2=1; //           beam+target     p+HI
         float cs_p_HI =pp_meson_cc(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
         NA1=1;    NA2=4; //           beam+target     p+He
         float cs_p_He =pp_meson_cc(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
         NA1=4;    NA2=1; //           beam+target     He+HI
         float cs_He_HI=pp_meson_cc(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1);
         NA1=4;    NA2=4; //           beam+target     He+He
         float cs_He_He=pp_meson_cc(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1);

         if(galaxy.n_spatial_dimensions==2)
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               { 
                  galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]+=       
                     (cs_p_HI +cs_p_He *galdef.He_H_ratio) *gcr[iprotons].cr_density.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip];
                  galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]+=       
                     (cs_He_HI+cs_He_He*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium ].A;
    
//cout<<"ir iz E_gamma bremss_emiss "<<ir<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "
//    <<galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]<<endl;  
               }//iz //ir //particle.n_spatial_dimensions==2
 
         if(galaxy.n_spatial_dimensions==3)
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {
                     galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]+=           
                        (cs_p_HI +cs_p_He *galdef.He_H_ratio) *gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip];
                     galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]+=            
                        (cs_He_HI+cs_He_He*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium ].A;                
             
//cout<<"ix iy  iz  E_gamma_pi0_decay_emiss "<<ix<<" "<<iy<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "
//    <<galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]<<endl; 
                  }//iz //iy //ix //particle.n_spatial_dimensions==3
         cout<<"gen_pi0_decay_emiss: iEgamma ip "<<iEgamma<<" "<<ip<<endl;
      } //ip //iEgamma
 
   double factor=1.0e-24*1.0e-3* log(galdef.Ekin_factor);
   galaxy.pi0_decay_emiss *= factor;

   if(galdef.verbose>=2) { cout<<"   pi0-decay emissivity "<<endl; galaxy.pi0_decay_emiss.print(); }
   cout<<" <<<< gen_pi0_decay_emiss"<<endl;

   return stat;
}

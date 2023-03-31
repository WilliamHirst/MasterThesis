
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propagate_particles.cc *                      galprop package * 5/03/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"

int propagate_particles()
{
   cout<<">>>>propagate_particles"<<endl;

   int i,net_iter;
   Particle particle;

   if(galdef.warm_start==0){                           //AWS20010121
     for(i=0; i<n_species; i++)gcr[i].cr_density=0.0;  //AWS20010121
   }                                                   //AWS20010121
   if(galdef.warm_start==1)read_gcr();                 //AWS20010121  

   particle.init();  // to signal that arrays not yet allocated

   particle=gcr[0];
   particle.create_transport_arrays();

   for (net_iter=1; net_iter<=galdef.network_iterations; net_iter++)
   {
      cout<<"    network iteration "<<net_iter<<endl;
//      for(i=0; i<n_species; i++)
      for(i=n_species-1; i>-1; i--)
      {
         particle=gcr[i];
         create_transport_arrays(particle);

// cout<<"particle.Dxx:"<<endl;particle.Dxx.print();

         cout<<">>>>generating secondary source for particle "<<i<<endl;
         if(gen_secondary_source(particle)!=0) return 1;
         if(galdef.verbose==10) particle.secondary_source_function.print();
         if(galdef.verbose>= 1) particle.print();
         if(galdef.verbose>= 0) cout<<"\n Network iteration "<<net_iter
            <<" species "<<i<<" "<<gcr[i].name<<" (Z,A) = ("<<gcr[i].Z
            <<","<<gcr[i].A<<")"<<endl; 

         if(propel(particle)!=0) return 1;
// cout<<"============particle.cr_density:"<<endl;particle.cr_density.print();

         gcr[i]=particle;
         if(galdef.verbose>= 1)
         {
            cout<<"\n Network iteration "<<net_iter
               <<" species "<<i<<"  "<<gcr[i].name<<endl;
            gcr[i].cr_density.print();
         }
      } 
//print_BC();
   } //nuc_iter

// convert density per momentum to flux per (KE/nucleon) for output NO done in store_gcr
// only for nuclei
// for(i=0; i<n_species; i++) if(gcr[i].A!=0) gcr[i].cr_density*=gcr[i].A;

   nuclei_normalize();
   if(galdef.primary_electrons==1) electrons_normalize();

   cr_luminosity();

   cout<<"<<<<propagate_particles"<<endl;
   return 0;
}



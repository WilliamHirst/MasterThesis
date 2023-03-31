
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * cr_luminosity.cc *                      galprop package * 10/03/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624


#include"galprop_classes.h"
#include"galproph.h"

// compute cosmic-ray luminosity of galaxy

/*
 

CR density  gcr.cr_density is in c/4pi * n(p) [ cm s^-1  * cm^-3 MeV^-1]
primary source function is    in c/4pi        [ cm s^-1  * cm^-3 MeV^-1 s-1]
luminosity (cm^-3)= Integral n(Ekin). A.Ekin. Ekin dlog(Ekin)
since each nucleus has KE=A.Ekin


The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
    UNITS OF c/4pi*n(p) = (1/A)flux(Ekin)
 =                (1/A)c/4pi*beta*n(Ekin)

so factor required = A.4pi/c/beta 
*/

int Galprop::cr_luminosity()
{
   cout<<"cr_luminosity"<<endl;
   cout<<"computing cr_luminosity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   int stat=0;
   double CR_luminosity=0;

// identify CR protons
   int iprotons=-1;
   for(int i=0;i<n_species;i++)if(gcr[i].A==1 && gcr[i].Z==1)iprotons=i;
   if(iprotons==-1){cout<<"CR protons not found!"<<endl;return 1;}
   cout<<"  CR protons found as species #"<<iprotons<<endl;
 
// identify CR Helium
   int iHelium =-1;
   for(int i=0;i<n_species;i++)if(gcr[i].A==4 && gcr[i].Z==2)iHelium =i;
   if(iHelium ==-1){cout<<"CR Helium  not found!"<<endl;return 1; }
   else cout<<"  CR Helium  found as species #"<<iHelium <<endl;

   Particle particle;
   particle.init();  // to signal that arrays not yet allocated

   particle=gcr[0];
   particle.create_transport_arrays();

   for (int iprHe = 1; iprHe<=2; iprHe++)
   {
     if(iprHe==1) particle=gcr[iprotons];
     if(iprHe==2) particle=gcr[iHelium ];
     create_transport_arrays(particle);
     for(int ip=0;ip<gcr[iprotons].n_pgrid;ip++)
     {
        if(galaxy.n_spatial_dimensions==2)
           for(int ir=0;ir<gcr[iprotons].n_rgrid;ir++)
              for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip] 
                    / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A  *particle.A 
                    *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
 
        if(galaxy.n_spatial_dimensions==3)
           for(int ix=0;ix<gcr[iprotons].n_xgrid;ix++)
              for(int iy=0;iy<gcr[iprotons].n_ygrid;iy++)
                 for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 {
// weights for non-symmetric case
                    double sym_weight=1.0;
// weights for fully symmetric case
                    if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
                    if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0
                    CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
                       / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A  *particle.A 
                       *galdef.dx*galdef.dy*galdef.dz *sym_weight;
                 }//iz//iy//ix//particle.n_spatial_dimensions==3
//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip
   }//iprHe

   CR_luminosity *= 4.0*Pi/c *1.0e6 *eV_to_erg *pow(kpc2cm,3) 
                    *log(gcr[iprotons].Ekin[1]/ gcr[iprotons].Ekin[0]) ;

//CR_luminosity *= gcr[iprotons].normalization_factor;//IMOS20030217

   cout<<endl;
   cout<<"====================================================="<<endl;
   cout<<" CR_luminosity="<<  CR_luminosity<<" erg s^-1"<<endl;
   cout<<"====================================================="<<endl;
   cout<<endl;

   cout<<" <<<< cr_luminosity"<<endl;
   return stat;
}

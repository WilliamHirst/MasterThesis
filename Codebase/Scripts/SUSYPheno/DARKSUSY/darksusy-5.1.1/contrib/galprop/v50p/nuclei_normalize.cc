
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nuclei_normalize.cc *                         galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// renormalization coefficient is defined from PRIMARY protons flux only;
// this renormalization coefficient then applied to all nuclei and secondary species

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galproph.h"

int Galprop::nuclei_normalize()
{
   cout<<">>>>nuclei_normalize"<<endl;

// identify the primary CR protons                            IMOS20000609
   int iprotons=-1;
   for(int i=n_species-1; i>=0; i--) if(gcr[i].Z==1&&gcr[i].A==1) {  iprotons=i;   break;  } // IMOS20010816
   if(iprotons==-1) { cout<<"CR protons not found!"<<endl; return 1; }
   if(strcmp(gcr[iprotons].name,"Hydrogen_1")!=0) { cout<<"CR protons not found!"<<endl; return 1; }
   cout<<"  CR protons found as species #"<<iprotons<<endl;

   double v1,v2,v3,v4,v5,v6;
   double r0=8.5; // solar Galactocentric radius, kpc

   int ip=(int)(log(galdef.proton_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 0.5);//IMOS20060420
   int iz=(int)((-galdef.z_min)/galdef.dz + 0.5);//IMOS20060420

   if(galdef.n_spatial_dimensions==2)
   {
      int ir=(int)((r0-galdef.r_min)/galdef.dr + 0.5);//IMOS20060420
      cout<<"Grid point for normalization: ir r[ir] iz z[iz] ip Ekin[ip] "<<ir<<" " <<gcr[iprotons].r[ir]
          <<" " <<iz <<" "<<gcr[iprotons].z[iz]<<" "<<ip<<" "<< gcr[iprotons].Ekin[ip]<<endl;

      v1=gcr[iprotons].cr_density.d2[ir  ][iz].s[ip];
      v2=gcr[iprotons].cr_density.d2[ir+1][iz].s[ip];
      v3=gcr[iprotons].cr_density.d2[ir  ][iz].s[ip+1];
      v4=gcr[iprotons].cr_density.d2[ir+1][iz].s[ip+1];
      v5=v1+(r0-gcr[iprotons].r[ir])/galdef.dr*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[iprotons].r[ir])/galdef.dr*(v4-v3); // r0 ip+1
   }

   if(galdef.n_spatial_dimensions==3)
   {
      int ix=(int)((r0-galdef.x_min)/galdef.dx + 0.5);//IMOS20060420
      int iy=(int)(( 0-galdef.y_min)/galdef.dy + 0.5);//IMOS20060420

      cout<<"Grid point for normalization: ix x[ix] iy y[iy] iz z[iz] ip Ekin[ip] "<<ix<<" " <<gcr[iprotons].x[ix]
          <<" "<<iy<<" "<< gcr[iprotons].y[iy]<<" "<<iz <<" "<<gcr[iprotons].z[iz]<<" "<<ip<<" "<< gcr[iprotons].Ekin[ip]<<endl;   //AWS20001121

      v1=gcr[iprotons].cr_density.d3[ix  ][iy][iz].s[ip];     //AWS20001121
      v2=gcr[iprotons].cr_density.d3[ix+1][iy][iz].s[ip];     //AWS20001121
      v3=gcr[iprotons].cr_density.d3[ix  ][iy][iz].s[ip+1];   //AWS20001121
      v4=gcr[iprotons].cr_density.d3[ix+1][iy][iz].s[ip+1];   //AWS20001121
      v5=v1+(r0-gcr[iprotons].x[ix])/galdef.dx*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[iprotons].x[ix])/galdef.dx*(v4-v3); // r0 ip+1
   }

   double vnorm=exp( log(v5)+log(galdef.proton_norm_Ekin/gcr[iprotons].Ekin[ip])/log(galdef.Ekin_factor)*log(v6/v5) );
   cout<<"v1 v2 v3 v4 v5 v6 vnorm  "<<v1<<" " <<v2<<" " <<v3 <<" "<<v4 <<" "<<v5<<" "<< v6<<" "<<vnorm <<endl;

   galdef.source_normalization *= galdef.proton_norm_flux/vnorm; // IMOS20030214
//   cout<<" nuclei_normalize >>> source_normalization= "<<galdef.source_normalization<<endl;

// normalize all species except primary electrons since these are normalized independently
   for(int i=0; i<n_species; i++) 
     {
       if(strcmp(gcr[i].name,"primary_electrons")==0) continue;// IMOS20050912

// don't re-normalize the DM annihilation products
       if(strcmp(gcr[i].name,"DM_positrons"  )==0) continue;  // IMOS20050912
       if(strcmp(gcr[i].name,"DM_electrons"  )==0) continue;  // IMOS20050912
       if(strcmp(gcr[i].name,"DM_antiprotons")==0) continue;  // IMOS20050912
       
       cout<<" normalizing "<<gcr[i].name<<endl;
       if(galdef.verbose>=1)
	 {
	   cout<<" nucleus "<<gcr[i].name<<" before normalization:"<<endl;
	   gcr[i].cr_density.print();
	 }
       
       gcr[i].cr_density         *=(galdef.proton_norm_flux/vnorm);
       gcr[i].normalization_factor=(galdef.proton_norm_flux/vnorm);//AWS20010121
       
       if(galdef.verbose>=1)
	 {
	   cout<<" nucleus "<<gcr[i].name<<" after normalization:"<<endl;
	   gcr[i].cr_density.print();
	 }
       
     }
   if(galdef.verbose>=1)
     {
       cout<<"primary protons after normalization:"<<endl;
       gcr[iprotons].cr_density.print();
     }  
   cout<<"<<<<nuclei_normalize"<<endl;
   return 0;
}

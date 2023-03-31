
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_isrf_energy_density.cc *                  galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>


#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"
#include"constants.h"

// generate ISRF energy density for all components




int gen_isrf_energy_density(){

  cout<<">>>>gen_isrf_energy_density"<<endl;
  

int stat=0;
int ir,ix,iy,iz;
float sum;

int i_comp;
 
galaxy.ISRF_energy_density=new Distribution[galaxy.n_ISRF_components]; 

 
/*
  ISRF is in Hz eV cm-3 Hz-1
  
  integral  energy density Hz-1  d(nu) = integral (nu* energy density Hz-1) d(log nu)
  d(log nu) is constant in this ISRF


factor=  LOG(nu(2)/nu(1)) 

*/


 


double factor=  log(galaxy.nu_ISRF[1]/ galaxy.nu_ISRF[0] ) ;

 if(galaxy.n_spatial_dimensions==2){

   cout<<"generating ISRF energy density for n_spatial_dimensions="<<galaxy.n_spatial_dimensions<<endl;

    for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
        galaxy.ISRF_energy_density[i_comp].init(galaxy.n_rgrid,galaxy.n_zgrid,1);

   for(ir=0;ir<galaxy.n_rgrid;ir++){
   for(iz=0;iz<galaxy.n_zgrid;iz++){
    
    for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++){

        
         sum=0.0;
         for(int inu=0;    inu<galaxy.ISRF[0].n_pgrid; inu++)
                 sum+=galaxy.ISRF[i_comp].d2[ir][iz].s[inu];

      
          galaxy.ISRF_energy_density[i_comp].d2[ir][iz].s[0]=sum*factor;
 

 

    }//ISRF_components

 
 }//iz
 }//ir
 
   }//particle.n_spatial_dimensions==2

 if(galaxy.n_spatial_dimensions==3){

   cout<<"generating ISRF energy density for n_spatial_dimensions="<<galaxy.n_spatial_dimensions<<endl;

    for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
        galaxy.ISRF_energy_density[i_comp].init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,1);

   for(ix=0;ix<galaxy.n_xgrid;ix++){
   for(iy=0;iy<galaxy.n_ygrid;iy++){
   for(iz=0;iz<galaxy.n_zgrid;iz++){
    
     // cout<<" "<<ix<<" "<<iy<<" "<<iz<<endl;

    for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++){  

         sum=0.0;
         for(int inu=0;    inu<galaxy.ISRF[0].n_pgrid;                inu++)
             sum+=galaxy.ISRF[i_comp].d3[ix][iy][iz].s[inu];

	 
          galaxy.ISRF_energy_density[i_comp].d3[ix][iy][iz].s[0]=sum*factor;
 

 

    }//ISRF_components

 
 }//iz
 }//iy
 }//ix
 
   }//particle.n_spatial_dimensions==3








if(galdef.verbose>=2)  
  for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++){
    cout<<"ISRF energy density for component #"<<i_comp<<endl;
galaxy.ISRF_energy_density[i_comp].print();}

cout<<" <<<< gen_isrf_energy_density"<<endl;
return stat;
}

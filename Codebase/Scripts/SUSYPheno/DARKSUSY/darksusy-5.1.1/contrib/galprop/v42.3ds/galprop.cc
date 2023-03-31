
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * galprop.cc *                                  galprop package * 5/04/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"

extern "C" int galprop_(char*); //IMOS20060221
int galprop(int,char**); //IMOS20060221

extern "C" int galprop_(char*header) //IMOS20060221
{
  //  char *a[]={"dstest",
  //	     "123456789012345678901234567890123456789012345678901234567890"};
  char *a0 = "dstest";
  char *a1 = &header[0];
  char *a[2] = { a0, a1 };
  galprop(2,a);
}

int galprop(int argc, char*argv[]) //IMOS20060221
{

//   skeleton(); for testing 
//   skeleton2(); for testing

  char version[]="42.3ds";  // IMOS20050912
  
  cout<<">>>>galprop version "<<version<<endl;
  
  if(argc!=2){cout<<"no galdef file specified!"<<endl;return -1;}
  cout<<argc<<" "<<argv[0]<<" "<<argv[1]<<endl;
  
  if(configure.init() !=0)return 1;
  cout<<configure.galdef_directory<<endl;// just a test
  
  if(galdef.read  (version,argv[1],configure.galdef_directory) !=0) return 1;
  
  if(galdef.test_suite>0) { test_suite(); return 0; }
  
//reading all isotopic cross sections, nuclear reaction network, nuclear fits
  read_nucdata();
  
//initialization of the Barashenkov & Polanski cross section code 
  sigtap_cc(-1);             // IMOS20010511 AWS20010620

//initialization of the Webber's routine
  set_sigma_cc();

  if(create_galaxy() !=0) return 1;
  if(create_gcr()    !=0) return 1;

//major routine
  if(propagate_particles() !=0) return 1;

//   if(print_BC() !=0) return 1;
  if(store_gcr() !=0) return 1;

  if(galdef.output_gcr_full!=0) if(store_gcr_full() !=0) return 1;

  if(galdef.synchrotron)
    {
     if(galdef.primary_electrons)
       {
         if (  gen_synch_emiss()  !=0)  return 1;
         if (  gen_synch_skymap() !=0)  return 1;
         if (store_synch_skymap() !=0)  return 1;
       } //galdef.primary_electrons
    }//galdef.synchrotron
  
  if(galdef.gamma_rays)
    {
      if(galdef.primary_electrons || galdef.secondary_electrons                       // IMOS20050912
	 || galdef.secondary_positrons || galdef.DM_positrons || galdef.DM_electrons) // IMOS20050912
	{
	  if (  gen_bremss_emiss()          !=0)  return 1;
	  if (store_bremss_emiss()          !=0)  return 1;
	  if (  gen_bremss_skymap        () !=0)  return 1;
	  if (store_bremss_skymap        () !=0)  return 1;
	  
	  if (  gen_bremss_ionized_skymap() !=0)  return 1;
	  if (store_bremss_ionized_skymap() !=0)  return 1;
	  
	  if (               gen_IC_emiss() !=0)  return 1;
	  if (               gen_IC_skymap()!=0)  return 1;
	  if (store_IC_skymap("isotropic")!=0)  return 1;
	  if(galdef.IC_anisotropic) if(store_IC_skymap("anisotropic")!=0) return 1; // AWS20010206
	  
	  if (store_IC_skymap_comp("isotropic")!=0)  return 1;
	  if(galdef.IC_anisotropic) if(store_IC_skymap_comp("anisotropic")!=0)  return 1;
	} //galdef.primary_electrons

      if (  gen_pi0_decay_emiss() !=0)  return 1;
      if (store_pi0_decay_emiss() !=0)  return 1;
      if (  gen_pi0_decay_skymap()!=0)  return 1;
      if (store_pi0_decay_skymap()!=0)  return 1;
   } //galdef.gamma_rays

  if(galdef.DM_gammas)                           // IMOS20050912
    {
      if (  gen_DM_emiss ()          !=0)  return 1;
      if (store_DM_emiss ()          !=0)  return 1;
      if (  gen_DM_skymap()          !=0)  return 1;
      if (store_DM_skymap()          !=0)  return 1;
    }

  if (  gen_ionization_rate()!=0)  return 1;
  if (store_ionization_rate()!=0)  return 1;
  
  cout<<"completed processing of galdef_"<<version<<"_"<<argv[1]<<endl;
  cout<<"<<<<galprop"<<endl;
  cout << '\a';                            // IMOS20010511
  
  return 0;
}

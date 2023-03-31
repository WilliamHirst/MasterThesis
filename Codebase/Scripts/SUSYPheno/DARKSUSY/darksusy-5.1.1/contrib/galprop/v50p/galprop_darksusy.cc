using namespace std;


#include"galprop_classes.h"
#include"galproph.h"

Galprop* galprop=0; // IMOS20060420

extern "C" int galprop_(char*header)
{
  char *a0 = "dstest";
  char *a1 = &header[0];
  char *a[2] = { a0, a1 };
  
  galprop = new Galprop;
  galprop->run(2, a);
  delete galprop;
  return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//IMOS20060420 transferred from fort_interface1.cc
// routine is called by FORTRAN routine emiss(r) (in file cfactor.f)
// hence underscore is appended and extern "C" supplied

extern "C" void isrf_energy_density_(float *r, float *z, float *energy_density)
{
   *energy_density=galprop->isrf_energy_density(*r,*z);
//   cout<<"energy_density = "<<*energy_density<<endl;
   return;
}

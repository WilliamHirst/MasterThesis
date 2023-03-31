using namespace std;


#include"galprop_classes.h"
#include"galproph.h"

Galprop* galprop=0; // IMOS20060420

main(int argc, char*argv[])
{
  galprop = new Galprop;
  galprop->run(argc, argv);
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

extern "C" void rho_darksusy__(double *x,double *y,double *z, double *rho)
{
  static int first = 1;
  *rho = -1.;
  if ( first ) {
    cerr<<"You have called rho_darksusy from GALPROP."<<endl<<
      "DarkSUSY is not present."<<endl<<"Did you really mean this?"<<endl;
    first = 0;
  }
  return;
}

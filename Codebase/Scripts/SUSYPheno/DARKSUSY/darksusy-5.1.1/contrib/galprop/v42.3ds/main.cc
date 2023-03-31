int galprop(int,char**);
extern "C" void rho_darksusy__(double*,double*,double*,double*);

int main(int argc, char*argv[])
{
  galprop(argc, argv);
}

extern "C" void rho_darksusy__(double *x,double *y,double *z, double *rho)
{
  *rho=-1.;
  return;
}

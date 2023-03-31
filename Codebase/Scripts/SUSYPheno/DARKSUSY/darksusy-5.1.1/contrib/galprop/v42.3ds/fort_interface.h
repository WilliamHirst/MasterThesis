
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * fort_interface.h *                            galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef fort_interface_h
#define fort_interface_h

// fort_interface1.cc prototypes

void set_sigma_cc();
double e_loss_compton_cc(double,double);
double bremss_spec_cc(double,double,int,int);
double pp_meson_cc(double,double,int,int,int);
int test_e_loss_compton_cc();
int test_bremss_spec_cc();
int test_pp_meson_cc();

// fort_interface2.cc prototypes

double wsigma_cc(int,int,int,int,double);                 // IMOS20020502
double yieldx_cc(int,int,int,int,float);                  // IMOS20020502
void sigtap_cc(int);                                      // IMOS20010511
double sighad_cc(int,double,double,double,double,double); // IMOS20020502
double cfactor_cc(int,double,double,double,double,double,double,double,double);
double antiproton_cc(int,double,double,int,int,int,int);  // IMOS20010511
double synchrotron_cc(double,double,double); 
int test_antiproton_cc();
int test_synchrotron_cc();
int test_cfactor();

#endif


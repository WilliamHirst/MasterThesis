
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * galprop.h *                                   galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// function prototypes

#ifndef galprop_h
#define galprop_h

#include "fort_interface.h"

void nucleon_cs(int,double,int,int,int,double*,double*,double*,double*,double*,double*);// IMOS20010511
void read_nucdata();
void Kcapture_cs(double,int,int,double*,double*);               // IMOS20010816
double isotope_cs(double,int,int,int,int,int,int*);
double nucdata(int,int,int,int,int,int,int*,int*,double*);      // IMOS20010816
double nucleon_loss(int,int,double,double,double,double,        // nucleon
          double*, double*);                                    // energy losses
double electron_loss(double,double,double,double,double,double, // electron
          double*,double*,double*,double*,double*,double*);     // energy losses

int create_gcr();
int create_transport_arrays(Particle&);
int create_galaxy();
int create_SNR();
int source_SNR_event(Particle&,double);
void decayed_cross_sections(int,int,int,int,double*,int,double*); // IMOS20010816
double sigma_boron_dec_heinbach_simon(int,int,int,int,double);
int gen_secondary_source(Particle&);
void read_nucdata();

// DM routines IMOS20050912
int gen_DM_source(Particle&);
int gen_DM_emiss();
double DM_profile(double, double, double);
double DM_profile_av(double,double,double,double,double);
double DM_profile_av(double,double,double,double,double,double,double);
int gen_DM_skymap();
int store_DM_emiss();
int store_DM_skymap();


int kinematic(int,int,char*,double&,double&,double&,double&,double&,double&,int);
int print_BC();
int propagate_particles();
int propel(Particle &particle);
int tridag    (float*,float*,float*,float*,float*,int);
int tridag_sym(float*,float*,float*,float*,float*,int);

void protri(Particle &particle,

Distribution  &alpha1_x,Distribution  &alpha1_y,Distribution  &alpha1_z,Distribution  &alpha1_p,
Distribution  &alpha2_x,Distribution  &alpha2_y,Distribution  &alpha2_z,Distribution  &alpha2_p,
Distribution  &alpha3_x,Distribution  &alpha3_y,Distribution  &alpha3_z,Distribution  &alpha3_p,

Distribution &Nx1_, Distribution &Ny1_, Distribution &Nz1_, Distribution &Np1_,
Distribution &Nx2_, Distribution &Ny2_, Distribution &Nz2_, Distribution &Np2_,
Distribution &Nx3_, Distribution &Ny3_, Distribution &Nz3_, Distribution &Np3_,

Distribution &total_source_function,

double dt, int nrept_outer ,double f_use);

int source_SNR_event_vec(Particle &particle,double t,
                         float *total_source_function_x);

int store_gcr();
int store_gcr_full();
int read_gcr();

int read_isrf();
int read_HIR();
int read_COR();


int gen_isrf_energy_density();
int He_to_H_CS(double,int,int,int,int,double*,double*);
double D_pp(double,double,double,double);

int e_KN_loss(Particle&);
double IC_cross_section(double,double,double);
double IC_anisotropy_factor(double,double,double,int);
double ionization_bethe(int Z, double beta);

double nHI (double,double);
double nH2 (double,double);
double nHII(double,double);

double nHI_av (double,double,double,double,double);
double nH2_av (double,double,double,double,double);
double nHII_av(double,double,double,double,double);


double source_distribution (double,double,double);
double        B_field_model(double,double,int);
double        B_field_model(double,double,double,int);

int propel_diagnostics
   (Particle &particle,
   Distribution &alpha1_r,
   Distribution &alpha1_z,
   Distribution &alpha1_p,
   Distribution &alpha2_r,
   Distribution &alpha2_z,
   Distribution &alpha2_p,
   Distribution &alpha3_r,
   Distribution &alpha3_z,
   Distribution &alpha3_p,
   Distribution &total_source_function,
   double dt
   );

int propel_diagnostics
   (Particle &particle,
   Distribution &alpha1_x,
   Distribution &alpha1_y,
   Distribution &alpha1_z,
   Distribution &alpha1_p,
   Distribution &alpha2_x,
   Distribution &alpha2_y,
   Distribution &alpha2_z,
   Distribution &alpha2_p,
   Distribution &alpha3_x,
   Distribution &alpha3_y,
   Distribution &alpha3_z,
   Distribution &alpha3_p,
   Distribution &total_source_function,
   double dt
   );

int electrons_normalize();
int nuclei_normalize();
int gen_IC_emiss();
int gen_IC_skymap();
int store_IC_skymap(char *IC_type);
int store_IC_skymap_comp(char *IC_type);
int gen_bremss_emiss();
int store_bremss_emiss();
int store_bremss_ionized_skymap();
int gen_bremss_ionized_skymap();
int gen_bremss_skymap();
int store_bremss_skymap();
int gen_pi0_decay_emiss();
int store_pi0_decay_emiss();
int gen_pi0_decay_skymap();
int store_pi0_decay_skymap();
int HIR_iRing(double RR); 

int gen_synch_emiss();
int gen_synch_skymap();
int store_synch_skymap();
int store_ionization_rate();

int gen_secondary_positron_source  (Particle &particle);
int gen_tertiary_antiproton_source (Particle &particle);   // IMOS20000605.15
int gen_secondary_antiproton_source(Particle &particle);
int gen_secondary_proton_source    (Particle &particle);   // IMOS20000605.16

int gen_ionization_rate();
int cr_luminosity();

float isrf_energy_density(float rr, float zz);

void skeleton();
void skeleton2();

double gauss(double mean, double sigma);

int test_sigma_boron_dec_heinbach_simon();
int test_kinematic();
int test_He_to_H_CS();
int test_nH();
int test_Distribution();
int test_Particle();
int test_isotope_cs();
int test_float_accuracy();
int test_source_SNR_event();
int test_suite();
int test_gauss(); //AWS20010410

#endif

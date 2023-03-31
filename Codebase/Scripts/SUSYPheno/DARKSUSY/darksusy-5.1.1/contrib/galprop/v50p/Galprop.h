
class Galprop
{
  
 public:
  int run(int argc, char*argv[]);
  
  int create_gcr();
  int create_transport_arrays(Particle&);
  int create_galaxy();
  int create_SNR();
  int cr_luminosity();
  int propagate_particles();
  int e_KN_loss(Particle &particle);
  
  double D_pp(double,double,double,double);
  int    D_xx (Particle&,int,int,int,int,int,int);
  static double fu(double);
  
// DM routines IMOS20050912
  int gen_DM_source(Particle&);
  int gen_DM_emiss();
  double DM_profile(double, double, double);
  double DM_profile_av(double,double,double,double,double);
  double DM_profile_av(double,double,double,double,double,double,double);
  int gen_DM_skymap();
  int store_DM_emiss();
  int store_DM_skymap();
  
  double IC_anisotropy_factor(double,double,double,double,int);//IMOS20060420
  void decayed_cross_sections(int iz,int ia,int jz,int ja, double *Ekin,int np,double *sigma);
  
  int nuclei_normalize();
  int electrons_normalize();
  
  int propel(Particle &particle);
  int propel_diagnostics(
			 Particle &particle, 
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
  int propel_diagnostics(
			 Particle &particle,
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
  
  void protri(Particle &particle,
	      Distribution  &alpha1_x,Distribution  &alpha1_y,Distribution  &alpha1_z,Distribution  &alpha1_p,
	      Distribution  &alpha2_x,Distribution  &alpha2_y,Distribution  &alpha2_z,Distribution  &alpha2_p,
	      Distribution  &alpha3_x,Distribution  &alpha3_y,Distribution  &alpha3_z,Distribution  &alpha3_p,
	      
	      Distribution &Nx1_, Distribution &Ny1_, Distribution &Nz1_, Distribution &Np1_,
	      Distribution &Nx2_, Distribution &Ny2_, Distribution &Nz2_, Distribution &Np2_,
	      Distribution &Nx3_, Distribution &Ny3_, Distribution &Nz3_, Distribution &Np3_,
	      
	      Distribution &total_source_function,
	      
	      double dt, int nrept_outer ,double f_use 
	      );
  double source_distribution (double x,double y,double z);
  int source_SNR_event(Particle &particle,double t);
  int source_SNR_event_vec(Particle &particle,double t,
			   float *total_source_function_x);
  
  int test_Particle();
  int test_source_SNR_event();
  int test_isotope_cs();
  int test_cfactor();
  
  int print_BC();
  float isrf_energy_density(float rr, float zz);
  int HIR_iRing(double RR);
  
  int read_gcr();
  int read_isrf();
  int read_HIR();
  int read_COR();
  
  int store_gcr();
  int store_gcr_full();
  int store_IC_skymap(char *IC_type);
  int store_IC_skymap_comp(char *IC_type);
  int store_bremss_emiss();
  int store_bremss_ionized_skymap();
  int store_bremss_skymap();
  int store_pi0_decay_emiss();
  int store_pi0_decay_skymap();
  int store_pi0_decay_HIR_skymap(); //AWS20041215
  int store_pi0_decay_H2R_skymap(); //AWS20041215
  int store_bremss_HIR_skymap();    //AWS20041215
  int store_bremss_H2R_skymap();    //AWS20041215
  int store_synch_skymap();
  int store_ionization_rate();
  
  int gen_secondary_source(Particle&);
  int gen_isrf_energy_density();
  int gen_IC_emiss();
  int gen_IC_skymap();
  int gen_bremss_emiss();
  int gen_bremss_ionized_skymap();
  int gen_bremss_skymap();
  int gen_pi0_decay_emiss();
  int gen_pi0_decay_skymap();
  int gen_synch_emiss();
  int gen_synch_skymap();
  int gen_secondary_positron_source  (Particle &particle);
  int gen_tertiary_antiproton_source (Particle &particle);   
  int gen_secondary_antiproton_source(Particle &particle);
  int gen_secondary_proton_source    (Particle &particle);   
  int gen_ionization_rate();
  
  int test_suite   ();
  
//////////////////////////////////////////////////
  Spectrum S;
  Distribution Dist;
  Particle P;
  
  int n_species;
  int isrf_energy_density_i_comp; // required for cfactor

///////////////////////////////////
  Configure configure;
  Galdef galdef;
  
  Particle *gcr; // all species
  Galaxy    galaxy;
};

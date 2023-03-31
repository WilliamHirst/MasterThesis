
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.h *                                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Galdef_h
#define Galdef_h

#include<stdio.h>
#include<string.h>

class Galdef
{
 public:

   char version[10];    // galprop version number
   char run_no [10];    // identifier of galprop run
   char galdef_ID[100]; // full identifier e.g. 01_123456

   int n_spatial_dimensions;             // 1,2 or 3
   double z_min,z_max,dz;                // for 1,2,3D    
   double r_min,r_max,dr;                // for 2D 
   double x_min,x_max,dx,y_min,y_max,dy; // for 3D 
   double p_min,p_max,p_factor;          // momentum start, end, factor
   double Ekin_min,Ekin_max,Ekin_factor; // kinetic energy per nucleon start, end, factor

   char p_Ekin_grid[5];                    // "p"||"Ekin": construct grid in p or Ekin

   double E_gamma_min,E_gamma_max,E_gamma_factor; // gamma-ray energy (MeV) start, end, factor
   double long_min,long_max;                      // gamma-ray skymap longitude min,max (degrees)
   double  lat_min, lat_max;                      // gamma-ray skymap latitude  min,max (degrees)
   double d_long,d_lat;                           // gamma-ray skymap longitude,latitude binsize (degrees)       

   double nu_synch_min,nu_synch_max,nu_synch_factor;// synchrotron frequency min,max (Hz), factor

   double D0_xx;                         // diffusion coefficient at break  rigidity   
   double D_rigid_br;                    // break     rigidity for diffusion coefficient in MV
   double D_g_1;                         // diffusion coefficient index below break  rigidity
   double D_g_2;                         // diffusion coefficient index above break  rigidity 

   int diff_reacc;                       // 1=include diffusive reacceleration
   double v_Alfven;                      // Alfven speed in km s-1

   int convection;                       // 1=include convection                       AWS20010323
   double   v0_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323
   double dvdz_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323

   double nuc_rigid_br;                  // break rigidity for nuclei injection  in MV
   double nuc_g_1;                       // spectral index below break  
   double nuc_g_2;                       // spectral index above break 

   char inj_spectrum_type[9];            // "rigidity"||"beta_rig"||"Etot": choose the nucleon injection spectrum IMOS20000613.1

   double electron_rigid_br;             // break rigidity for electron injection in MV
   double electron_g_1;                  // spectral index below break  
   double electron_g_2;                  // spectral index above break 

   double He_H_ratio;                    // He/H of ISM, by number
   double X_CO;                          // conversion factor CO integrated temperature -> H2 column density   AWS20010126
   int fragmentation;                    // 1=include fragmentation
   int momentum_losses;                  // 1=include momentum losses
   int radioactive_decay;                // 1=include radioactive_decay
   int K_capture;                        // 1=include K_capture                        AWS20010731

   double start_timestep;                // start time step in years
   double   end_timestep;                //   end time step in years
   double       timestep_factor;         //   factor to multiply timestep
   int          timestep_repeat;         //   number of times to repeat for factor
   int          timestep_repeat2;        //   number of  times to repeat in timestep_mode=2
   int          timestep_print ;         //   number of timesteps between printings
   int          timestep_diagnostics;    //   number of timesteps between propel diagnostics
   int           control_diagnostics;    //   control details of propel diagnostics

   int          network_iterations;      //   number of iterations of the entire network


   int prop_r;                           // for 2D: 1=propagate in r;
   int prop_x,prop_y,prop_z;             // for 2D: 1=propagate in z;for 3D 1=propagate in x,y,z
   int prop_p;                           // propagate in momentum
   int use_symmetry;                     // xyz symmetry (3D)

   int vectorized;                       // 0=unvectorized   code, 1=vectorized   code

   int source_specification;             // 0,1,2

   int  max_Z;                           // maximum   nucleus  Z in galdef file
   int *use_Z;                           // 1=use this nucleus Z

   double **isotopic_abundance;          // isotopic abundances

   int total_cross_section;              // total inel. cross-sections option AWS20010620
   int cross_section_option;             // controls which cross-sections used

   double t_half_limit;                  // lower limit on radioactive half-life for explicit inclusion  AWS20010731

   int primary_electrons;                // 1=propagate primary electrons
   int secondary_positrons;              // 1=propagate secondary positrons
   int secondary_electrons;              // 1=propagate secondary electrons
   int tertiary_antiprotons;             // 1=propagate tertiary antiprotons   IMOS20000605.13
   int secondary_antiprotons;            // 1=propagate secondary antiprotons
   int secondary_protons;                // 1=propagate secondary protons      IMOS20000605.14

   int gamma_rays;                       // 1=compute gamma-ray emission
   int IC_anisotropic;                   // 1=compute anisotropic inverse Compton
   int synchrotron;                      // 1=compute synchrotron emission


// DM parameters IMOS20050912
   int DM_positrons;                     // 1=compute positrons from DM
   int DM_electrons;                     // 1=compute electrons from DM
   int DM_antiprotons;                   // 1=compute antiprotons from DM
   int DM_gammas;                        // 1=compute gamma rays from DM
   double                                // user-defined params of DM (double)
     DM_double0, DM_double1, DM_double2, DM_double3, DM_double4,
     DM_double5, DM_double6, DM_double7, DM_double8, DM_double9;
   int                                   // user-defined params of DM (int)
     DM_int0, DM_int1, DM_int2, DM_int3, DM_int4,
     DM_int5, DM_int6, DM_int7, DM_int8, DM_int9;


   int source_model;                     //1= 2= 3= 5= 6=S&Mattox with cutoff
   double source_parameters[10];         // for source_model=1

   int   n_cr_sources;                   // number of pointlike cosmic-ray sources
   double *cr_source_x;                  // source x positions
   double *cr_source_y;                  // source y positions
   double *cr_source_z;                  // source z positions
   double *cr_source_w;                  // source width sigma in kpc
   double *cr_source_L;                  // source luminosity in TBD units

   int    SNR_events;                    // handle stochastic SNR events
   double SNR_interval;                  // time in years between SNRs in 1 kpc^3
   double SNR_livetime;                  // CR-producing live-time in years of an SNR
   double SNR_electron_sdg;              // delta electron source index Gaussian sigma
   double SNR_nuc_sdg;                   // delta nucleus  source index Gaussian sigma
   double SNR_electron_dgpivot;          // delta electron source index pivot rigidity (MeV)   AWS20010410
   double SNR_nuc_dgpivot;               // delta nucleus  source index pivot rigidity (MeV)   AWS20010410 

   int B_field_model;                    // >1000=parameterized model

   double   proton_norm_Ekin;            // proton   kinetic energy for normalization (MeV)
   double   proton_norm_flux;            // flux of protons   at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)
   double electron_norm_Ekin;            // electron kinetic energy for normalization (MeV)
   double electron_norm_flux;            // flux of electrons at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)

   int output_gcr_full;                  // output full 2D or 3D gcr array
   int warm_start;                       // read nuclei file and continue run   AWS20010121

   int verbose;                          // verbosity: 0=min 10=max 
   int test_suite;                       // run test suit instead of normal run 

//interface functions prototypes
   int read(char *version_,char  *run_no_,char *galdef_directory);
   int read_galdef_parameter(char *filename,char *parstring,double *value);
   int read_galdef_parameter(char *filename,char *parstring,int    *value);
   int read_galdef_parameter(char *filename,char *parstring,char *value);
   void print();
};

#endif











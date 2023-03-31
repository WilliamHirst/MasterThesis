
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include<string.h>
#include"Galdef.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::read(char *version_, char *run_no_, char *galdef_directory)
{
   int stat=0;
   cout<<" >>>>galdef read"<<endl;

   strcpy(version,version_);
   strcpy(run_no,run_no_);

   n_spatial_dimensions=2;

   r_min=1.0;
   r_max=15.0;
   dr  =2.0;

   x_min=-5.0;
   x_max=+6.0;
   dx  =2.0; 

   y_min=-7.0;
   y_max=+8.0;
   dy  =2.5; 

   z_min=-4.0;
   z_max=+4.0;
   dz  =1.5; 

   p_min=  10.0;
   p_max=100.0;
   p_factor=4. ;

   start_timestep=1.e5;
   end_timestep=1.e3;
   timestep_factor=0.5;
   timestep_repeat=20000;

   strcpy(galdef_ID,version);
   strcat(galdef_ID,"_");
   strcat(galdef_ID,run_no);

   char galdef_file[100];
   strcpy(galdef_file,galdef_directory);
   strcat(galdef_file,"galdef_");
   strcat(galdef_file,galdef_ID);
   cout<<galdef_file<<endl;


   char  filename[100];
   strcpy(filename,galdef_file);

   FILE *ft;
   ft=fopen(filename,"r");
   if(ft==NULL)
   {
      cout<<"no galdef file called "<<filename<<endl; return -1;
   }
   fclose(ft);

   char parstring[100];
  
   strcpy(parstring,                              "n_spatial_dimensions"  );
   stat= read_galdef_parameter(filename,parstring,&n_spatial_dimensions   );

   if(stat!=0)
   {
      cout<<"no spatial dimensions specified in "<<filename<<endl; return -1;
   }

//radial grid
   strcpy(parstring,                              "r_min"   );
   stat= read_galdef_parameter(filename,parstring,&r_min    );
   strcpy(parstring,                              "r_max"   );
   stat= read_galdef_parameter(filename,parstring,&r_max    );
   strcpy(parstring,                              "dr"      );
   stat= read_galdef_parameter(filename,parstring,&dr       );

//z-grid
   strcpy(parstring,                              "z_min"   );
   stat= read_galdef_parameter(filename,parstring,&z_min    );
   strcpy(parstring,                              "z_max"   );
   stat= read_galdef_parameter(filename,parstring,&z_max    );
   strcpy(parstring,                              "dz"      );
   stat= read_galdef_parameter(filename,parstring,&dz       );

//x-grid
   strcpy(parstring,                              "x_min"   );
   stat= read_galdef_parameter(filename,parstring,&x_min    );
   strcpy(parstring,                              "x_max"   );
   stat= read_galdef_parameter(filename,parstring,&x_max    );
   strcpy(parstring,                              "dx"      );
   stat= read_galdef_parameter(filename,parstring,&dx       );

//y-grid
   strcpy(parstring,                              "y_min"   );
   stat= read_galdef_parameter(filename,parstring,&y_min    );
   strcpy(parstring,                              "y_max"   );
   stat= read_galdef_parameter(filename,parstring,&y_max    );
   strcpy(parstring,                              "dy"      );
   stat= read_galdef_parameter(filename,parstring,&dy       );

//momentum grid
   strcpy(parstring,                              "p_min"   );
   stat= read_galdef_parameter(filename,parstring,&p_min    );
   strcpy(parstring,                              "p_max"   );
   stat= read_galdef_parameter(filename,parstring,&p_max    );
   strcpy(parstring,                              "p_factor");
   stat= read_galdef_parameter(filename,parstring,&p_factor );

//kinetic energy grid
   strcpy(parstring,                              "Ekin_min"   );
   stat= read_galdef_parameter(filename,parstring,&Ekin_min    );
   strcpy(parstring,                              "Ekin_max"   );
   stat= read_galdef_parameter(filename,parstring,&Ekin_max    );
   strcpy(parstring,                              "Ekin_factor");
   stat= read_galdef_parameter(filename,parstring,&Ekin_factor );

//p||Ekin option
   strcpy(parstring,                              "p_Ekin_grid");
   stat= read_galdef_parameter(filename,parstring, p_Ekin_grid );

//gamma-ray energy grid
   strcpy(parstring,                              "E_gamma_min"   );
   stat= read_galdef_parameter(filename,parstring,&E_gamma_min    );
   strcpy(parstring,                              "E_gamma_max"   );
   stat= read_galdef_parameter(filename,parstring,&E_gamma_max    );
   strcpy(parstring,                              "E_gamma_factor");
   stat= read_galdef_parameter(filename,parstring,&E_gamma_factor );

//synchrotron grid
   strcpy(parstring,                              "nu_synch_min"   );
   stat= read_galdef_parameter(filename,parstring,&nu_synch_min    );
   strcpy(parstring,                              "nu_synch_max"   );
   stat= read_galdef_parameter(filename,parstring,&nu_synch_max    );
   strcpy(parstring,                              "nu_synch_factor");
   stat= read_galdef_parameter(filename,parstring,&nu_synch_factor );

//longitude-latitude grid
   strcpy(parstring,                              "long_min"       );
   stat= read_galdef_parameter(filename,parstring,&long_min        );
   strcpy(parstring,                              "long_max"       );
   stat= read_galdef_parameter(filename,parstring,&long_max        );
   strcpy(parstring,                              "lat_min"        );
   stat= read_galdef_parameter(filename,parstring,&lat_min         );
   strcpy(parstring,                              "lat_max"        );
   stat= read_galdef_parameter(filename,parstring,&lat_max         );

   strcpy(parstring,                              "d_long"        );
   stat= read_galdef_parameter(filename,parstring,&d_long         );
   strcpy(parstring,                              "d_lat"         );
   stat= read_galdef_parameter(filename,parstring,&d_lat          );

//space diffusion
   strcpy(parstring,                              "D0_xx"      );
   stat= read_galdef_parameter(filename,parstring,&D0_xx       );

// diffusion coefficient parametrization
   strcpy(parstring,                              "D_rigid_br" );
   stat= read_galdef_parameter(filename,parstring,&D_rigid_br  );
   strcpy(parstring,                              "D_g_1"      );
   stat= read_galdef_parameter(filename,parstring,&D_g_1       );
   strcpy(parstring,                              "D_g_2"      );
   stat= read_galdef_parameter(filename,parstring,&D_g_2       );

//diffusive reacceleration: Alfven speed
   strcpy(parstring,                              "v_Alfven"      );
   stat= read_galdef_parameter(filename,parstring,&v_Alfven       );
   strcpy(parstring,                              "diff_reacc"      );
   stat= read_galdef_parameter(filename,parstring,&diff_reacc       );

// convection
   strcpy(parstring,                              "convection"      );  //AWS20010323
   stat= read_galdef_parameter(filename,parstring,&convection       );  //AWS20010323
   strcpy(parstring,                              "v0_conv"         );  //AWS20010323
   stat= read_galdef_parameter(filename,parstring,&v0_conv          );  //AWS20010323
   strcpy(parstring,                              "dvdz_conv"         );//AWS20010323
   stat= read_galdef_parameter(filename,parstring,&dvdz_conv          );//AWS20010323

//injection spectra parametrization (nucleons)
   strcpy(parstring,                              "nuc_rigid_br" );
   stat= read_galdef_parameter(filename,parstring,&nuc_rigid_br  );
   strcpy(parstring,                              "nuc_g_1"      );
   stat= read_galdef_parameter(filename,parstring,&nuc_g_1       );
   strcpy(parstring,                              "nuc_g_2"      );
   stat= read_galdef_parameter(filename,parstring,&nuc_g_2       );

//rigidity||beta_rig||Etot  option
   strcpy(parstring,                              "inj_spectrum_type"); // IMOS20000613.2
   stat= read_galdef_parameter(filename,parstring, inj_spectrum_type ); // IMOS20000613.3

//injection spectra parametrization (electrons)
   strcpy(parstring,                              "electron_rigid_br" );
   stat= read_galdef_parameter(filename,parstring,&electron_rigid_br  );
   strcpy(parstring,                              "electron_g_1"      );
   stat= read_galdef_parameter(filename,parstring,&electron_g_1       );
   strcpy(parstring,                              "electron_g_2"      );
   stat= read_galdef_parameter(filename,parstring,&electron_g_2       );

//other parameters
   strcpy(parstring,                              "He_H_ratio"   );
   stat= read_galdef_parameter(filename,parstring,&He_H_ratio    );

   strcpy(parstring,                              "X_CO"         );
   stat= read_galdef_parameter(filename,parstring,&X_CO          );

   strcpy(parstring,                              "fragmentation");
   stat= read_galdef_parameter(filename,parstring,&fragmentation );
   strcpy(parstring,                              "momentum_losses");
   stat= read_galdef_parameter(filename,parstring,&momentum_losses );
   strcpy(parstring,                              "radioactive_decay");
   stat= read_galdef_parameter(filename,parstring,&radioactive_decay );
   strcpy(parstring,                              "K_capture");          //AWS20010731
   stat= read_galdef_parameter(filename,parstring,&K_capture );          //AWS20010731


//time grid
   strcpy(parstring,                              "start_timestep");
   stat= read_galdef_parameter(filename,parstring,&start_timestep );
   strcpy(parstring,                              "end_timestep");
   stat= read_galdef_parameter(filename,parstring,&end_timestep );
   strcpy(parstring,                              "timestep_factor");
   stat= read_galdef_parameter(filename,parstring,&timestep_factor );
   strcpy(parstring,                              "timestep_repeat");
   stat= read_galdef_parameter(filename,parstring,&timestep_repeat );
   strcpy(parstring,                              "timestep_repeat2");
   stat= read_galdef_parameter(filename,parstring,&timestep_repeat2 );

   strcpy(parstring,                              "timestep_print" );
   stat= read_galdef_parameter(filename,parstring,&timestep_print  );
   strcpy(parstring,                              "timestep_diagnostics" );
   stat= read_galdef_parameter(filename,parstring,&timestep_diagnostics  );
   strcpy(parstring,                              "control_diagnostics" );
   stat= read_galdef_parameter(filename,parstring,&control_diagnostics  );

//other parameters
   strcpy(parstring,                              "network_iterations" );
   stat= read_galdef_parameter(filename,parstring,&network_iterations  );

   strcpy(parstring,                              "prop_r"        );
   stat= read_galdef_parameter(filename,parstring,&prop_r         );
   strcpy(parstring,                              "prop_x"        );
   stat= read_galdef_parameter(filename,parstring,&prop_x         );
   strcpy(parstring,                              "prop_y"        );
   stat= read_galdef_parameter(filename,parstring,&prop_y         );
   strcpy(parstring,                              "prop_z"        );
   stat= read_galdef_parameter(filename,parstring,&prop_z         );

   strcpy(parstring,                              "prop_p"        );
   stat= read_galdef_parameter(filename,parstring,&prop_p         );

   strcpy(parstring,                              "use_symmetry"  );
   stat= read_galdef_parameter(filename,parstring,&use_symmetry   );

   strcpy(parstring,                              "vectorized"    );
   stat= read_galdef_parameter(filename,parstring,&vectorized     );


   strcpy(parstring,                              "source_specification"        );
   stat= read_galdef_parameter(filename,parstring,&source_specification         );

//nuclei to include & abundances
   strcpy(parstring,                              "max_Z"         );
   stat= read_galdef_parameter(filename,parstring,&max_Z          );

   use_Z = new int[max_Z+1];
   for (int iZ=1; iZ<=max_Z; iZ++)
   {
      sprintf(parstring,"use_Z_%d",iZ);
      stat= read_galdef_parameter(filename,parstring,&use_Z[iZ]   );
   }

   isotopic_abundance = new double*[max_Z+1];
   for(int iZ=0; iZ<=max_Z; iZ++)
   {
      isotopic_abundance[iZ]=new double[max_Z*3];
      for(int iA=0; iA<max_Z*3; isotopic_abundance[iZ][iA++]=0.0); 
   }
   for (int iZ=1; iZ<=max_Z; iZ++)
      if(use_Z[iZ]==1)
         for (int iA=iZ; iA<iZ*3; iA++)
         {
            sprintf(parstring,"iso_abundance_%02d_%03d",iZ,iA); // eg _03_007
            cout<<parstring<<endl;
            stat= read_galdef_parameter(filename,parstring,&isotopic_abundance[iZ][iA]);
         }

   strcpy(parstring,                              "total_cross_section");//AWS20010620
   stat= read_galdef_parameter(filename,parstring,&total_cross_section );//AWS20010620


   strcpy(parstring,                              "cross_section_option");
   stat= read_galdef_parameter(filename,parstring,&cross_section_option );

   strcpy(parstring,                              "t_half_limit");       //AWS20010731
   stat= read_galdef_parameter(filename,parstring,&t_half_limit );       //AWS20010731

//CR species to include
   strcpy(parstring,                              "primary_electrons"    );
   stat= read_galdef_parameter(filename,parstring,&primary_electrons     );
   strcpy(parstring,                              "secondary_positrons"  );
   stat= read_galdef_parameter(filename,parstring,&secondary_positrons   );
   strcpy(parstring,                              "secondary_electrons"  );
   stat= read_galdef_parameter(filename,parstring,&secondary_electrons   );
   strcpy(parstring,                              "secondary_antiproton" );
   stat= read_galdef_parameter(filename,parstring,&secondary_antiprotons );
   strcpy(parstring,                              "tertiary_antiproton" );  // IMOS20000605.1
   stat= read_galdef_parameter(filename,parstring,&tertiary_antiprotons );  // IMOS20000605.2
   strcpy(parstring,                              "secondary_protons" );  // IMOS20000605.3
   stat= read_galdef_parameter(filename,parstring,&secondary_protons  );  // IMOS20000605.4

   strcpy(parstring,                              "gamma_rays"           );
   stat= read_galdef_parameter(filename,parstring,&gamma_rays            );

   strcpy(parstring,                              "IC_anisotropic"       );
   stat= read_galdef_parameter(filename,parstring,&IC_anisotropic        );

   strcpy(parstring,                              "synchrotron"         );
   stat= read_galdef_parameter(filename,parstring,&synchrotron          );

//CR source parameters
   strcpy(parstring,                              "source_model"       );
   stat= read_galdef_parameter(filename,parstring,&source_model        );
   strcpy(parstring,                              "source_parameters_1");
   stat= read_galdef_parameter(filename,parstring,&source_parameters[1]);
   strcpy(parstring,                              "source_parameters_2");
   stat= read_galdef_parameter(filename,parstring,&source_parameters[2]);
   strcpy(parstring,                              "source_parameters_3");
   stat= read_galdef_parameter(filename,parstring,&source_parameters[3]);

   strcpy(parstring,                              "n_cr_sources"       );
   stat= read_galdef_parameter(filename,parstring,&n_cr_sources        );
   cr_source_x         =new double[n_cr_sources];
   cr_source_y         =new double[n_cr_sources];
   cr_source_z         =new double[n_cr_sources];
   cr_source_w         =new double[n_cr_sources];
   cr_source_L         =new double[n_cr_sources];

   for(int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++)
   {
      sprintf(parstring,"cr_source_x_%02d", i_cr_source+1  ); // eg _01
      stat= read_galdef_parameter(filename,parstring,&cr_source_x[i_cr_source]  );

      sprintf(parstring,"cr_source_y_%02d", i_cr_source+1  ); 
      stat= read_galdef_parameter(filename,parstring,&cr_source_y[i_cr_source]  );

      sprintf(parstring,"cr_source_z_%02d", i_cr_source+1  ); 
      stat= read_galdef_parameter(filename,parstring,&cr_source_z[i_cr_source]  );

      sprintf(parstring,"cr_source_w_%02d", i_cr_source+1  ); 
      stat= read_galdef_parameter(filename,parstring,&cr_source_w[i_cr_source]  );

      sprintf(parstring,"cr_source_L_%02d", i_cr_source+1  ); 
      stat= read_galdef_parameter(filename,parstring,&cr_source_L[i_cr_source]  );
   }

//SNR
   strcpy(parstring,                              "SNR_events"          );
   stat= read_galdef_parameter(filename,parstring,&SNR_events           );
   strcpy(parstring,                              "SNR_interval"        );
   stat= read_galdef_parameter(filename,parstring,&SNR_interval         );
   strcpy(parstring,                              "SNR_livetime"        );
   stat= read_galdef_parameter(filename,parstring,&SNR_livetime         );

   strcpy(parstring,                              "SNR_electron_sdg"        );  //AWS20010410
   stat= read_galdef_parameter(filename,parstring,&SNR_electron_sdg         );  //AWS20010410
   strcpy(parstring,                              "SNR_nuc_sdg"             );  //AWS20010410
   stat= read_galdef_parameter(filename,parstring,&SNR_nuc_sdg              );  //AWS20010410

   strcpy(parstring,                              "SNR_electron_dgpivot"    );  //AWS20010410
   stat= read_galdef_parameter(filename,parstring,&SNR_electron_dgpivot     );  //AWS20010410
   strcpy(parstring,                              "SNR_nuc_dgpivot"         );  //AWS20010410
   stat= read_galdef_parameter(filename,parstring,&SNR_nuc_dgpivot          );  //AWS20010410


//magnetic field
   strcpy(parstring,                              "B_field_model"       );
   stat= read_galdef_parameter(filename,parstring,&B_field_model        );

//spectra normalization
   strcpy(parstring,                              "proton_norm_Ekin"       );
   stat= read_galdef_parameter(filename,parstring,&proton_norm_Ekin        );
   strcpy(parstring,                              "proton_norm_flux"       );
   stat= read_galdef_parameter(filename,parstring,&proton_norm_flux        );

   strcpy(parstring,                              "electron_norm_Ekin"       );
   stat= read_galdef_parameter(filename,parstring,&electron_norm_Ekin        );
   strcpy(parstring,                              "electron_norm_flux"       );
   stat= read_galdef_parameter(filename,parstring,&electron_norm_flux        );

//dark matter (DM) parameters  IMOS20050912
   strcpy(parstring,                              "DM_positrons"  );
   stat= read_galdef_parameter(filename,parstring,&DM_positrons   );
   strcpy(parstring,                              "DM_electrons"  );
   stat= read_galdef_parameter(filename,parstring,&DM_electrons   );
   strcpy(parstring,                              "DM_antiprotons");
   stat= read_galdef_parameter(filename,parstring,&DM_antiprotons );
   strcpy(parstring,                              "DM_gammas"     );
   stat= read_galdef_parameter(filename,parstring,&DM_gammas      );

   strcpy(parstring,                              "DM_double0");
   stat= read_galdef_parameter(filename,parstring,&DM_double0 );
   strcpy(parstring,                              "DM_double1");
   stat= read_galdef_parameter(filename,parstring,&DM_double1 );
   strcpy(parstring,                              "DM_double2");
   stat= read_galdef_parameter(filename,parstring,&DM_double2 );
   strcpy(parstring,                              "DM_double3");
   stat= read_galdef_parameter(filename,parstring,&DM_double3 );
   strcpy(parstring,                              "DM_double4");
   stat= read_galdef_parameter(filename,parstring,&DM_double4 );
   strcpy(parstring,                              "DM_double5");
   stat= read_galdef_parameter(filename,parstring,&DM_double5 );
   strcpy(parstring,                              "DM_double6");
   stat= read_galdef_parameter(filename,parstring,&DM_double6 );
   strcpy(parstring,                              "DM_double7");
   stat= read_galdef_parameter(filename,parstring,&DM_double7 );
   strcpy(parstring,                              "DM_double8");
   stat= read_galdef_parameter(filename,parstring,&DM_double8 );
   strcpy(parstring,                              "DM_double9");
   stat= read_galdef_parameter(filename,parstring,&DM_double9 );

   strcpy(parstring,                              "DM_int0");
   stat= read_galdef_parameter(filename,parstring,&DM_int0 );
   strcpy(parstring,                              "DM_int1");
   stat= read_galdef_parameter(filename,parstring,&DM_int1 );
   strcpy(parstring,                              "DM_int2");
   stat= read_galdef_parameter(filename,parstring,&DM_int2 );
   strcpy(parstring,                              "DM_int3");
   stat= read_galdef_parameter(filename,parstring,&DM_int3 );
   strcpy(parstring,                              "DM_int4");
   stat= read_galdef_parameter(filename,parstring,&DM_int4 );
   strcpy(parstring,                              "DM_int5");
   stat= read_galdef_parameter(filename,parstring,&DM_int5 );
   strcpy(parstring,                              "DM_int6");
   stat= read_galdef_parameter(filename,parstring,&DM_int6 );
   strcpy(parstring,                              "DM_int7");
   stat= read_galdef_parameter(filename,parstring,&DM_int7 );
   strcpy(parstring,                              "DM_int8");
   stat= read_galdef_parameter(filename,parstring,&DM_int8 );
   strcpy(parstring,                              "DM_int9");
   stat= read_galdef_parameter(filename,parstring,&DM_int9 );


//output controls
   strcpy(parstring,                              "output_gcr_full"       );
   stat= read_galdef_parameter(filename,parstring,&output_gcr_full       );

  
   strcpy(parstring,                              "warm_start"       );//AWS20010121
   stat= read_galdef_parameter(filename,parstring,&warm_start        );//AWS20010121

   strcpy(parstring,                              "verbose"       );
   stat= read_galdef_parameter(filename,parstring,&verbose        );
   strcpy(parstring,                              "test_suite"       );
   stat= read_galdef_parameter(filename,parstring,&test_suite        );

   print();
   cout<<" <<<<galdef read"<<endl;
   return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::read_galdef_parameter(char *filename, char *parstring, double *value)
{

   FILE *ft;
   ft=fopen(filename,"r"); if (ft==NULL) return -1;
   char *flag;
   int   found;   
  
   char dumstring[100];
   char input[200];
   int stat=0;

   found=-1;
   while(found!=0)
   {
      flag=fgets(input,1000,ft); // read string until newline (Schildt p.222)
      if(flag==0)
      {
         cout<<parstring<<" not found in galdef file!"<<endl;
         fclose(ft); 
         stat=1;
         return stat;
      }
//  printf("%s",input);       // string is \0 terminated
      sscanf(input,"%s"  ,dumstring  );
      found=strcmp(dumstring,parstring); // search for parstring  in input
  
      if(found==0)
      {
         sscanf(input,"%22c%le"  ,dumstring, value);  // le for double
//  cout<<parstring<<"    found!"<<endl;
//  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
      }
   }

   fclose(ft); 

//  cout<<" <<<< read_galdef_parameter"<<endl;
   return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::read_galdef_parameter(char *filename, char *parstring, int *value)
{
   FILE *ft;
   ft=fopen(filename,"r"); if (ft==NULL) return -1;
   char *flag;
   int   found; 
   char dumstring[100];
   char input[200];
   int stat=0;

   found=-1;
   while(found!=0)
   {
      flag=fgets(input,1000,ft);   // read string until newline (Schildt p.222)
      if(flag==0)
      {
         cout<<parstring<<" not found in galdef file!"<<endl;
         fclose(ft);
         stat=1;
         return stat;
      }
//  printf("%s",input);       // string is \0 terminated
      sscanf(input,"%s"  ,dumstring  );
      found=strcmp(dumstring,parstring); // search for parstring  in input
  
      if(found==0)
      {
         sscanf(input,"%22c%d"   ,dumstring, value);  // d for int    
//  cout<<parstring<<"    found!"<<endl;
//  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
      }  
   }
   fclose(ft); 
//  cout<<" <<<< read_galdef_parameter"<<endl;
   return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::read_galdef_parameter(char *filename, char *parstring, char *value)
{
   FILE *ft;
   ft=fopen(filename,"r"); if (ft==NULL) return -1;
   char *flag;
   int   found;   
   char dumstring[100];
   char input[200];
   int stat=0;

   found=-1; 
   while(found!=0)
   {
      flag=fgets(input,1000,ft);   // read string until newline (Schildt p.222)
      if(flag==0)
      {
         cout<<parstring<<" not found in galdef file!"<<endl;
         fclose(ft); 
         stat=1;
         return stat;
      }
//  printf("%s",input);       // string is \0 terminated
      sscanf(input,"%s"  ,dumstring  );
      found=strcmp(dumstring,parstring); // search for parstring  in input
   
      if(found==0)
      {
	 sscanf(input,"%22c%s"   ,dumstring, value      ); // s for string 
// cout<<parstring<<"    found!"<<endl;
// cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
      }  
   }
   fclose(ft); 
//  cout<<" <<<< read_galdef_parameter"<<endl;
   return stat;
}

void Galdef::print()
{
   cout<<"  ======= galdef: "<<galdef_ID  <<endl;
   cout<<"  version  "<<version <<endl;
   cout<<"  run_no   "<<run_no      <<endl;

   cout<<"  n_spatial_dimensions "<<n_spatial_dimensions   <<endl;

   cout<<"  r_min    "<<r_min   <<endl;
   cout<<"  r_max    "<<r_max   <<endl; 
   cout<<"  dr       "<<dr      <<endl;
   cout<<"  z_min    "<<z_min   <<endl;
   cout<<"  z_max    "<<z_max   <<endl;
   cout<<"  dz       "<<dz      <<endl;
   cout<<"  p_min    "<<p_min   <<endl;
   cout<<"  p_max    "<<p_max   <<endl;
   cout<<"  p_factor "<<p_factor<<endl;
   cout<<"  Ekin_min "<<   Ekin_min      <<endl;
   cout<<"  Ekin_max "<<   Ekin_max      <<endl;
   cout<<"  Ekin_factor "<<Ekin_factor   <<endl;

   cout<<"  p_Ekin_grid "<<p_Ekin_grid   <<endl;

   cout<<"  E_gamma_min   "<<   E_gamma_min      <<endl;
   cout<<"  E_gamma_max   "<<   E_gamma_max      <<endl;
   cout<<"  E_gamma_factor"<<   E_gamma_factor   <<endl;

   cout<<" nu_synch_min    "<<  nu_synch_min      <<endl;
   cout<<" nu_synch_max    "<<  nu_synch_max      <<endl;
   cout<<" nu_synch_factor "<<  nu_synch_factor   <<endl;

   cout<<"  long_min      "<<   long_min         <<endl;
   cout<<"  long_max      "<<   long_max         <<endl;
   cout<<"   lat_min      "<<    lat_min         <<endl;
   cout<<"   lat_max      "<<    lat_max         <<endl;
   cout<<"  d_long        "<<   d_long           <<endl;
   cout<<"  d_lat         "<<   d_lat            <<endl;

   cout<<"  D0_xx       "<<D0_xx         <<endl;
   cout<<"  D_rigid_br  "<<D_rigid_br    <<endl;
   cout<<"  D_g_1       "<<D_g_1         <<endl;
   cout<<"  D_g_2       "<<D_g_2         <<endl;
   cout<<"  diff_reacc  "<<diff_reacc    <<endl;
   cout<<"  v_Alfven    "<<v_Alfven      <<endl;
   cout<<"  convection  "<<convection    <<endl; //AWS20010323
   cout<<"  v0_conv     "<<v0_conv       <<endl; //AWS20010323
   cout<<"  dvdz_conv   "<<dvdz_conv     <<endl; //AWS20010323


   cout<<"  nuc_rigid_br  "<<nuc_rigid_br    <<endl;
   cout<<"  nuc_g_1       "<<nuc_g_1         <<endl;
   cout<<"  nuc_g_2       "<<nuc_g_2         <<endl;

   cout<<"  inj_spectrum_type "<<inj_spectrum_type<<endl; // IMOS20000613.4

   cout<<"  electron_rigid_br  "<<electron_rigid_br    <<endl;
   cout<<"  electron_g_1       "<<electron_g_1         <<endl;
   cout<<"  electron_g_2       "<<electron_g_2         <<endl;

   cout<<"  He_H_ratio        "<<He_H_ratio          <<endl;
   cout<<"  X_CO              "<<X_CO                <<endl; //AWS20010126
   cout<<"  fragmentation     "<<fragmentation       <<endl;
   cout<<"  momentum_losses   "<<momentum_losses     <<endl;
   cout<<"  radioactive_decay "<<radioactive_decay   <<endl;
   cout<<"  K_capture         "<<K_capture           <<endl;

   cout<<"  x_min    "<<x_min   <<endl;
   cout<<"  x_max    "<<x_max   <<endl;
   cout<<"  dx       "<<dx      <<endl;
   cout<<"  y_min    "<<y_min   <<endl;
   cout<<"  y_max    "<<y_max   <<endl;
   cout<<"  dy       "<<dy      <<endl;

   cout<<"  start_timestep   "<< start_timestep<<endl;
   cout<<"    end_timestep   "<<   end_timestep<<endl;
   cout<<"  timestep_factor  "<<timestep_factor<<endl;
   cout<<"  timestep_repeat  "<<timestep_repeat<<endl;
   cout<<"  timestep_repeat2 "<<timestep_repeat2<<endl;
   cout<<"  timestep_print   "<<timestep_print <<endl;
   cout<<"  timestep_diagnostics  "<<timestep_diagnostics <<endl;
   cout<<"  control_diagnostics  " << control_diagnostics <<endl;

   cout<<"  network_iterations "<<network_iterations<<endl;

   cout<<"  prop_r   "<<prop_r  <<endl;
   cout<<"  prop_x   "<<prop_x  <<endl;
   cout<<"  prop_y   "<<prop_y  <<endl;
   cout<<"  prop_z   "<<prop_z  <<endl;
   cout<<"  prop_p   "<<prop_p  <<endl;

   cout<<"  use_symmetry  "<<use_symmetry<<endl;
   cout<<"  vectorized    "<<vectorized  <<endl;

   cout<<"  source_specification  "<<source_specification <<endl;
   cout<<"  source_model          "<<source_model         <<endl;
   cout<<"  source_parameters_1   "<<source_parameters[1] <<endl;
   cout<<"  source_parameters_2   "<<source_parameters[2] <<endl;
   cout<<"  source_parameters_3   "<<source_parameters[3] <<endl;

   cout<<"  n_cr_sources          "<<n_cr_sources         <<endl;
   for(int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++)
   {
      cout<<"  cr_source_x["<<i_cr_source<<"]  "<<   cr_source_x[i_cr_source]<<endl;
      cout<<"  cr_source_y["<<i_cr_source<<"]  "<<   cr_source_y[i_cr_source]<<endl;
      cout<<"  cr_source_z["<<i_cr_source<<"]  "<<   cr_source_z[i_cr_source]<<endl;
      cout<<"  cr_source_w["<<i_cr_source<<"]  "<<   cr_source_w[i_cr_source]<<endl;
      cout<<"  cr_source_L["<<i_cr_source<<"]  "<<   cr_source_L[i_cr_source]<<endl;
   }

   cout<<"  SNR_events            "<<SNR_events           <<endl;
   cout<<"  SNR_interval          "<<SNR_interval         <<endl;
   cout<<"  SNR_livetime          "<<SNR_livetime         <<endl;

   cout<<"  SNR_electron_sdg      "<<SNR_electron_sdg     <<endl;
   cout<<"  SNR_nuc_sdg           "<<SNR_nuc_sdg          <<endl;
   cout<<"  SNR_electron_dgpivot  "<<SNR_electron_dgpivot <<endl;
   cout<<"  SNR_nuc_dgpivot       "<<SNR_nuc_dgpivot      <<endl;


   cout<<"  B_field_model         "<<B_field_model        <<endl;

   cout<<"    proton_norm_Ekin    "<<  proton_norm_Ekin   <<endl;
   cout<<"    proton_norm_flux    "<<  proton_norm_flux   <<endl;
   cout<<"  electron_norm_Ekin    "<<electron_norm_Ekin   <<endl;
   cout<<"  electron_norm_flux    "<<electron_norm_flux   <<endl;

   cout<<"  max_Z                 "<<max_Z                <<endl;
   for(int iZ=1; iZ<=max_Z; iZ++) cout<<"use_Z_"<<iZ<<" "<<use_Z[iZ]<<endl;

   for (int iZ=1; iZ<=max_Z; iZ++)
      for (int iA=iZ; iA<max_Z*3; iA++)
      {
         if(isotopic_abundance[iZ][iA]!=0.0)
            cout<<"isotopic abundance for (Z,A) ("
               <<iZ<<","<<iA<<") ="<<isotopic_abundance[iZ][iA]<<endl ;
      }

   cout<<"  total_cross_section  " <<total_cross_section  <<endl;  //AWS20010620

   cout<<"  cross_section_option  "<<cross_section_option <<endl;
   cout<<"  t_half_limit          "<<t_half_limit         <<endl;  //AWS20010731


   cout<<"  primary_electrons     "<<primary_electrons    <<endl;
   cout<<"  secondary_positrons   "<<secondary_positrons  <<endl;
   cout<<"  secondary_electrons   "<<secondary_electrons  <<endl;
   cout<<"  secondary_antiprotons "<<secondary_antiprotons<<endl;
   cout<<"  tertiary_antiprotons  "<<tertiary_antiprotons <<endl;  // IMOS20000605.5
   cout<<"  secondary_protons     "<<secondary_protons    <<endl;  // IMOS20000605.6
   cout<<"  gamma_rays            "<<gamma_rays           <<endl;
   cout<<"  IC_anisotropic        "<<IC_anisotropic       <<endl;
   cout<<"  synchrotron           "<<synchrotron          <<endl;

// DM: IMOS20050912
   cout<<"  DM_positrons          "<<DM_positrons         <<endl;
   cout<<"  DM_electrons          "<<DM_electrons         <<endl;
   cout<<"  DM_antiprotons        "<<DM_antiprotons       <<endl;
   cout<<"  DM_gammas             "<<DM_gammas            <<endl;

   cout<<"  DM_double0            "<<DM_double0           <<endl;
   cout<<"  DM_double1            "<<DM_double1           <<endl;
   cout<<"  DM_double2            "<<DM_double2           <<endl;
   cout<<"  DM_double3            "<<DM_double3           <<endl;
   cout<<"  DM_double4            "<<DM_double4           <<endl;
   cout<<"  DM_double5            "<<DM_double5           <<endl;
   cout<<"  DM_double6            "<<DM_double6           <<endl;
   cout<<"  DM_double7            "<<DM_double7           <<endl;
   cout<<"  DM_double8            "<<DM_double8           <<endl;
   cout<<"  DM_double9            "<<DM_double9           <<endl;

   cout<<"  DM_int0               "<<DM_int0              <<endl;
   cout<<"  DM_int1               "<<DM_int1              <<endl;
   cout<<"  DM_int2               "<<DM_int2              <<endl;
   cout<<"  DM_int3               "<<DM_int3              <<endl;
   cout<<"  DM_int4               "<<DM_int4              <<endl;
   cout<<"  DM_int5               "<<DM_int5              <<endl;
   cout<<"  DM_int6               "<<DM_int6              <<endl;
   cout<<"  DM_int7               "<<DM_int7              <<endl;
   cout<<"  DM_int8               "<<DM_int8              <<endl;
   cout<<"  DM_int9               "<<DM_int9              <<endl;



   cout<<"  output_gcr_full      "<<output_gcr_full       <<endl;
   cout<<"  warm_start           "<<warm_start            <<endl;

   cout<<"  verbose    "<<verbose    <<endl;
   cout<<"  test_suite "<<test_suite <<endl;
}













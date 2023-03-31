
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_IC_skymap_comp.cc *                     galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"
#include "fitsio.h"

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// Stores IC skymap by components
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::store_IC_skymap_comp(char* IC_type)
{
  int stat=0;
  cout<<" >>>> store_IC_skymap_comp"<<endl;
  fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
  int status, ii, jj;
  long  fpixel = 1, naxis = 4, nelements, exposure;
  long naxes[4];
  
  naxes[0]=galaxy.n_long;
  naxes[1]=galaxy.n_lat;             
  naxes[2]=galaxy.n_E_gammagrid;
  naxes[3]=galaxy.n_ISRF_components;
  nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

  float *array;
  array=new float[nelements];
  char outfile[100];
  strcpy(outfile,"!"); // create new file or overwrite existing one
  strcat(outfile,configure.fits_directory);
  
  if(strcmp(IC_type,"isotropic")!=0 && strcmp(IC_type,"anisotropic")!=0)
    { cout<<" store_IC_skymap_comp: invalid IC_type "<<IC_type<<endl; return 1; }
  if(strcmp(IC_type,  "isotropic")==0) strcat(outfile,"ics_skymap_comp_");
  if(strcmp(IC_type,"anisotropic")==0) strcat(outfile,"ics_anisotropic_skymap_comp_");
  
  strcat(outfile,galdef.galdef_ID);
  cout<<" storing IC "<<IC_type<< " skymap components in file "<<outfile<<endl;  

  status = 0;         // initialize status before calling fitsio routines
  fits_create_file(&fptr, outfile, &status);   // create new file or overwrite existing one
//cout<<"create_file"<<endl;
  fits_report_error(stderr, status);  // print out any error messages
  
// Create the primary array image (32-bit float pixels)
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
//cout<<"create_img"<<endl;
  fits_report_error(stderr, status);  // print out any error messages

// Write a keyword; must pass the ADDRESS of the value
  exposure = 1500;
  fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);
//cout<<"update_key"<<endl;
  fits_report_error(stderr, status);  // print out any error messages
    
  int i=0;
  for (int i_comp=0; i_comp<naxes[3]; i_comp++)
    for (int ip=0; ip<naxes[2]; ip++)
      for (int ib=0; ib<naxes[1]; ib++)
	for (int il=0; il<naxes[0]; il++)
	  {
	    if( strcmp(IC_type,  "isotropic")==0) array[i]=galaxy.IC_iso_skymap  [i_comp].d2[il][ib].s[ip];
	    if( strcmp(IC_type,"anisotropic")==0) array[i]=galaxy.IC_aniso_skymap[i_comp].d2[il][ib].s[ip];
	    array[i]*=pow(galaxy.E_gamma[ip],2);
	    i++;
	  }
// Write the array of floats to the image
  fits_write_img(fptr, TFLOAT, fpixel, nelements, array, &status);
//cout<<"write_img"<<endl;
  fits_report_error(stderr, status);  // print out any error messages

// Write basic FITS keywords
  double crval1,crval2,crval3,crval4;
  double cdelt1,cdelt2,cdelt3,cdelt4;
  
  crval1=galaxy.long_min;
  crval2=galaxy. lat_min;
  crval3=log10(galaxy.E_gamma_min);
  crval4=1;
  
  cdelt1=galaxy.d_long;
  cdelt2=galaxy.d_lat;
  cdelt3=log10(galaxy.E_gamma_factor);
  cdelt4=1;
  
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);

//cout<<"update_key"<<endl;
  fits_report_error(stderr, status);  // print out any error messages
  fits_close_file(fptr, &status);     // close the file
  fits_report_error(stderr, status);  // print out any error messages
  return( status );
  cout<<" <<<< store_IC_skymap_comp"<<endl;
  return stat;
}

/*
    // write keywords describing nuclei
char keyword[20];
char comment[40];

 for(int i_nucleus=0;i_nucleus<n_species;i_nucleus++){
   sprintf(keyword,"NUCZ%03d",       i_nucleus+1 ); // e.g. NUCZ012
   sprintf(comment,"Z of nucleus %d",i_nucleus+1 ); 
   cout<<keyword<<" "<<gcr[i_nucleus].Z<<endl;
   fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].Z,comment, &status);

   sprintf(keyword,"NUCA%03d",       i_nucleus+1 ); // e.g. NUCA012
   sprintf(comment,"A of nucleus %d",i_nucleus+1 ); 
      cout<<keyword<<" "<<gcr[i_nucleus].A<<endl;
   fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].A,comment, &status);


 }
*/


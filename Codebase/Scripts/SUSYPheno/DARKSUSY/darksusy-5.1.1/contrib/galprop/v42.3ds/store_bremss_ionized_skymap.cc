
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_ionized_skymap.cc *              galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"

#include "fitsio.h" 
int store_bremss_ionized_skymap(){
int stat;

  cout<<" >>>> store_bremss_ionized_skymap"<<endl;
stat=0;
   fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, ii, jj;
    long  fpixel = 1, naxis = 4, nelements, exposure;
    long naxes[4]  ; 

  
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;             
    naxes[2]=galaxy.n_E_gammagrid;
    naxes[3]=1;

    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
    float *array;          
    array=new float[nelements];

    char outfile[100];
    strcpy(outfile,"!"); /* create new file or overwrite existing one */
    strcat(outfile,configure.fits_directory);
    strcat(outfile,"bremss_ionized_skymap_");
    strcat(outfile,galdef.galdef_ID);
    cout<<" storing nuclei in file "<<outfile<<endl;

    status = 0;         /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, outfile, &status);   /* create new file or overwrite existing one */

    /* Create the primary array image (32-bit float pixels */
    fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

    /* Write a keyword; must pass the ADDRESS of the value */
    exposure = 1500;
    fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status);

    
 

  
   int i=0;
  
   for (int ip       =0;        ip<naxes[2];       ip++){
   for (int ib       =0;        ib<naxes[1];       ib++){
   for (int il       =0;        il<naxes[0];       il++){
     array[i]=galaxy.bremss_ionized_skymap.d2[il][ib].s[ip];i++;
   }
   }
   }
   
    

 

    /* Write the array of floats to the image */
    fits_write_img(fptr, TFLOAT, fpixel, nelements, array, &status);


    // write basic FITS keywords

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

    fits_close_file(fptr, &status);            /* close the file */

    fits_report_error(stderr, status);  /* print out any error messages */
    return( status );


  cout<<" <<<< store_bremss_ionized_skymap"<<endl;
return stat;
}

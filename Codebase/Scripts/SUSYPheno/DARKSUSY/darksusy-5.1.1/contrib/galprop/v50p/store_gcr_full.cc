
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_gcr_full.cc *                           galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

#include "fitsio.h" 
int Galprop::store_gcr_full()
{
int stat;

  cout<<" >>>> store_gcr_full"<<endl;
stat=0;
   fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, ii, jj;
    long  fpixel = 1, naxis = 4, nelements, exposure;
    long naxes[5]  ; 

    if(gcr[0].n_spatial_dimensions==2){
     naxis=4;
     naxes[0]=gcr[0].n_rgrid;
     naxes[1]=gcr[0].n_zgrid;
     naxes[2]=gcr[0].n_pgrid;
     naxes[3]=n_species;
     nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
    }

    if(gcr[0].n_spatial_dimensions==3){
     naxis=5;
     naxes[0]=gcr[0].n_xgrid;
     naxes[1]=gcr[0].n_ygrid;
     naxes[2]=gcr[0].n_zgrid;
     naxes[3]=gcr[0].n_pgrid;
     naxes[4]=n_species;
     nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3]*naxes[4];
    }

    
    float *array;          
    array=new float[nelements];

    char outfile[100];
    strcpy(outfile,"!"); /* create new file or overwrite existing one */
    strcat(outfile,configure.fits_directory);
    strcat(outfile,"nuclei_full_");
    strcat(outfile,galdef.galdef_ID);
    cout<<" storing full nuclei array in file "<<outfile<<endl;

    status = 0;         /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, outfile, &status);   /* create new file or overwrite existing one */

    /* Create the primary array image (32-bit float pixels */
    fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

    /* Write a keyword; must pass the ADDRESS of the value */
    exposure = 1500;
    fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status);

    
  

    if(gcr[0].n_spatial_dimensions==2){
   int i=0;
   for (int i_species=0;i_species< naxes[3];i_species++){ 
   for (int ip       =0;        ip<naxes[2];       ip++){
   for (int iz       =0;        iz<naxes[1];       iz++){
   for (int ir       =0;        ir<naxes[0];       ir++){

     if(gcr[i_species].A!=0)
     array[i]=    gcr[i_species].cr_density.d2[ir]    [iz].s[ip] 
             *    gcr[i_species].A
             *pow(gcr[i_species].Ekin[ip],2);

     if(gcr[i_species].A==0)  // electrons, positrons
     array[i]=    gcr[i_species].cr_density.d2[ir]    [iz].s[ip] 
             *pow(gcr[i_species].Ekin[ip],2);


           i++;
   }//ir
   }//iz
   }//ip
   }//i_species
    }//n_spatial_dimensions==2

    if(gcr[0].n_spatial_dimensions==3){
 
   int i=0;
   for (int i_species=0;i_species< naxes[4];i_species++){ 
   for (int ip       =0;        ip<naxes[3];       ip++){
   for (int iz       =0;        iz<naxes[2];       iz++){
   for (int iy       =0;        iy<naxes[1];       iy++){
   for (int ix       =0;        ix<naxes[0];       ix++){

     if(gcr[i_species].A!=0)
     array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz].s[ip] 
             *    gcr[i_species].A
             *pow(gcr[i_species].Ekin[ip],2);


     if(gcr[i_species].A==0)     // electrons, positrons
     array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz].s[ip]
             *pow(gcr[i_species].Ekin[ip],2);

           i++;
   }//ix
   }//iy
   }//iz
   }//ip
   }//i_species
    }//n_spatial_dimensions==3

    /* Write the array of floats to the image */
    fits_write_img(fptr, TFLOAT, fpixel, nelements, array, &status);


    // write basic FITS keywords

    double crval1,crval2,crval3,crval4,crval5;
    double cdelt1,cdelt2,cdelt3,cdelt4,cdelt5;

    if(gcr[0].n_spatial_dimensions==2){
    crval1=gcr[0].r[0];
    crval2=gcr[0].z[0];
    crval3=log10(gcr[0].Ekin[0]);
    crval4=1;
    cdelt1=gcr[0].dr;  
    cdelt2=gcr[0].dz;    
    cdelt3=log10(gcr[0].Ekin_factor   );
    cdelt4=1;
    }

    if(gcr[0].n_spatial_dimensions==3){
    crval1=gcr[0].x[0];
    crval2=gcr[0].y[0];
    crval3=gcr[0].z[0];
    crval4=log10(gcr[0].Ekin[0]);
    crval5=1;
    cdelt1=gcr[0].dx;
    cdelt2=gcr[0].dy;    
    cdelt3=gcr[0].dz;    
    cdelt4=log10(gcr[0].Ekin_factor   );
    cdelt5=1;
      }


    fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
    fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
    fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
    fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
    if(gcr[0].n_spatial_dimensions==3)
    fits_update_key(fptr, TDOUBLE, "CRVAL5", &crval5,"Start of axis 5", &status);

    fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
    fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
    fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
    fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
    if(gcr[0].n_spatial_dimensions==3)
    fits_update_key(fptr, TDOUBLE, "CDELT5", &cdelt5,"Increment of axis 5", &status);

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

   sprintf(keyword,"NUCK%03d",       i_nucleus+1 ); // e.g. NUCK012                      //AWS20010731                
   sprintf(comment,"K-electrons of nucleus %d",i_nucleus+1 );                            //AWS20010731
      cout<<keyword<<" "<<gcr[i_nucleus].K_electron<<endl;                               //AWS20010731
   fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].K_electron,comment, &status);//AWS20010731
 }


 // write normalization keywords   AWS20010121
 for(int i_nucleus=0;i_nucleus<n_species;i_nucleus++){
    if(gcr[i_nucleus].A==1 && gcr[i_nucleus].Z==1)
        fits_update_key(fptr, TDOUBLE, "NUCNORM",
	 &galdef.source_normalization,"Nucleon  norm. factor", &status); // IMOS20030217
//         &gcr[i_nucleus].normalization_factor,"Nucleon  norm. factor", &status);

    if(gcr[i_nucleus].A==0 && gcr[i_nucleus].Z==-1)
        fits_update_key(fptr, TDOUBLE, "ELENORM",
         &galdef.electron_source_normalization,"Electron norm. factor", &status); // IMOS20031016
//         &gcr[i_nucleus].normalization_factor,"Electron norm. factor", &status);

 }

    fits_close_file(fptr, &status);            /* close the file */

    fits_report_error(stderr, status);  /* print out any error messages */
    return( status );


  cout<<" <<<< store_gcr_full"<<endl;
return stat;
}

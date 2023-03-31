
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_gcr.cc *                           galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"global.h"

#include "fitsio.h" 
int read_gcr(){
int stat;

double NUCNORM,ELENORM;
char   comment[100];

  cout<<" >>>> read_gcr"<<endl;
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

    char  infile[100];
  
    strcpy( infile,configure.fits_directory);
 
    strcat( infile,"nuclei_full_");
    strcat( infile,galdef.galdef_ID);
    cout<<"  reading  nuclei array from file "<<infile<<endl;

    status = 0;         /* initialize status before calling fitsio routines */
    fits_open_file(&fptr,  infile,READONLY, &status);   
    cout<<"  fits open status = "<<status<<endl;



    /* Read the array of floats  */
    float nulval=0;
    int anynul;
    fits_read_img(fptr, TFLOAT, fpixel, nelements, &nulval,array, &anynul,&status);
 


    if(gcr[0].n_spatial_dimensions==2){
   int i=0;
   for (int i_species=0;i_species< naxes[3];i_species++){ 
   for (int ip       =0;        ip<naxes[2];       ip++){
   for (int iz       =0;        iz<naxes[1];       iz++){
   for (int ir       =0;        ir<naxes[0];       ir++){


      /*this is how the arrays are output, for reference
     if(gcr[i_species].A!=0)  
     array[i]=    gcr[i_species].cr_density.d2[ir]    [iz].s[ip] 
             *    gcr[i_species].A
             *pow(gcr[i_species].Ekin[ip],2);

     if(gcr[i_species].A==0)  // electrons, positrons
     array[i]=    gcr[i_species].cr_density.d2[ir]    [iz].s[ip] 
             *pow(gcr[i_species].Ekin[ip],2);
	     */


     if(gcr[i_species].A!=0)
         gcr[i_species].cr_density.d2[ir]    [iz].s[ip]=array[i]
             /    gcr[i_species].A
             /pow(gcr[i_species].Ekin[ip],2);

     if(gcr[i_species].A==0)  // electrons, positrons
       gcr[i_species].cr_density.d2[ir]    [iz].s[ip]= array[i]
             /pow(gcr[i_species].Ekin[ip],2);

           i++;
   }//ir
   }//iz
   }//ip
     if(galdef.verbose==1)gcr[i_species].cr_density.print();
   }//i_species
    }//n_spatial_dimensions==2






    if(gcr[0].n_spatial_dimensions==3){
 
   int i=0;
   for (int i_species=0;i_species< naxes[4];i_species++){ 
   for (int ip       =0;        ip<naxes[3];       ip++){
   for (int iz       =0;        iz<naxes[2];       iz++){
   for (int iy       =0;        iy<naxes[1];       iy++){
   for (int ix       =0;        ix<naxes[0];       ix++){

     /* this is how the arrays are output, for reference
     if(gcr[i_species].A!=0)
     array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz].s[ip] 
             *    gcr[i_species].A
             *pow(gcr[i_species].Ekin[ip],2);


     if(gcr[i_species].A==0)     // electrons, positrons
     array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz].s[ip]
             *pow(gcr[i_species].Ekin[ip],2);
	     */


     if(gcr[i_species].A!=0)
          gcr[i_species].cr_density.d3[ix][iy][iz].s[ip]=array[i]
             /    gcr[i_species].A
             /pow(gcr[i_species].Ekin[ip],2);


     if(gcr[i_species].A==0)     // electrons, positrons
          gcr[i_species].cr_density.d3[ix][iy][iz].s[ip]=array[i]
             /pow(gcr[i_species].Ekin[ip],2);

           i++;
   }//ix
   }//iy
   }//iz
   }//ip
     if(galdef.verbose==1)gcr[i_species].cr_density.print();
   }//i_species


    }//n_spatial_dimensions==3


    // remove normalization factor to allow further processing
    // for a warm start

 

   for (int i_species=0;i_species< n_species;i_species++){ 
     //nuclei 
    if(gcr[i_species].A!=0){
     fits_read_key(fptr,TDOUBLE,"NUCNORM" ,&NUCNORM ,comment,&status);  
     cout<<"read NUCNORM status= "<<status<<"  NUCNORM="<<NUCNORM<<endl;
     gcr[i_species].cr_density/=NUCNORM;
    }
    // positrons
    if(gcr[i_species].A==0&&gcr[i_species].Z==+1){
     fits_read_key(fptr,TDOUBLE,"NUCNORM" ,&NUCNORM ,comment,&status);  
     cout<<"read NUCNORM status= "<<status<<"  NUCNORM="<<NUCNORM<<endl;
     gcr[i_species].cr_density/=NUCNORM;
    }    
    // electrons
    if(gcr[i_species].A==0&&gcr[i_species].Z==-1){
     fits_read_key(fptr,TDOUBLE,"ELENORM" ,&ELENORM ,comment,&status);  
     cout<<"read ELENORM status= "<<status<<" ELENORM="<<ELENORM<<endl;
     gcr[i_species].cr_density/=ELENORM;
    }
    }//i_species
 

    fits_close_file(fptr, &status);            /* close the file */

    fits_report_error(stderr, status);  /* print out any error messages */


    delete[] array; //AWS20010216

  cout<<" <<<< read_gcr"<<endl;

    return( status );
return stat;
}


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_COR.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"
#include"fitsio.h"


int Galprop::read_COR()
{
   cout<<" >>>> read_COR    "<<endl;
   int status=0;

   fitsfile *fptr;
   char COR_filename[200];
   int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
   float CRVAL1,CRVAL2,CRVAL3;
   float CDELT1,CDELT2,CDELT3;
   char comment[100];

   Distribution COR_input; // CO rings as read from FITS file
   Distribution nuse;      // number of cells used in rebinned map


   strcpy(COR_filename,configure.fits_directory);
   //   strcat(COR_filename,"COMPASS.COR.M0010453");// 6 rings Ro=10 kpc
   if(galdef.CO_survey==8)                     //AWS20050913
   strcat(COR_filename,"COMPASS.COR.M0010465");// 8 rings Ro=8.5 kpc Digel
   if(galdef.CO_survey==9)                     //AWS20050913
// strcat(COR_filename,"COMPASS.COR.M0010465_9rings");// 8 rings Ro=8.5 kpc Digel, modified to 9 ring format
// strcat(COR_filename,"rbands_co2.fits"            );// 9 rings Ro=8.5 kpc Digel 13 Jan 2006  //AWS20060113
// strcat(COR_filename,"rbands_co3.fits"            );// 9 rings Ro=8.5 kpc Digel 16 Jan 2006  //AWS20060116
   strcat(COR_filename,"rbands_co4.fits"            );// 9 rings Ro=8.5 kpc Digel 18 Jan 2006  //AWS20060118

   cout<<" reading COR from "<<COR_filename<<endl;

   if( fits_open_file(&fptr,COR_filename,READONLY,&status) ) cout<<"read COR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read COR status= "<<status<<endl;


   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read COR status= "<<status<<endl;
 
   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read COR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read COR status= "<<status<<endl;
 

   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2   = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2   = "<<CDELT1<<" "<<CDELT2<<endl;

   long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
   float *COR_in=new float[nelements];
   float nulval=0;
   int anynul;

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,COR_in,&anynul,&status) )
      cout<<"#read COR status= "<<status<<endl;

   // for(int i=0; i<nelements; i++) cout<<COR_in[i]<<" ";

   cout<<"generating galaxy.COR:"<<endl;


         int i_long,i_lat,i_Ring;
         int i_long_in,i_lat_in;
         int n_long_in=NAXIS1;
         int n_lat_in =NAXIS2;
         int n_Ring=NAXIS3;

         
      
         COR_input.init(n_long_in,n_lat_in,n_Ring,1);

         cout<<" COR_input initialized"<<endl;


         int i=0;

         for  ( i_Ring=0;  i_Ring<n_Ring; i_Ring++){
	 
	  for ( i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            

                  COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= COR_in[i];
                  i++;

          
            }  //  i_long_in
           }   //  i_lat_in
	 }     //  i_Ring
   



	 if(galdef.verbose== -301 )//selectable debug
      {
      
                COR_input.print();
      }
  



double l,b;

// create column density map corresponding to required gamma skymaps
// by rebinning input map
// An extra ring is added for compatibility with HI rings (6 instead of the input 5) [COR-100453 only]
// and this CO ring is left at zero.

      galaxy.COR.init(galaxy.n_long,galaxy.n_lat,n_Ring,1);

      nuse.init      (galaxy.n_long,galaxy.n_lat,n_Ring,1);


         for  (i_Ring=0;     i_Ring   <n_Ring; i_Ring++)   {
       	  for (i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
           for(i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            
                  l=CRVAL1+(i_long_in)*CDELT1;
                  b=CRVAL2+(i_lat_in )*CDELT2;
                  i_long= (int) ((l-galaxy.long_min)/galaxy.d_long); // IMOS20010816
                  i_lat = (int) ((b-galaxy.lat_min )/galaxy.d_lat);  // IMOS20010816

                  if(i_long>=0 && i_long<galaxy.n_long&&i_lat >=0 && i_lat<galaxy.n_lat){

                  galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]+= COR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
                  nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;


		  //       cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
		  }
                  

          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring

	 // normalize by number of cells used in rebinning
        for  (i_Ring=0;  i_Ring<       n_Ring; i_Ring++){
	 for (i_lat =0;  i_lat <galaxy.n_lat;  i_lat++ ){
          for(i_long=0;  i_long<galaxy.n_long; i_long++){

            if(nuse.d3      [i_long][i_lat][i_Ring].s[0]>0)
               galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]/=nuse.d3      [i_long][i_lat][i_Ring].s[0];
          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring



       delete[] COR_in;
       COR_input.delete_array();  
       nuse.delete_array();




	if(galdef.verbose==-303)// selectable debug
      {
        cout<<"read_COR: galaxy.COR:"<<endl;
                galaxy.COR.print();
      }
   cout<<" <<<< read_COR    "<<endl;
   return status;
}























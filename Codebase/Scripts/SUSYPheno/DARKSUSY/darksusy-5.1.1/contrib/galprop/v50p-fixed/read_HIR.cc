
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_HIR.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"
#include"fitsio.h"

//* JE FIX: Following line(s) added for Snow Leopard
#include<cstring>
#include<cstdlib>

int Galprop::read_HIR()
{
   cout<<" >>>> read_HIR    "<<endl;
   int status=0;

   fitsfile *fptr;
   char HIR_filename[200];
   int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
   float CRVAL1,CRVAL2,CRVAL3;
   float CDELT1,CDELT2,CDELT3;
   char comment[100];

   Distribution HIR_input; // HI rings as read from FITS file
   Distribution nuse;      // number of cells used in rebinned map

   //   cout<<"galaxy.n_long,n_lat  "<<galaxy.n_long<<" "<<galaxy.n_lat<<endl;

   strcpy(HIR_filename,configure.fits_directory);
   // strcat(HIR_filename,"COMPASS.HIR.M0010051");// 6 rings Ro=10 kpc
   if(galdef.HI_survey==8)                                                //AWS20050913
      strcat(HIR_filename,"COMPASS.HIR.M0010058");// 8 rings Digel Ro=8.5 kpc
   if(galdef.HI_survey==9)
//    strcat(HIR_filename,"rbands_hi2_32.fits");// 9 rings Digel Aug  2005  AWS20050826
//    strcat(HIR_filename,"rbands_hi5.fits"   );// 9 rings Digel Jan  2006  AWS20060112
      strcat(HIR_filename,"rbands_hi7.fits"   );// 9 rings Digel Jan  2006  AWS20060116

   cout<<" reading HIR from "<<HIR_filename<<endl;

   if( fits_open_file(&fptr,HIR_filename,READONLY,&status) ) cout<<"read HIR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read HIR status= "<<status<<endl;
 

   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read HIR status= "<<status<<endl;

   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read HIR status= "<<status<<endl;


   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3   = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2   = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2   = "<<CDELT1<<" "<<CDELT2<<endl;

   long nelements=NAXIS1*NAXIS2*NAXIS3, felement=1;
   float *HIR_in=new float[nelements];
   float nulval=0;
   int anynul;

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,HIR_in,&anynul,&status) )
      cout<<"#read HIR status= "<<status<<endl;

   // for(int i=0; i<nelements; i++) cout<<HIR_in[i]<<" ";

   cout<<"generating galaxy.HIR:"<<endl;


         int i_long,i_lat,i_Ring;
         int i_long_in,i_lat_in;
         int n_long_in=NAXIS1;
         int n_lat_in =NAXIS2;
         int n_Ring=NAXIS3;

         
      
         HIR_input.init(n_long_in,n_lat_in,n_Ring,1);

         cout<<" HIR_input initialized"<<endl;


         int i=0;

         for  ( i_Ring=0;  i_Ring<n_Ring; i_Ring++){
	 
	  for ( i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            

                  HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= HIR_in[i];
                  i++;

          
            }  //  i_long_in
           }   //  i_lat_in
	 }     //  i_Ring
   

         HIR_input*= 1.0e20;   // applies to HIR.10058, not HIR.10051   

	 if(galdef.verbose== -201) // selectable debug
      {
      
                HIR_input.print();
      }
  



double l,b;

// create column density map corresponding to required gamma skymaps
// by rebinning input map
// cout<<"galay.HIR.init( "<<galaxy.n_long<<" "<<galaxy.n_lat<<" "<<n_Ring<<",1)"<<endl;

      galaxy.HIR.init(galaxy.n_long,galaxy.n_lat,n_Ring,1);

      nuse.init      (galaxy.n_long,galaxy.n_lat,n_Ring,1);


         for  (i_Ring=0;     i_Ring   <n_Ring; i_Ring++)   {
       	  for (i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
           for(i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            
                  l=CRVAL1+(i_long_in)*CDELT1;
                  b=CRVAL2+(i_lat_in )*CDELT2;
                  i_long= (int) ((l-galaxy.long_min)/galaxy.d_long); // IMOS20010816
                  i_lat = (int) ((b-galaxy.lat_min )/galaxy.d_lat);  // IMOS20010816

                  if(i_long>=0 && i_long<galaxy.n_long&&i_lat >=0 && i_lat<galaxy.n_lat){

                  galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+= HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0];
                  nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;


		  //       cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
		  }
                  

          
            }  //  i_long
           }   //  i_lat
	 }     //  i_Ring



	 

	/////////////////////////////////// |b|>39.75

  if(galdef.HI_survey==8)// only for 8 ring survey; 9 ring survey is all-sky  AWS20050914
  {

   strcpy(HIR_filename,configure.fits_directory);
   //   strcat(HIR_filename,"COMPASS.HIR.M0010049");// |b|>24.75, *0.5
   strcat(HIR_filename,"COMPASS.HIR.M0010055");// |b|>39.75  S.Digel
     
   cout<<" reading HIR from "<<HIR_filename<<endl;

   if( fits_open_file(&fptr,HIR_filename,READONLY,&status) ) cout<<"read HIR open status= "<<status<<endl;

   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read HIR status= "<<status<<endl;
 

   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read HIR status= "<<status<<endl;
;
   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read HIR status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read HIR status= "<<status<<endl;


   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3   = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<endl;
   cout<<" CRVAL1,2     = "<<CRVAL1<<" "<<CRVAL2<<endl;
   cout<<" CDELT1,2     = "<<CDELT1<<" "<<CDELT2<<endl;

        nelements=NAXIS1*NAXIS2*NAXIS3;
          HIR_in=new float[nelements];
    

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,HIR_in,&anynul,&status) )
      cout<<"#read HIR status= "<<status<<endl;

         n_long_in=NAXIS1;
         n_lat_in =NAXIS2;
      
         HIR_input.init(n_long_in,n_lat_in,n_Ring,1);

         cout<<" HIR_input initialized for high latitudes"<<endl;


         i=0;

 
	 //         i_Ring=2            ;// local ring (0-4-8-10-12-15-30)//HIR-10051

	                       //                0   1   2   3V    4    5    6    7
          i_Ring=3            ;// local ring (1.5-3.5-5.5-7.5 - 9.5-11.5-13.5-15.5-50) HIR-10055

	  for ( i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
             
	   for( i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            

                  HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]= HIR_in[i];
                  i++;

          
            }  //  i_long_in
           }   //  i_lat_in
 
          HIR_input*= 1.0e20;   // applies to HIR.10055, not HIR.10049   

	  if(galdef.verbose== -202) // selectable debug
      {
      
                HIR_input.print();
      }

    

       	  for (i_lat_in =0;  i_lat_in <n_lat_in;  i_lat_in++ ){
           for(i_long_in=0;  i_long_in<n_long_in; i_long_in++){
            
                  l=CRVAL1+(i_long_in)*CDELT1;
                  b=CRVAL2+(i_lat_in )*CDELT2;

                  if(b<-39.75 || b>+39.75){
		  i_long= (int) ((l-galaxy.long_min)/galaxy.d_long); // IMOS20010816
                  i_lat = (int) ((b-galaxy.lat_min )/galaxy.d_lat);  // IMOS20010816

                  if(i_long>=0 && i_long<galaxy.n_long&&i_lat >=0 && i_lat<galaxy.n_lat){

	                                                            
                  galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+= HIR_input.d3[i_long_in][i_lat_in][i_Ring].s[0]    ;
                  nuse.d3      [i_long][i_lat][i_Ring].s[0]+=1;

		  }
		  //	         cout<<"l b i_long_in i_lat_in i_long i_lat i_Ring "<<l<<" "<< b<<" "<< i_long_in<<" "<< i_lat_in<<" "<< i_long<<" "<< i_lat<<" "<< i_Ring<<endl;
		  }
                  

          
            }  //  i_long
           }   //  i_lat
       

  }// if (HI_survey==8)

	//////////////////////////////////////////////////

	 

	 // normalize by number of cells used in rebinning
        for  (i_Ring=0;  i_Ring<       n_Ring; i_Ring++){
	 for (i_lat =0;  i_lat <galaxy.n_lat;  i_lat++ ){
          for(i_long=0;  i_long<galaxy.n_long; i_long++){

            if(nuse.d3      [i_long][i_lat][i_Ring].s[0]>0)
               galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]/=nuse.d3      [i_long][i_lat][i_Ring].s[0];
          
            }  //  i_long
           }   //  i_lat
	}//i_Ring

  



       delete[] HIR_in;
       HIR_input.delete_array();  
       nuse.delete_array();


	//////////////////////////////////////////////////

       if(galdef.verbose== -203)// selectable debug
      {
        cout<<"read_HIR: galaxy.HIR:"<<endl;
                galaxy.HIR.print();
      }
   cout<<" <<<< read_HIR    "<<endl;
   return status;
}

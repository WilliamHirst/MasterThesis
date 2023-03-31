
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_isrf.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"
#include"fitsio.h"


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//indexing of input ISRF array as read into 1D array

int isrf_index(int ir,int iz,int inu,int icomp,int nr,int nz,int nnu,int ncomp)
{
   return   icomp  *nnu *nz*nr
          + inu         *nz*nr
          + iz             *nr
          + ir;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/* 
ISRF is in Hz eV cm-3 Hz-1 (or micron eV cm^-3 micron^-1)
integral  energy density Hz-1  d(nu) = integral (nu* energy density Hz-1) d(log nu)
d(log nu) is constant in this ISRF
factor=  LOG(nu(2)/nu(1)) 
*/

int Galprop::read_isrf()
{
   cout<<" >>>> read_isrf    "<<endl;
   int status=0;

   fitsfile *fptr;
   char ISRF_filename[200];
   int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
   float CRVAL1,CRVAL2,CRVAL3;
   float CDELT1,CDELT2,CDELT3;
   char comment[100];

   strcpy(ISRF_filename,configure.fits_directory);
   //strcat(ISRF_filename,"isrf_interp_04_000015");
   //strcat(ISRF_filename,"porter_ISRF.fits");           //AWS20050225
   //strcat(ISRF_filename,"porter_RFScattering10kpc.fits");//AWS20050301
   strcat(ISRF_filename,galdef.ISRF_file);//AWS20050301


   cout<<" reading ISRF from "<<ISRF_filename<<endl;

   if( fits_open_file(&fptr,ISRF_filename,READONLY,&status) ) cout<<"read isrf open status= "<<status<<endl;
   if(fptr==NULL)exit(-2);


   if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) cout<<"0read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"1read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"2read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) cout<<"3read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TINT,"NAXIS4",&NAXIS4,comment,&status) ) cout<<"4read isrf status= "<<status<<endl;

   if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"5read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"6read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) cout<<"7read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"8read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"9read isrf status= "<<status<<endl;
   if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) cout<<"/read isrf status= "<<status<<endl;

   cout<<" NAXIS = "<<NAXIS <<endl;
   cout<<" NAXIS1,2,3,4 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<" "<<NAXIS4<<endl;
   cout<<" CRVAL1,2,3 = "<<CRVAL1<<" "<<CRVAL2<<" "<<CRVAL3<<endl;
   cout<<" CDELT1,2,3 = "<<CDELT1<<" "<<CDELT2<<" "<<CDELT3<<endl;

   long nelements=NAXIS1*NAXIS2*NAXIS3*NAXIS4, felement=1;
   float *isrf_in=new float[nelements];
   float nulval=0;
   int anynul;

   if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,isrf_in,&anynul,&status) )
      cout<<"#read isrf status= "<<status<<endl;

// for(int i=0; i<nelements; i++) cout<<isrf_in[i]<<" ";

   cout<<"generating galaxy.ISRF:"<<endl;

   galaxy.n_ISRF_components    =NAXIS4;
   galaxy.ISRF=new Distribution[NAXIS4];

   for(int i=0; i<NAXIS4; i++)
   {
      if(galdef.n_spatial_dimensions==2)              // ==== 2D ====
      {
         galaxy.ISRF[i].init(galaxy.n_rgrid,galaxy.n_zgrid,NAXIS3);
         cout<<" galaxy.ISRF initialized with frequency axis dimension="
            <<galaxy.ISRF[i].n_pgrid<<endl;

         for(int ir=0; ir<galaxy.n_rgrid; ir++)
         {    
            for(int iz=0; iz<galaxy.n_zgrid; iz++)
            {
	      int irr=(int)((     galaxy.r[ir] -CRVAL1) /CDELT1+0.5);//IMOS20060420
	      int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
               if(irr>NAXIS1-2) irr=NAXIS1-2;
               if(izz>NAXIS2-2) izz=NAXIS2-2;
               float rr=CRVAL1+irr*CDELT1;
               float zz=CRVAL2+izz*CDELT2;

// cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;

               for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++)
               {
                  float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                  float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                  float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                  float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                  float v5=v1+(v2-v1)*(galaxy.r[ir]-rr)/CDELT1;
                  float v6=v3+(v4-v3)*(galaxy.r[ir]-rr)/CDELT1;
                  float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
                  if(value<0.0) value=0.0;
// reverse scale from wavelength to frequency
                  galaxy.ISRF[i].d2[ir][iz].s[galaxy.ISRF[i].n_pgrid-1-inu]= value;
    
// cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"   "<<irr<<" "<<izz<<"   "<<rr<<" "<<zz<<" "<< v5+(v6-v5)*(galaxy.z[iz]-zz)/CDELT2<<endl;
               }  //  inu
            }  //  iz
         }  //  ir
      }  //  2D

      if(galdef.n_spatial_dimensions==3)              // ==== 3D ====
      {
         galaxy.ISRF[i].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid,NAXIS3);
         for(int ix=0; ix<galaxy.n_xgrid; ix++)
         {  
            for(int iy=0; iy<galaxy.n_ygrid; iy++)
            {
               for(int iz=0; iz<galaxy.n_zgrid; iz++)
               {
                  float r=sqrt(pow(galaxy.x[ix],2)+pow(galaxy.y[iy],2));
                  int irr=(int)((            r     -CRVAL1) /CDELT1+0.5);//IMOS20060420
                  int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
                  if(irr>NAXIS1-2) irr=NAXIS1-2;
                  if(izz>NAXIS2-2) izz=NAXIS2-2;

                  float rr=CRVAL1+irr*CDELT1;
                  float zz=CRVAL2+izz*CDELT2;

// cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;

                  for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++)
                  {
                     float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                     float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                     float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                     float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
                     float v5=v1+(v2-v1)*(       r    -rr)/CDELT1;
                     float v6=v3+(v4-v3)*(       r    -rr)/CDELT1;
                     float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
                     if(value<0.0) value=0.0;
// reverse scale from wavelength to frequency
                     galaxy.ISRF[i].d3[ix][iy][iz].s[galaxy.ISRF[i].n_pgrid-1-inu]= value;
                  }  //  inu   
               }  //  ix
            }  //  iy
         }  //  iz
      }  //  3D

      if(galdef.verbose>=10)
      {
         cout<<"ISRF component "<<i+1<<endl;
         galaxy.ISRF[i].print();
      }
   }  // isrf component i

// Create the array of ISRF frequencies
// using wavelength in microns for axis 3 of input ISRF on log10 scale.
// Reverse scale so that frequency increases.

   galaxy.nu_ISRF=new double[galaxy.ISRF[0].n_pgrid];

// microns -> cm; nu=c/lambda
   for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
      galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid-1-inu]=c/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4);

   cout<<" ISRF frequency grid (Hz):"<<endl;
   for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) 
      cout<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu]<<endl;

   delete[] isrf_in;//AWS20010216

   cout<<" <<<< read_isrf    "<<endl;
   return status;
}

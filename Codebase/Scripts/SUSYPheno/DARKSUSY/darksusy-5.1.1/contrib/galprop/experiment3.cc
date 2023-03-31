//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// This file contains 4 routines on c++ which read the nuclear fits files,
// read all headers of a fits file, and print out cfitsio error messages,
// prepares a file to plot with gnuplot (v.3.7 and higher).
// The program provides a simple dialog, the file(s) to plot should be given
// in a command line.
//                            ### Igor Moskalenko, NASA/GSFC ###  9/24/2001 ###
// readimage    -reads nuclear fits files created by galprop package;
// make_gnufile -creates gnuplot file "tmp.gnu" for plotting;
// readheader   -reads all headers of a fits file;
// printerror   -prints out cfitsio error messages and exit program;
// To compile:
// g++ experiment3.cc ../libs/*.a -I/software/lheasoft/release/Linux_2.2_i686/include -I../libs/
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <fstream.h>
#include <iostream.h>
#include <math.h>
#include "fitsio.h"  // fitsio header file

void readheader( char* );
void readimage( char*,FILE*,int );
void printerror( int status);
void make_gnufile(int,char**,int);
float add_nuclei(int,long,long,float,float);

float Ekin=2000.,       // MeV/nucleon -kinetic energy
      R=8.5,            // kpc -solar position
      Phi=500.,Phi1=0., // MV -modulation potential
      ymax=0., add_nuclei_fraction, add_nuclei_energy, add_nuclei_escale, norm_nuc;
int z[10]={0},a[10]={0},z1[10]={0},a1[10]={0}, add_nuclei_case;
int datakey=0, spec_type=0, spec_units=0;
char tmp[50];

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

main(int argc, char*argv[])
{
   FILE *fout;
   int i,key=9;

   if(argc<2) 
   {
      cout<<"\nexample: plot_nuclei 41p_123456 41p_000111\n";
      cout<<"\nno input file specified: exit\n";
      exit(1);
   }

   cout<<"\nWhat's your wish ? (choose one option):\n\n";
   cout<<"*++++++++++++++++++++++++++++++++++++++*\n"
       <<"+ Plot all nuclei in one plot      = 0 +\n"
       <<"+ Plot spectra of certain species  = 1 +\n"
       <<"+ Plot ratio of certain species    = 2 +\n"
       <<"+ Plot electrons and positrons     = 3 +\n"
       <<"+ Exit now                         = 9 +\n"
       <<"*++++++++++++++++++++++++++++++++++++++*\n\n";
   cin>>key;

   switch(key)
   {
      case 0:
         cout<<"\nYour choise: Plot all nuclei\n";
         cout<<"\nChoose one option\n\n";
         cout<<"*++++++++++++++++++++++++++++*\n"
             <<"+ All isotopes           = 0 +\n"
             <<"+ All elements           = 1 +\n"
             <<"+ Isotopes of an element = 2 +\n"
             <<"*++++++++++++++++++++++++++++*\n";
         cin>>datakey;
         if(datakey != 0 && datakey != 1 && datakey != 2) datakey=0;
         if(datakey == 2)
	 {
	    cout<<"give the nuclear charge Z\n";
            cin>>z[0];
         }
         cout<<"give modulation potential, Phi [MV]; =0 for interstellar medium\n";
         cin>>Phi;
         cout<<"give kinetic energy after modulation, Ekin [MeV/nucleon] \n";
         cin>>Ekin;
         if(Ekin<0. || Phi<0. ) 
         {
	    cout<<"\nNegative number: Phi ="<<Phi<<" Ekin ="<<Ekin<<"\nexit \n";
            exit(1);
         }
         break;

      case 1:
         cout<<"\nYour choise: Plot spectra of certain species\n";
         cout<<"specify your choise: Z A  Z A etc. (0 0 -to stop)\n";
         cout<<"for a plot of an element type Z 0 0 0; for antiprotons type -1 1 0 0\n";
         cout<<"example1: 8 16 8 17 8 18 0 0  -to plot Oxygen \n";
         cout<<"example2: 8 0 0 0             -to plot Oxygen \n";
         for(i=0; i<10; i++)
         { 
            cin>>z[i]>>a[i];
            if(z[i] ==0 && a[i] ==0) break;
         }
         if(i==0)
         {
            cout<<"Wrong numbers: z= "<<z[0]<<" a= "<<a[0]<<endl;
            exit(1);
         }
         datakey = z[0];
         cout<<"give modulation potential, Phi [MV]; =0 for interstellar medium\n";
         cin>>Phi;
         cout<<"Spectrum presentation: 0 -Flux; 2 -E^2 Flux\n";
         cin>>spec_type;
         if(spec_type)  cout<<"Units: 0  -MeV/(cm2 s sr) vs. MeV\n"
                            <<"       1  -GeV/(m2  s sr) vs. GeV\n";
         else  cout<<"Units: 0  -1/(cm2 s sr MeV) vs. MeV\n"
                   <<"       1  -1/(m2  s sr GeV) vs. GeV\n";
         cin>>spec_units;
	 break;

      case 2:
         cout<<"\nYour choise: Plot a ratio (Z.A+...)/(Z1.A1+...)\n";
         cout<<"give numerator nuclei: Z A  Z A etc. (0 0 -to stop)\n";
         cout<<"for a plot of a certain ratio type:\n";
         cout<<"*+++++++++++++++++++++++++*\n"
             <<"+ -1 0 0 0    = pbar/p    +\n"
             <<"+  2 0 0 0    = He3/He4   +\n"
             <<"+  4 0 0 0    = Be10/Be9  +\n"
             <<"+  5 0 0 0    = B/C       +\n"
             <<"+ 13 0 0 0    = Al26/Al27 +\n"
             <<"+ 17 0 0 0    = Cl36/Cl   +\n"
             <<"+ 25 0 0 0    = Mn54/Mn   +\n"
             <<"+ 26 0 0 0    = Al26/Si28 +\n"
             <<"+ 36 0 0 0    = Cl36/Fe56 +\n"
             <<"+ 54 0 0 0    = Mn54/Fe56 +\n"
             <<"+ 56 0 0 0    = subFe/Fe  +\n"
             <<"*+++++++++++++++++++++++++*\n\n";
         for(i=0; i<10; i++)
         { 
            cin>>z[i]>>a[i];
            if(z[i] ==0 && a[i] ==0) break;
         }
         if(i==0)
         {
            cout<<"Wrong numbers: z= "<<z[0]<<" a= "<<a[0]<<endl;
            exit(1);
         }
         datakey = z[0];
         if(a[0]!=0)
         {
            cout<<"give denominator nuclei: Z1 A1  Z1 A1 etc. (0 0 -to stop)\n";
            for(i=0; i<10; i++)
            { 
               cin>>z1[i]>>a1[i];
               if(z1[i] ==0 || a1[i] ==0) break;
            }
            datakey = 0;
            if(i==0)
            {
               cout<<"Wrong numbers: z1= "<<z1[0]<<" a1= "<<a1[0]<<endl;
               exit(1);
            }
         } else 
         {
	    if(-1 ==z[0]) // pbar/p
	    {
	                   a [0]=  1;  z [1]= -1;  a [1]=  1;
               z1[0]=  1;  a1[0]=  1;  z1[1]=  1;  a1[1]=  1;
            }
	    if( 2 ==z[0]) // He3/He4
	    {
	                   a [0]=  3;
               z1[0]=  2;  a1[0]=  4;
            }
	    if( 4 ==z[0]) // Be10/Be9
	    {
	                   a [0]= 10;
               z1[0]=  4;  a1[0]=  9;
            }
	    if( 5 ==z[0]) // B/C
	    {
	                   a [0]= 10;  z [1]=  5;  a [1]= 11;
               z1[0]=  6;  a1[0]= 12;  z1[1]=  6;  a1[1]= 13;
            }
	    if(13 ==z[0]) // Al26/Al27
	    {
	                   a [0]= 26;
               z1[0]= 13;  a1[0]= 27;
            }
	    if(17 ==z[0]) // Cl36/Cl
	    {
	                   a [0]= 36;
               z1[0]= 17;  a1[0]= 35;  z1[1]= 17;  a1[1]= 36;  z1[2]= 17;  a1[2]= 37;
            }
	    if(25 ==z[0]) // Mn54/Mn
	    {
	                   a [0]= 54;
               z1[0]= 25;  a1[0]= 53;  z1[1]= 25;  a1[1]= 54;  z1[2]= 25;  a1[2]= 55;
            }
	    if(26 ==z[0]) // Al26/Si28
	    {
	       z [0]= 13;  a [0]= 26;
               z1[0]= 14;  a1[0]= 28;
            }
	    if(36 ==z[0]) // Cl36/Fe56
	    {
	       z [0]= 17;  a [0]= 36;
               z1[0]= 26;  a1[0]= 56;
            }
	    if(54 ==z[0]) // Mn54/Fe56
	    {
	       z [0]= 25;  a [0]= 54;
               z1[0]= 26;  a1[0]= 56;
            }
	    if(56 ==z[0]) // subFe/Fe 
	    {
	       z [0]= 21;  a [0]= 45; 
               z [1]= 22;  a [1]= 44;  z [2]= 22;  a [2]= 46;  z [3]= 22;  a [3]= 47;  
               z [4]= 22;  a [4]= 48;  z [5]= 22;  a [5]= 49;  z [6]= 22;  a [6]= 50;
	       z [7]= 23;  a [7]= 49;  z [8]= 23;  a [8]= 50;  z [9]= 23;  a [9]= 51; 
               z1[0]= 26;  a1[0]= 54;  z1[1]= 26;  a1[1]= 55;  z1[2]= 26;  a1[2]= 56;  
               z1[3]= 26;  a1[3]= 57;  z1[4]= 26;  a1[4]= 58; // a1=60 is not included
            }
	 }
         cout<<"give modulation potential, Phi [MV]; =0 for interstellar medium\n";
         cin>>Phi;
	 if(-1 ==z[0]) // pbar/p
	 {
            cout<<"give modulation potential for protons\n";
            cin>>Phi1;
         }
         cout<<"Units: 0  -MeV\n"
             <<"       1  -GeV\n";
         cin>>spec_units;
         break;

      case 3:
         cout<<"\nYour choise:  Plot electrons and positrons\n";
         cout<<"You have the following options:\n";
         cout<<"*+++++++++++++++*\n"
             <<"* electrons  =-1 +\n"
             <<"* positrons  =+1 +\n"
             <<"* e+/e_tot   = 0 + under construction \n"
             <<"*+++++++++++++++*\n\n";
         cin>>z[0];
         datakey = z[0];
         if(z[0]<-1 || z[0]>1)
         {
            cout<<"Wrong number: "<<z[0]<<endl;
            exit(1);
         }
         if(z[0]==0)
	 {
            cout<<"\nThis option is under construction: exit\n";
	    exit(0);
         }
         cout<<"give modulation potential, Phi [MV]; =0 for interstellar medium\n";
         cin>>Phi;
         cout<<"Spectrum presentation: 0 -Flux; 2 -E^2 Flux\n";
         cin>>spec_type;
         if(spec_type)  cout<<"Units: 0  -MeV/(cm2 s sr) vs. MeV\n"
                            <<"       1  -GeV/(m2  s sr) vs. GeV\n";
         else  cout<<"Units: 0  -1/(cm2 s sr MeV) vs. MeV\n"
                   <<"       1  -1/(m2  s sr GeV) vs. GeV\n";
         cin>>spec_units;
         break;

      case 9:
      default:
         cout<<"\nYour choise: Exit now\n";
	 exit(0);
   }
   if(key<3)
   {
      cout<<"experimental: give an option (1,2,3)\n"
          <<"              a fraction (1.2),\n"
          <<"              O16 norm. energy [150 MeV/nucleon],\n" 
          <<"              and energy scale [450 MeV/nucleon]\n";
      cin>>add_nuclei_case>>add_nuclei_fraction>>add_nuclei_energy>>add_nuclei_escale;
   }

// read fits files and prepare output files
   for(i=0; i<argc-1; i++)
   {  
      fout = fopen(argv[i+1],"w");
      strcpy(tmp,"../FITS/nuclei_");
//      readheader(strcat(tmp,argv[i+1])); exit(0);
      readimage(strcat(tmp,argv[i+1]),fout,key);
      fclose(fout);
   }

// make a gnuplot driver file
   make_gnufile(argc,argv,key);

   cout<<"Done...\n";
   exit(0);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void readimage( char* filename, FILE *fout, int key)
{
   fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
   int i,j,k,m,n, ir,ip,iz,status=0, iprotons=0,nfound, anynull;
   long naxis, *naxes,  *Z, *A, *K, felement=1, nelements;
   float *array, tmp,tmp1,nullval=0.,E_2, case2[10]={0.},fraction;
   double *crval, *cdelt, T,T0, Mp=939., // MeV, nucleon rest mass
          Mee=0.511, *numerator, *denominator, *modnum, *moddenom, sum;
   double abund_elem[30];
   char comment[FLEN_COMMENT], card[FLEN_CARD], *incl[]={"NUCZ0??","NUCA0??","NUCK0??"};

   if ( fits_open_file(&fptr, filename, READONLY, &status) )
      printerror( status );

// read the NAXIS keyword to get number of axis
   if ( fits_read_key(fptr, TLONG, "NAXIS", &naxis, comment, &status) )
      printerror( status );

   naxes= new long[naxis];
   crval= new double[naxis];
   cdelt= new double[naxis];

// read the NAXES keywords to get dimensions
   if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxes, &nfound, &status) )
      printerror( status );
// read the CRVAL keywords to get starting points
   if ( fits_read_keys_dbl(fptr, "CRVAL", 1, naxis, crval, &nfound, &status) )
      printerror( status );
// read the CDELT keywords to get increments
   if ( fits_read_keys_dbl(fptr, "CDELT", 1, naxis, cdelt, &nfound, &status) )
      printerror( status );

//   for(i=0;i<naxis;i++) printf(" %d ",naxes[i]); printf("\n");

   Z= new long[naxes[3]];
   A= new long[naxes[3]];
   K= new long[naxes[3]]; for(i=0; i<naxes[3]; K[i++]=0);
// read the Z keywords to get nuclear charges
   if ( fits_read_keys_lng(fptr, "NUCZ00", 1, 9, Z, &nfound, &status) )
      printerror( status );
   if ( fits_read_keys_lng(fptr, "NUCZ0", 10, naxes[3]-9, &Z[9], &nfound, &status) )
      printerror( status );
// read the A keywords to get atomic numbers
   if ( fits_read_keys_lng(fptr, "NUCA00", 1, 9, A, &nfound, &status) )
      printerror( status );
   if ( fits_read_keys_lng(fptr, "NUCA0", 10, naxes[3]-9, &A[9], &nfound, &status) )
      printerror( status );
// read the K keywords to get the number of K-electrons (=0,1)
   if ( fits_read_keys_lng(fptr, "NUCK00", 1, 9, K, &nfound, &status) )
      printerror( status ); i=nfound;
   if ( fits_read_keys_lng(fptr, "NUCK0", 10, naxes[3]-9, &K[9], &nfound, &status) )
      printerror( status );
   if(i+nfound==0) for(i=0; i<naxes[3]; K[i++]=0);
//   printf("%2d.%2d%5d%5d%5d\n",Z[1],A[1],K[1],i+nfound, status);

   for(nelements=1, i=0; i<naxis; i++) nelements*= naxes[i];
   array = new float[nelements];
   if ( fits_read_img(fptr, TFLOAT, felement, nelements, &nullval,
      array, &anynull, &status) ) printerror( status );
   if ( fits_close_file(fptr, &status) ) printerror( status );

   ir = (int) ((R-crval[0])/cdelt[0]);
   for(k=0; k<naxes[3]; k++) if(101 <= 100*Z[k]+A[k]) break; // find nuclei
   for(n=1, i=0; i<naxis-1; i++) n*= naxes[i];
   for(i=0; i<Z[naxes[3]-1]; i++) abund_elem[i]=0.;

// experimental
   if(key<3)
   {
      add_nuclei(0,0,0,0.,0.);
      for(i=0; i<naxes[3]; i++) if(816 == 100*Z[i]+A[i]) break; // find O16
      ip = (int) ((log10(add_nuclei_energy)-crval[2])/cdelt[2]);
      fraction= i==naxes[3]? 0.: array[i*n+ip*naxes[0]+ir]*add_nuclei_fraction
         *pow(add_nuclei_energy,-2)/add_nuclei(add_nuclei_case,Z[i],A[i],add_nuclei_energy,add_nuclei_escale);
   } //printf("#>> %3d  %.3e\n",i,fraction);
//

// DERIVE ISOTOPIC ABUNDANCES (key = 0)

   if(key==0)
   {
      for(sum=0., m=0, i=k; i<naxes[3]; i++)
      {
         T0 = Ekin+Z[i]*Phi/A[i];
         ip = (int) ((log10(T0)-crval[2])/cdelt[2]);
 	 if(ip<0)
         {
	    cout<<"\nThe kinetic energy entered is too low, the lowest possible is "
	        <<pow(10,crval[2])-Z[i]*Phi/A[i]<<" MeV/nucleon \nexit\n";
            exit(1);
         }
	 tmp  = array[i*n+ ip   *naxes[0]+ir];             // interpolation in R @ p1
         tmp +=(array[i*n+ ip   *naxes[0]+ir+1]-tmp )/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
         tmp *= pow(10,-2* (crval[2]+ ip   *cdelt[2]));    // *E1^-2
         tmp1 = array[i*n+(ip+1)*naxes[0]+ir];             // interpolation in R @ p2
         tmp1+=(array[i*n+(ip+1)*naxes[0]+ir+1]-tmp1)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
         tmp1*= pow(10,-2* (crval[2]+(ip+1)*cdelt[2]));    // *E2^-2
         if(tmp<=0) { cout<<Z[i]<<"."<<A[i]<<": negative value in spectrum ! "<<tmp<<endl;  tmp=0.; }
// Modulation: Gleeson & Axford (1968) force-field approximation
	 else 
// experimental
//            tmp  = pow(10., log10(tmp)
//               +(log10(T0)-crval[2]-ip*cdelt[2]) *log10(tmp1/tmp)/cdelt[2] )
//               *Ekin*(Ekin+2*Mp)/T0/(T0+2*Mp);
            tmp  = (fraction*add_nuclei(add_nuclei_case,Z[i],A[i],T0,add_nuclei_escale)+pow(10., log10(tmp)
	       +(log10(T0)-crval[2]-ip*cdelt[2]) *log10(tmp1/tmp)/cdelt[2] ) )
               *Ekin*(Ekin+2*Mp)/T0/(T0+2*Mp);

// Normalization to hydrogen=10^6
         if(i==k) { if(101 == 100*Z[i]+A[i]) T=1.e-6*tmp; else T = 1; }

//cout<<" tmp= "<<tmp<<" T= "<<T<<endl;
         j=i;
         if(datakey == 2) if(z[0] == Z[j]) { a[m]=A[j]; case2[m++]+=tmp/T; if(K[j]) m--; }
         if(datakey == 0)
	 {
   	    if(i==k)
            {
               sum = tmp;
               continue;
            }
	    if(i>k)
               if(100*Z[i-1]+A[i-1]==100*Z[i]+A[i]) 
               {
		  sum+= tmp;
                  if(i<naxes[3]-1) continue;
               }
            j--;
         }
         else
	 {
   	    if(i==k)
            {
               sum = tmp;
               continue;
            }
	    if(i>k)
               if(Z[i-1]==Z[i]) 
               {
		  sum+= tmp;
                  if(i<naxes[3]-1) continue;
               }
            j--;
         }

// Normalization to Si28 = 10^3
         if(1429 == 100*Z[i]+A[i]) norm_nuc=sum/T;

         if(datakey == 0) fprintf(fout,"  %3d  %10.3e  %3d %2d  %3d\n",Z[j],sum/T, Z[j],A[j],j+1);
         else abund_elem[Z[j]-1]=sum/T;
         sum = tmp;
      }
      T = 1.e-3 * abund_elem[7];                 // normalization to Oxigen
      if(Ekin < 400) T = 1.e-2 * abund_elem[13]; // normalization to Silicon
      if(datakey == 1)
	for(j=2; j<Z[naxes[3]-1]; j++)
            if(abund_elem[j]) 
               fprintf(fout,"  %3d  %10.3e  %3d %2d  %3d\n",j+1,abund_elem[j]/T,j+1,0.,0.);
      if(datakey == 2) 
	 for(i=0; i<m; i++) 
            fprintf(fout,"  %3d  %10.3e  %3d %2d  %3d\n",a[i],case2[i]/abund_elem[z[0]-1],z[0],a[i],0);
   }


// DERIVE PARTICLE SPECTRA

   if(key==1)
   {
      numerator  = new double[naxes[2]];   // interstellar spectra
      modnum     = new double[naxes[2]];   // modulated
      denominator= new double[naxes[2]];   // secondary p & tertiary pbars
      moddenom   = new double[naxes[2]];   // modulated secondary p & tertiary pbars
      for(i=0; i<naxes[2]; i++) numerator[i]=denominator[i]=modnum[i]=moddenom[i]=0.;
      iz = 0;
//      for(i=0;i<naxes[3]; i++) printf("%2d.%2d\n",Z[i],A[i]);
      for(j=0; z[j]!=0; j++)
//      for(j=0; z[j]!=0 || a[j]!=0; j++)
      {
	 if(z[j]!=0 && a[j]==0)  // sum of all isotopes of an element
	 {
            for(k=iz; k<naxes[3]; k++) if(z[j] == Z[k] && 0!=A[k]) break;
            z[j+1] = z[j];
            a[j+1] = 0;
            iz = k+1;
         }
	 else for(k=iz; k<naxes[3]; k++) if(100*z[j]+a[j] == 100*Z[k]+A[k]) break;
         if(K[k]==1&&A[k]>1)
	 {
	    iz = k+1; 
            j--;
         }  
	 if(iprotons==1)  // to account for prim+sec protons and sec+tert pbars
	 {
            for(i=0; i<naxes[2]; i++) denominator[i]=numerator[i];
            for(i=0; i<naxes[2]; i++) moddenom[i]   =modnum[i];
         }
         if(abs(z[j])==1) iprotons++;
         if(k==naxes[3])
	 {
            if(a[j]==0) 
            {
	       z[j] = 0;
               goto label1;
            }
	    printf(" nucleus (%d,%d) is unstable or was not included in the run: exit\n",z[j],a[j]);
            exit(1);
         }
         if(abs(z[j])!=1) a[j] = A[k];  // all nuclei except protons and pbars
         for(i=0; i<naxes[2]; i++)
         {
            T  = pow(10., crval[2]+i*cdelt[2]);
            tmp = array[k*n+i*naxes[0]+ir];         // interpolation in R
            tmp+=(array[k*n+i*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
	    numerator[i]+= tmp*pow(T,spec_type-2);
            numerator[i]+= fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T,add_nuclei_escale)*pow(T,spec_type);
// Modulation: Gleeson & Axford (1968) force-field approximation
            T0 = T+abs(Z[k])*Phi/A[k];
            ip = (int) ((log10(T0)-crval[2])/cdelt[2]);
            if(ip>=naxes[2]-1)
	    { 
	       modnum[i] = numerator[i];
               continue;
            }
            tmp  = array[k*n+ ip   *naxes[0]+ir];   // interpolation in R
            tmp +=(array[k*n+ ip   *naxes[0]+ir+1]-tmp )/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp *= pow(10,-2* (crval[2]+ ip   *cdelt[2]));
            tmp1 = array[k*n+(ip+1)*naxes[0]+ir];   // interpolation in R
            tmp1+=(array[k*n+(ip+1)*naxes[0]+ir+1]-tmp1)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp1*= pow(10,-2* (crval[2]+(ip+1)*cdelt[2]));
	    if(tmp>0.&&tmp1>0.) modnum[i]+= ( fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T0,add_nuclei_escale) +pow(10., // interp. in E and modulation
	       log10(tmp) +(log10(T0)-crval[2]-ip*cdelt[2])/cdelt[2] *log10(tmp1/tmp) ) )
               *T*(T+2*Mp)/T0/(T0+2*Mp) *pow(T,spec_type);
         }
      }
label1:
      for(i=0; i<naxes[2]; i++)
      {
	 T= pow(10.,crval[2]+i*cdelt[2]);
	 if(spec_units) 
	 {
	    T             *= 1.e-3;
	    if(spec_type)      // transfer MeV/cm^2 to GeV/m^2
	    {
	       numerator[i]  *= 10.;   modnum[i]     *= 10.;
	       denominator[i]*= 10.;   moddenom[i]   *= 10.;
            } else             // transfer [cm^2 MeV]^-1 to [m^2 GeV]^-1
	    {
	       numerator[i]  *= 1.e7;  modnum[i]     *= 1.e7;
	       denominator[i]*= 1.e7;  moddenom[i]   *= 1.e7;
            }
         }
	 ymax =  ymax < numerator[i] ? numerator[i]: ymax;
	 fprintf(fout,"  %12.4g   %12.4g   %12.4g   %12.4g   %12.4g\n",T,
            numerator[i], modnum[i], denominator[i], moddenom[i]);
      }
      if(ymax == 0.) ymax=1.;
   }

// DERIVE ISOTOPIC RATIO

   if(key==2)
   {
      numerator  = new double[naxes[2]];   // interstellar spectra
      denominator= new double[naxes[2]];   // interstellar spectra
      modnum     = new double[naxes[2]];   // modulated
      moddenom   = new double[naxes[2]];   // modulated
      for(i=0; i<naxes[2]; i++) 
         numerator[i] =denominator[i] =modnum[i] = moddenom[i] = 0.;
      for(iz=0, j=0; z[j]!=0 || a[j]!=0; j++)
      {
	 for(k=iz; k<naxes[3]; k++) if(100*z[j]+a[j] == 100*Z[k]+A[k]) break;
         if(K[k])
	 {
	    iz = k+1; 
            j--;
         }  
	 if(z[j]==-1 && a[j]==1) iz = k+1; // antiprotons: find secondary & tertiary
         if(k==naxes[3])
	 {
	    printf(" nucleus (%d,%d) is unstable or was not included in the run: exit\n",z[j],a[j]);
            continue;
         }
//cout<<endl;
//cout<<Z[k]<<" "<<A[k]<<endl;
         for(i=0; i<naxes[2]; i++)
         {
            tmp = array[k*n+i*naxes[0]+ir];
            tmp+=(array[k*n+i*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
	    numerator[i]+= tmp;
//cout<<numerator[i]<<" ";
            T  = pow(10., crval[2]+i*cdelt[2]);
            numerator[i]+= fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T,add_nuclei_escale)*pow(T,2);
// Modulation: Gleeson & Axford (1968) force-field approximation
            T0 = T+abs(Z[k])*Phi/A[k];
            ip = (int) ((log10(T0)-crval[2])/cdelt[2]);
            if(ip>=naxes[2]-1)
	    { 
	       modnum[i] = numerator[i];
               continue;
            }
	    E_2  = pow(10,-2* (crval[2]+ip*cdelt[2]));
            tmp  = array[k*n+ip*naxes[0]+ir];
            tmp +=(array[k*n+ip*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp *= E_2;
	    E_2  = pow(10,-2* (crval[2]+(ip+1)*cdelt[2]));
            tmp1 = array[k*n+(ip+1)*naxes[0]+ir];
            tmp1+=(array[k*n+(ip+1)*naxes[0]+ir+1]-tmp1)/cdelt[0]
               *(R-crval[0]-ir*cdelt[0]);
            tmp1*= E_2;
	    modnum[i]+= ( fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T0,add_nuclei_escale)+ pow(10., log10(tmp)
	       +(log10(T0)-crval[2]-ip*cdelt[2])/cdelt[2] *log10(tmp1/tmp) ) )
               *T*(T+2*Mp)/T0/(T0+2*Mp);
         }
      }
//      for(i=0; i<naxes[2]; i++) {T= pow(10.,crval[2]+i*cdelt[2])/1000.; cout<<T<<"  "<<numerator[i]<<endl;}

      if(-1 ==z[0]) Phi=Phi1;  // pbar/p

      for(iz=0, j=0; z1[j]!=0 || a1[j]!=0; j++)
      {
    	 for(k=iz; k<naxes[3]; k++) if(100*z1[j]+a1[j] == 100*Z[k]+A[k]) break;
         if(K[k])
	 {
	    iz = k+1; 
            j--;
         }  
	 if(z1[j]==1 && a1[j]==1) iz = k+1; // protons: find primary and secondary
         if(k==naxes[3])
	 {
	    printf(" nucleus (%d,%d) is unstable or was not included in the run: exit\n",z1[j],a1[j]);
            continue;
         }
//cout<<endl;
//cout<<Z[k]<<" "<<A[k]<<endl;
         for(i=0; i<naxes[2]; i++)
         {
            tmp = array[k*n+i*naxes[0]+ir];
            tmp+=(array[k*n+i*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
	    denominator[i]+= tmp;
//cout<<denominator[i]<<" ";
            T  = pow(10., crval[2]+i*cdelt[2]);
            denominator[i]+= fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T,add_nuclei_escale)*pow(T,2);
// Modulation: Gleeson & Axford (1968) force-field approximation
            T0 = T+abs(Z[k])*Phi/A[k];
            ip = (int) ((log10(T0)-crval[2])/cdelt[2]);
            if(ip>=naxes[2]-1)
	    {
	       moddenom[i] = denominator[i];
               continue;
            }
	    E_2  = pow(10,-2* (crval[2]+ip*cdelt[2]));
            tmp  = array[k*n+ip*naxes[0]+ir];
            tmp +=(array[k*n+ip*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp *= E_2;
	    E_2  = pow(10,-2* (crval[2]+(ip+1)*cdelt[2]));
            tmp1 = array[k*n+(ip+1)*naxes[0]+ir];
            tmp1+=(array[k*n+(ip+1)*naxes[0]+ir+1]-tmp1)/cdelt[0]
               *(R-crval[0]-ir*cdelt[0]);
            tmp1*= E_2;
            moddenom[i]+= ( fraction*add_nuclei(add_nuclei_case,Z[k],A[k],T0,add_nuclei_escale)+ pow(10, log10(tmp)
	       +(log10(T0)-crval[2]-ip*cdelt[2])/cdelt[2] *log10(tmp1/tmp) ) )
               *T*(T+2*Mp)/T0/(T0+2*Mp);
         }
      }
//      for(i=0; i<naxes[2]; i++) {T= pow(10.,crval[2]+i*cdelt[2])/1000.; cout<<T<<"  "<<denominator[i]<<endl;}

      for(i=0; i<naxes[2]; i++)
      { 
	 T= pow(10.,crval[2]+i*cdelt[2]);
	 if(spec_units) T *= 1.e-3; // transfer MeV -> GeV
         if(!denominator[i]) continue;
	 ymax = (ymax<modnum[i]/moddenom[i]) ? modnum[i]/moddenom[i]: ymax;
	 fprintf(fout,"  %.3g   %.3g   %.3g\n",T,
            numerator[i]/denominator[i], modnum[i]/moddenom[i]);
      }
      if(ymax <= 0.) ymax=1;
   }

// DERIVE SPECTRA OF ELECTRONS AND POSITRONS

   if(key==3)
   {
      numerator  = new double[naxes[2]];   // interstellar spectra
      modnum     = new double[naxes[2]];   // modulated
      denominator= new double[naxes[2]];   // secondary e-
      moddenom   = new double[naxes[2]];   // modulated secondary e-
      for(i=0; i<naxes[2]; i++) numerator[i]=denominator[i]=modnum[i]=moddenom[i]=0.;
      iz = 0;
      for(j=0; z[j]!=0; j++)
      {
         for(k=iz; k<naxes[3]; k++) if(100*z[j]+a[j] == 100*Z[k]+A[k]) break;
         z[j+1] = z[j];
         iz = k+1;
	 if(abs(z[j])==1 && iprotons==1)  // to account for prim+sec e- and sec e+
	 {
            for(i=0; i<naxes[2]; i++) denominator[i]=numerator[i];
            for(i=0; i<naxes[2]; i++) moddenom[i]   =modnum[i];
         }
         iprotons++;
         if(k==naxes[3]) goto label2;
         for(i=0; i<naxes[2]; i++)
         {
            T  = pow(10., crval[2]+i*cdelt[2]);
            tmp = array[k*n+i*naxes[0]+ir];
            tmp+=(array[k*n+i*naxes[0]+ir+1]-tmp)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
	    numerator[i]+= tmp *pow(T,spec_type-2);
// Modulation: Gleeson & Axford (1968) force-field approximation
            T0 = T+abs(Z[k])*Phi;
            ip = (int) ((log10(T0)-crval[2])/cdelt[2]);
            if(ip>=naxes[2]-1)
	    { 
	       modnum[i] = numerator[i];
               continue;
            }
            tmp  = array[k*n+ ip   *naxes[0]+ir];
            tmp +=(array[k*n+ ip   *naxes[0]+ir+1]-tmp )/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp *= pow(10,-2* (crval[2]+ ip   *cdelt[2]));
            tmp1 = array[k*n+(ip+1)*naxes[0]+ir];
            tmp1+=(array[k*n+(ip+1)*naxes[0]+ir+1]-tmp1)/cdelt[0]*(R-crval[0]-ir*cdelt[0]);
            tmp1*= pow(10,-2* (crval[2]+(ip+1)*cdelt[2]));
	    if(tmp>0.) modnum[i]+= pow(10., 
               log10(tmp) +(log10(T0)-crval[2]-ip*cdelt[2])/cdelt[2] *log10(tmp1/tmp) )
               *T*(T+2*Mee)/T0/(T0+2*Mee) *pow(T,spec_type);
         }
      }
label2:
      for(i=0; i<naxes[2]; i++)
      {
	 T= pow(10.,crval[2]+i*cdelt[2]);
	 if(spec_units) 
	 {
	    T             *= 1.e-3;
	    if(spec_type)      // transfer MeV/cm^2 to GeV/m^2
	    {
	       numerator[i]  *= 10.;   modnum[i]     *= 10.;
	       denominator[i]*= 10.;   moddenom[i]   *= 10.;
            } else             // transfer [cm^2 MeV]^-1 to [m^2 GeV]^-1
	    {
	       numerator[i]  *= 1.e7;  modnum[i]     *= 1.e7;
	       denominator[i]*= 1.e7;  moddenom[i]   *= 1.e7;
            }
         }
	 ymax =  (ymax<modnum[i]) ? modnum[i]: ymax;
	 fprintf(fout,"  %12.3g   %12.3g   %12.3g   %12.3g   %12.3g\n",T,
            numerator[i], modnum[i], denominator[i], moddenom[i]);
      }
      if(ymax == 0.) ymax=1.;
   }

   delete naxes;
   delete crval;
   delete cdelt;
   delete Z;
   delete A;
   delete K;
   delete array;
   if(key !=0)
   {
      delete numerator;
      delete modnum;
      delete denominator;
      delete moddenom;
   }
   cout<<"Arrays deleted"<<endl;
   return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void make_gnufile(int argc, char*argv[], int key)
{
   FILE *fgnu,*fin;
   int i,j,Z,A,A1=0;
   float llimit=1., ulimit=1., yscale, abund;
   char s[2], gev[]="_gev";
   fgnu = fopen("tmp.gnu","w");

   if(!spec_units) gev[0]='\0';

   fprintf(fgnu,"set terminal post port\n");
   fprintf(fgnu,"set output \"tmp.ps\"\n");
   fprintf(fgnu,"set nogrid\n");
//   fprintf(fgnu,"set key title \"models     \"\n");
   fprintf(fgnu,"set tics in\n");
   fprintf(fgnu,"set ticslevel 0.5\n");
//   fprintf(fgnu,"set ticscale 1 0.5\n");
//   fprintf(fgnu,"set x2label \"\" 0.000000,0.000000\n");
//   fprintf(fgnu,"set x2range [*: *] noreverse nowriteback  # (currently [-10:10] )\n");
//   fprintf(fgnu,"set y2label \"\" 0.000000,0.000000  ""\n");
//   fprintf(fgnu,"set y2range [*:*] noreverse nowriteback  # (currently [-10:10] )\n");
   fprintf(fgnu,"set zero 1e-15\n");

// PLOT ISOTOPIC ABUNDANCES
   if(key==0)
   {
      fprintf(fgnu,"set logscale y 10\n");
      fprintf(fgnu,"set size 1,1\n");
      if (datakey == 1) fprintf(fgnu,"set title \"ELEMENTAL ABUNDANCES @ %.0f MeV/nucleon\" 0.,0.\n",Ekin);
      else fprintf(fgnu,"set title \"ISOTOPIC ABUNDANCES @ %.0f MeV/nucleon\" 0.,0.\n",Ekin);
      fprintf(fgnu,"set xlabel \"Nucleus charge, Z\" 0.,0.\n");
      if (!datakey) fprintf(fgnu,"set ylabel \"Relative abundances (H=%.0g)\" 0.,0.\n",1.e6);
      if(datakey == 1)
      {
         if (fabs(Ekin-0.18e3)<0.07e3) fprintf(fgnu,"set ylabel \"Relative abundances (Si=%.0f)\" 0.,0.\n",1.e2);
         else fprintf(fgnu,"set ylabel \"Relative abundances (O=%.0f)\" 0.,0.\n",1.e3);
      }
      fin = fopen(argv[1],"r");
      while(fscanf(fin,"%d%f%d%d%d",&Z,&abund,&i,&A,&j)!=EOF)
      {
	 if(!A1) A1 = A;
	 if(abund>3.e4 || abund <= 0.) continue;
	 if(datakey) { if(abund<llimit) llimit=abund/3.; }
         else llimit= 0.1;
	 if(abund>ulimit) ulimit=abund*1.3;
         if(abund<llimit) continue;
         if (!datakey) fprintf(fgnu,"set label \"%d\" at %f,%g\n",A,Z+.1,abund*1.1);
      }
      if(datakey==0) { llimit=0.3; ulimit=1.e4; }
      if(datakey==1 && llimit<0.2) llimit=0.2;
      fprintf(fgnu,"set label \"Phi = %.0f MV\" at %.1f,%.3f\n",Phi,5.,llimit*3.);
      if (datakey == 2)
      {
         fprintf(fgnu,"set nologscale y\n");
         fprintf(fgnu,"set size 0.7,0.5\n");
         fprintf(fgnu,"set xlabel \"Atomic number, A\" 0.,0.\n");
         fprintf(fgnu,"set ylabel \"Relative abundances\" 0.,0.\n");
         fprintf(fgnu,"set nokey\n"); 
         fprintf(fgnu,"set label \"Z = %d\" at %.1f,%.1f\n",z[0],A1+0.2,0.9);
         fprintf(fgnu,"set xtics %d,%d\n",A1,1);
         fprintf(fgnu,"plot [%f:%f][0:1]\\\n",A1-0.1,A+0.5);
      }
      else fprintf(fgnu,"plot [2.1:29][%.2f:%.0f]\\\n",llimit,ulimit);
      for(j=1, i=0; i<argc-1; i++)
      {
 	 if(i==2) j++;
	 if(i!=0) fprintf(fgnu,"\\\n, ");
         if (datakey == 2) fprintf(fgnu,"\"%s\" w p %d,\"%s\" w i %d",argv[i+1],2*i+7,argv[i+1],3);
         else fprintf(fgnu,"\"%s\" w p %d",argv[i+1],j++);
         if (datakey == 1) fprintf(fgnu,",\"%s\" w l %d",argv[i+1],j);
      }

      if (datakey == 0) 
         if(fabs(Ekin-0.20e3)<0.03e3)
	 {
	   fprintf(fgnu,"\\\n, \"ACE_isot_Wiedenbeck.dat\" u ($1-.2):($3/1000.*%f):($4/1000.*%f) w e %d",norm_nuc,norm_nuc,1);
           fprintf(fgnu,"\\\n, \"ACE_isot_Wiedenbeck.dat\" u ($1-.2):($3/1000.*%f) w p %d",norm_nuc,72);
         }

      if (datakey == 1)
      { 
         fprintf(fgnu,"\\\n");
         if(fabs(Ekin-0.20e3)<0.2e3) fprintf(fgnu,",\"ACE_elem_Wiedenbeck.dat\"  u ($1-.4):2:3 w e %d",-1);
         if(fabs(Ekin-0.18e3)<0.07e3) fprintf(fgnu,",\"abund_ulysses96.dat\"   u 1: 2: 3 w e %d",-1);
         if(fabs(Ekin-0.62e3)<0.08e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1: 2: 3 w e %d",-1);
         if(fabs(Ekin-0.80e3)<0.10e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1: 4: 5 w e %d",-1);
         if(fabs(Ekin-1.00e3)<0.10e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1: 6: 7 w e %d",-1);
         if(fabs(Ekin-1.25e3)<0.15e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1: 8: 9 w e %d",-1);
         if(fabs(Ekin-1.60e3)<0.20e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:10:11 w e %d",-1);
         if(fabs(Ekin-2.10e3)<0.25e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:12:13 w e %d",-1);
         if(fabs(Ekin-2.65e3)<0.30e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:14:15 w e %d",-1);
         if(fabs(Ekin-3.35e3)<0.40e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:16:17 w e %d",-1);
         if(fabs(Ekin-4.30e3)<0.55e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:18:19 w e %d",-1);
         if(fabs(Ekin-5.60e3)<0.75e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:20:21 w e %d",-1);
         if(fabs(Ekin-7.50e3)<1.25e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:22:23 w e %d",-1);
         if(fabs(Ekin-10.6e3)<1.85e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:24:25 w e %d",-1);
         if(fabs(Ekin-16.2e3)<3.00e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:26:27 w e %d",-1);
         if(fabs(Ekin-35.0e3)<5.00e3) fprintf(fgnu,",\"abund_engelmann90.dat\" u 1:28:29 w e %d",-1);
      }
      if (datakey == 2)
      {
	 if(z[0] ==  3) // Li
	 { 
            fprintf(fgnu,",\"Li_iso.dat\" u 1: 2: 3 w e %d",4);
         }
	 if(z[0] ==  4) // Be
	 { 
            fprintf(fgnu,",\"Be_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Be_iso.dat\" u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] ==  5) // B
	 { 
            fprintf(fgnu,",\"B_iso.dat\"  u 1: 2: 3 w e %d",4);
         }
	 if(z[0] ==  6) // C
	 { 
            fprintf(fgnu,",\"C_iso.dat\"  u 1: 2: 3 w e %d,\\\n\"C_iso.dat\"  u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] ==  7) // N
	 { 
            fprintf(fgnu,",\"N_iso.dat\"  u 1: 2: 3 w e %d,\\\n\"N_iso.dat\"  u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] ==  8) // O
	 { 
            fprintf(fgnu,",\"O_iso.dat\"  u 1: 2: 3 w e %d,\\\n\"O_iso.dat\"  u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] == 10) // Ne
	 { 
            fprintf(fgnu,",\"Ne_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Ne_iso.dat\" u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] == 12) // Mg
	 { 
            fprintf(fgnu,",\"Mg_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Mg_iso.dat\" u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] == 13) // Al
	 { 
            fprintf(fgnu,",\"Al_iso.dat\"  u 1: 2: 3 w e %d",4);
         }
	 if(z[0] == 14) // Si
	 { 
            fprintf(fgnu,",\"Si_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Si_iso.dat\" u 4: 5: 6 w e %d\\\n",4,6);
            fprintf(fgnu,",\"Si_iso.dat\" u 7: 8: 9 w e %d",8);
         }
	 if(z[0] == 16) // S
	 { 
            fprintf(fgnu,",\"S_iso.dat\"  u 1: 2: 3 w e %d,\\\n\"S_iso.dat\"  u 4: 5: 6 w e %d",4,6); 
         }
	 if(z[0] == 17) // Cl
	 { 
            fprintf(fgnu,",\"Cl_iso.dat\"  u 1: 2: 3 w e %d",4);
         }
	 if(z[0] == 20) // Ca
            fprintf(fgnu,",\"Ca_iso.dat\" u 1: 2: 3 w e %d",4);
	 if(z[0] == 22) // Ti
            fprintf(fgnu,",\"Ti_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Ti_iso.dat\" u 4: 5: 6 w e %d\\\n",4,6);
	 if(z[0] == 23) // V
            fprintf(fgnu,",\"V_iso.dat\"  u 1: 2: 3 w e %d,\\\n\"V_iso.dat\" u 4: 5: 6 w e %d\\\n",4,6);
	 if(z[0] == 24) // Cr
            fprintf(fgnu,",\"Cr_iso.dat\" u 1: 2: 3 w e %d",4);
	 if(z[0] == 25) // Mn
	 { 
            fprintf(fgnu,",\"Mn_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Mn_iso.dat\" u 4: 5: 6 w e %d",4,6);
         }
	 if(z[0] == 26) // Fe
	 { 
            fprintf(fgnu,",\"Fe_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Fe_iso.dat\" u 4: 5: 6 w e %d\\\n",4,6);
            fprintf(fgnu,",\"Fe_iso.dat\" u 7: 8: 9 w e %d",8);
         }
	 if(z[0] == 27) // Co
            fprintf(fgnu,",\"Co_iso.dat\" u 1: 2: 3 w e %d",4);
	 if(z[0] == 28) // Ni
            fprintf(fgnu,",\"Ni_iso.dat\" u 1: 2: 3 w e %d,\\\n\"Ni_iso.dat\" u 4: 5: 6 w e %d",4,6);
      }
      fprintf(fgnu,"\n");
      fclose(fin);
   }

// PLOT PARTICLE SPECTRA

   if(key==1)
   {
      llimit=90; ulimit = 5.e6;
      if(!spec_type)  ulimit/=10.;
      if(z[0] == -1)  { llimit =10.;  ulimit =1.e5; }
      if(z[0] ==  1 || z[0] ==  2)    ulimit =3.e8;
      if(spec_units)  { llimit /=1000.; ulimit /=1000.; }
      fprintf(fgnu,"set logscale xy 10\n");
      fprintf(fgnu,"set size 1,0.7\n");
      fprintf(fgnu,"set title \"PARTICLE SPECTRA: ");
      for(i=0; z[i]!=0 || a[i]!=0; i++)
      {
	 if(i>0) fprintf(fgnu,"+");
         if(abs(z[0])==1) a[i]=1;  // protons & antiprotons 
         fprintf(fgnu,"%2d.%2d", z[i],a[i]);
         if(abs(z[0])==1) break;   // protons & antiprotons 
      }
      fprintf(fgnu,"\n");
      yscale = pow(10,(int) log10(ymax*1.e-3));
      if(z[0] ==  1 || z[0] ==  2) yscale/=10.;
      if(spec_type)                // E^2 Flux
	 if(spec_units)            // GeV/m^2
         { 
            fprintf(fgnu,"set xlabel \"Kinetic energy, GeV/nucleon\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"E^2 Flux, GeV/nucleon m^-2 s^-1 sr^-1\" 0.,0.\n");
         } else                    // MeV/cm^2
         { 
            fprintf(fgnu,"set xlabel \"Kinetic energy, MeV/nucleon\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"E^2 Flux, MeV/nucleon cm^-2 s^-1 sr^-1\" 0.,0.\n");
         }
      else                         // Flux
	 if(spec_units)            // 1/(m^2 GeV)
         {
            fprintf(fgnu,"set xlabel \"Kinetic energy, GeV/nucleon\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"Flux, m^-2 s^-1 sr^-1 (GeV/nucleon)^-1\" 0.,0.\n");
         } else                    // 1/(cm^2 MeV)
         {
            fprintf(fgnu,"set xlabel \"Kinetic energy, MeV/nucleon\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"Flux, cm^-2 s^-1 sr^-1 (MeV/nucleon)^-1\" 0.,0.\n");
         }
//      if(z[0]==-1)  yscale = pow(10,(int) log10(ymax*1.e-3));  // Antiprotons
// Protons & He
      if((z[0]==1 || z[0]==2) && !spec_type && spec_units) { yscale = 1.e-3; fprintf(fgnu,"set size 0.9,1\n"); }
      if(z[0]==26)  yscale = pow(10,(int) log10(ymax*1.e-2));   // Iron
      fprintf(fgnu,"set label \"Phi = %.0f MV\" at %.1g,%.1g\n",Phi,10*llimit, 10*yscale);

      if(datakey) fprintf(fgnu,"set nokey\n");
      fin = fopen(argv[1],"r");
      fprintf(fgnu,"plot [%.1g:%.1g][%.1g:%.2g] ",llimit, ulimit, yscale, ymax*1.2 ); 
      for(i=0; i<argc-1; i++)
      {
	 if(i!=0) fprintf(fgnu,", ");
         fprintf(fgnu,"\"%s\" w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:3 w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:4 w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:5 w l %d ",argv[i+1],i+1);
      }
      if(datakey ==-1)  // antiprotons
      {
	 if(spec_type)
	 {
            fprintf(fgnu,", \"pbar_stochaj00%s.dat\" w l %d, \"pbar_stochaj00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"pbar_orito00%s.dat\" w l %d, \"pbar_orito00%s.dat\" u 3:4 w p %d",gev,1,gev,7);
         } else
         {
            fprintf(fgnu,", \"pbar0_stochaj00.dat\" w l %d, \"pbar0_stochaj00.dat\" u 3:4 w p %d",1,5);
            fprintf(fgnu,", \"pbar0_orito00.dat\" w l %d, \"pbar0_orito00.dat\" u 3:4 w p %d",1,7);
            fprintf(fgnu,", \"pbar0_maeno01.dat\" w l %d, \"pbar0_maeno01.dat\" u 3:4 w p %d",1,6);
            fprintf(fgnu,", \"pbar0_boezio01.dat\" w l %d, \"pbar0_boezio01.dat\" u 3:4 w p %d",1,8);
         }
      }
      if(datakey == 1)  // protons
      {
	 if(spec_type)
	 {
            fprintf(fgnu,", \"H_boezio99%s.dat\" w l %d, \"H_boezio99%s.dat\" u 3:4 w p %d",gev,1,gev,2);
            fprintf(fgnu,", \"H_menn00%s.dat\" w l %d,   \"H_menn00%s.dat\" u 3:4 w p %d",gev,1,gev,4);
//            fprintf(fgnu,", f(x,1) w l %d, f(x,-1) w l %d",-1,-1);            //Menn interstellar
            fprintf(fgnu,", \"H_sanuki00%s.dat\" w l %d, \"H_sanuki00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"H_alcaraz00a%s.dat\" w l %d, \"H_alcaraz00a%s.dat\" u 3:4 w p %d",gev,1,gev,7);
            fprintf(fgnu,", \"sokol%s.dat\" w e %d",gev,6);
            fprintf(fgnu,", 1.6e4*x**(-0.75) w l 4");
         } else if(spec_units)
	 {
            fprintf(fgnu,", \"H.boezio99\" u 3:4:5 w e %d, \"H.boezio99\" u 3:4 w p %d",1,2);
            fprintf(fgnu,", \"H.menn00\" u 3:4:5 w e %d,   \"H.menn00\" u 3:4 w p %d",1,4);
            fprintf(fgnu,", \"H.sanuki00\" u 3:4:5 w e %d, \"H.sanuki00\" u 3:4 w p %d",1,5);
//            fprintf(fgnu,", \"H.alcaraz00a\" u ?? w l %d, \"H.alcaraz00a\" u 3:4 w p %d",1,7);
            fprintf(fgnu,", 1.6e4*x**(-2.75) w l 4");
         }
      }
      if(datakey == 2)  // Helium
      {
	 if(spec_type)
	 {
            fprintf(fgnu,", \"He_boezio99%s.dat\" w l %d, \"He_boezio99%s.dat\" u 3:4 w p %d",gev,1,gev,2);
            fprintf(fgnu,", \"He_menn00%s.dat\" w l %d,   \"He_menn00%s.dat\" u 3:4 w p %d",gev,1,gev,4);
            fprintf(fgnu,", \"He_sanuki00%s.dat\" w l %d, \"He_sanuki00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"He_alcaraz00%s.dat\" w l %d, \"He_alcaraz00%s.dat\" u 3:4 w p %d",gev,1,gev,6);
            fprintf(fgnu,", \"He_jacee98%s.dat\" w l %d, \"He_jacee98%s.dat\" u 3:4 w p %d",gev,1,gev,7);
            fprintf(fgnu,", \"He_sokol%s.dat\" u 3:4 w p %d",gev,8);
            fprintf(fgnu,", 641*x**(-0.72) w l 4");
         } else if(spec_units)
	 {
            fprintf(fgnu,", \"He.boezio99\" u 3:4:5 w e %d, \"He.boezio99\" u 3:4 w p %d",1,2);
            fprintf(fgnu,", \"He.menn00\" u 3:4:5 w e %d,   \"He.menn00\" u 3:4 w p %d",1,4);
            fprintf(fgnu,", \"He.sanuki00\" u 3:4:5 w e %d, \"He.sanuki00\" u 3:4 w p %d",1,5);
            fprintf(fgnu,", \"He.alcaraz00\" u 3:4:5 w e %d, \"He.alcaraz00\" u 3:4 w p %d",1,6);
            fprintf(fgnu,", \"He0_jacee98.dat\" u 1:4:5 w e %d, \"He0_jacee98.dat\" u 3:4 w p %d",1,7);
            fprintf(fgnu,", \"He.sokol\" u 1:2 w p %d",8);
            fprintf(fgnu,", 641*x**(-2.72) w l 4");
         }
      }
      if(datakey == 5)  // Boron
      {
         fprintf(fgnu,", \"B.dat\" u  1: 2: 3 w e %d, \"B.dat\" u  1: 2 w p %d", 1, 1);
         fprintf(fgnu,", \"B.dat\" u  4: 5: 6 w e %d, \"B.dat\" u  4: 5 w p %d", 1, 2);
      }
      if(datakey == 6)  // Carbon
      {
//         fprintf(fgnu,", \"C_1.dat\" u  1: 2: 3 w e %d, \"C_1.dat\" u  1: 2 w p %d", 1, 1);
         fprintf(fgnu,", \"C_1.dat\" u  4: 5: 6 w e %d, \"C_1.dat\" u  4: 5 w p %d", 1, 2);
         fprintf(fgnu,", \"C_1.dat\" u  7: 8: 9 w e %d, \"C_1.dat\" u  7: 8 w p %d", 1, 3);
         fprintf(fgnu,", \"C_1.dat\" u 10:11:12 w e %d, \"C_1.dat\" u 10:11 w p %d", 1, 4);
         fprintf(fgnu,", \"C_1.dat\" u 13:14:15 w e %d, \"C_1.dat\" u 13:14 w p %d", 1, 5);
         fprintf(fgnu,", \"C_1.dat\" u 16:17:18 w e %d, \"C_1.dat\" u 16:17 w p %d", 1, 6);
         fprintf(fgnu,", \"C_1.dat\" u 19:20:21 w e %d, \"C_1.dat\" u 19:20 w p %d", 1, 7);
         fprintf(fgnu,", \"C_1.dat\" u 22:23:24 w e %d, \"C_1.dat\" u 22:23 w p %d", 1, 8);
         fprintf(fgnu,", \"C_1.dat\" u 25:26:27 w e %d, \"C_1.dat\" u 25:26 w p %d", 1, 9);
         fprintf(fgnu,", \"C_1.dat\" u 28:29:30 w e %d, \"C_1.dat\" u 28:29 w p %d", 1,10);
         fprintf(fgnu,", \"C_2.dat\" u  1: 2: 3 w e %d, \"C_2.dat\" u  1: 2 w p %d", 1, 1);
      }
      if(datakey == 8)  // Oxygen
      {
         fprintf(fgnu,", \"O.dat\" u  1: 2: 3 w e %d, \"O.dat\" u  1: 2 w p %d", 1, 1);
         fprintf(fgnu,", \"O.dat\" u  4: 5: 6 w e %d, \"O.dat\" u  4: 5 w p %d", 1, 2);
         fprintf(fgnu,", \"O.dat\" u  7: 8: 9 w e %d, \"O.dat\" u  7: 8 w p %d", 1, 3);
         fprintf(fgnu,", \"O.dat\" u 10:11:12 w e %d, \"O.dat\" u 10:11 w p %d", 1, 4);
         fprintf(fgnu,", \"O.dat\" u 13:14:15 w e %d, \"O.dat\" u 13:14 w p %d", 1, 5);
         fprintf(fgnu,", \"O.dat\" u 16:17:18 w e %d, \"O.dat\" u 16:17 w p %d", 1, 6);
         fprintf(fgnu,", \"O.dat\" u 19:20:21 w e %d, \"O.dat\" u 19:20 w p %d", 1, 7);
         fprintf(fgnu,", \"O.dat\" u 22:23:24 w e %d, \"O.dat\" u 22:23 w p %d", 1, 8);
         fprintf(fgnu,", \"O.dat\" u 25:26:27 w e %d, \"O.dat\" u 25:26 w p %d", 1, 9);
         fprintf(fgnu,", \"O.dat\" u 28:29:30 w e %d, \"O.dat\" u 28:29 w p %d", 1,10);
         fprintf(fgnu,", \"O.dat\" u 31:32:33 w e %d, \"O.dat\" u 31:32 w p %d", 1,11);
      }
      if(datakey ==26)  // Iron
      {
         fprintf(fgnu,", \"Fe.dat\" u  1: 2: 3 w e %d, \"Fe.dat\" u  1: 2 w p %d", 1, 1);
         fprintf(fgnu,", \"Fe.dat\" u  4: 5: 6 w e %d, \"Fe.dat\" u  4: 5 w p %d", 1, 2);
         fprintf(fgnu,", \"Fe.dat\" u  7: 8: 9 w e %d, \"Fe.dat\" u  7: 8 w p %d", 1, 3);
         fprintf(fgnu,", \"Fe.dat\" u 10:11:12 w e %d, \"Fe.dat\" u 10:11 w p %d", 1, 4);
         fprintf(fgnu,", \"Fe.dat\" u 13:14:15 w e %d, \"Fe.dat\" u 13:14 w p %d", 1, 5);
         fprintf(fgnu,", \"Fe.dat\" u 16:17:18 w e %d, \"Fe.dat\" u 16:17 w p %d", 1, 6);
         fprintf(fgnu,", \"Fe.dat\" u 19:20:21 w e %d, \"Fe.dat\" u 19:20 w p %d", 1, 7);
      }
      fprintf(fgnu," \\\n");
      fclose(fin);
   }

// PLOT ISOTOPIC RATIO

   if(key==2)
   {
      llimit=10; ulimit = 1.e6;
      if(spec_units)  { llimit /=1000.; ulimit /=1000.;}
      fprintf(fgnu,"set logscale x 10\n");
      fprintf(fgnu,"set size 1,0.7\n");
      fprintf(fgnu,"set title \"ISOTOPIC RATIO: (");
      for(i=0; z[i]!=0 || a[i]!=0; i++)
      {
	 if(i>0) fprintf(fgnu,"+");
         fprintf(fgnu,"%2d.%2d", z[i],a[i]);
      }
      fprintf(fgnu,") / (");
      for(i=0; z1[i]>0 || a1[i]>0; i++)
      {
	 if(i>0) fprintf(fgnu,"+");
         fprintf(fgnu,"%2d.%2d", z1[i],a1[i]);
      }
      fprintf(fgnu,")\" 0.,0.\n");
      if(datakey) fprintf(fgnu,"set nokey\n");
      if(!spec_units) fprintf(fgnu,"set xlabel \"Kinetic energy, MeV/nucleon\" 0.,0.\n");
      else fprintf(fgnu,"set xlabel \"Kinetic energy, GeV/nucleon\" 0.,0.\n");
      fprintf(fgnu,"set ylabel \"\" 0.,0.\n");

      fin = fopen(argv[1],"r");
      if(datakey ==-1)  //  pbar/p
      {
         fprintf(fgnu,"set label \"Phi = %.0f MV\" at %.1f,%.1g\n",Phi,llimit*300, pow(10,(int) log10(ymax*1.e-3)));
         fprintf(fgnu,"set logscale y 10\n");
         fprintf(fgnu,"plot [%.1g:%.1g][%.1g:%.1g] ",llimit,ulimit, pow(10,(int) log10(ymax*1.e-4)),pow(10,(int)log10(ymax*1.1))); 
      } else
      {
         fprintf(fgnu,"set label \"Phi = %.0f MV\" at %.1g,%.3g\n",Phi,llimit*300,ymax/4);
         fprintf(fgnu,"plot [%.1g:%.1g][0:%.2g] \\\n",llimit,ulimit,ymax*1.33); 
      }
      for(i=0; i<argc-1; i++)
      {
	 if(i!=0) fprintf(fgnu,"\\\n, ");
         fprintf(fgnu,"\"%s\" w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:3 w l %d",argv[i+1],i+1);
      }
      if(datakey ==-1) // pbar/p
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,",\"pbar2p_hof96%s.dat\" w l %d,\"pbar2p_hof96%s.dat\" u 1:3 w p %d\\\n",gev,1,gev,2);
         fprintf(fgnu,",\"pbar2p_bergstroem00%s.dat\" w l %d,\"pbar2p_bergstroem00%s.dat\" u 1:3 w p %d\\\n",gev,1,gev,4);
         fprintf(fgnu,",\"pbar2p_stochaj00%s.dat\" w l %d,\"pbar2p_stochaj00%s.dat\" u 1:3 w p %d\\\n",gev,1,gev,5);
         fprintf(fgnu,",\"pbar2p_orito00%s.dat\" w l %d,\"pbar2p_orito00%s.dat\" u 1:3 w p %d",gev,1,gev,6);
         fprintf(fgnu,",\"pbar2p_maeno01%s.dat\" w l %d,\"pbar2p_maeno01%s.dat\" u 1:3 w p %d",gev,1,gev,8);
      }
      if(datakey == 2) // He3/He4
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"He_rat%s.dat\" u  2: 3: 4 w e %d, \"He_rat%s.dat\" u  2: 3 w p %d, \"He_rat%s.dat\" u  1: 3 w l %d\\\n",gev,1,gev,4,gev,1);
      }
      if(datakey == 4) // Be10/Be9
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"Be_rat%s.dat\" u  2: 3: 4 w e %d, \"Be_rat%s.dat\" u  2: 3 w p %d, \"Be_rat%s.dat\" u  1: 3 w l %d\\\n",gev,1,gev,4,gev,1);
         fprintf(fgnu,", \"Be_rat%s.dat\" u  6: 7: 8 w e %d, \"Be_rat%s.dat\" u  6: 7 w p %d, \"Be_rat%s.dat\" u  5: 7 w l %d\\\n",gev,1,gev,7,gev,1);
         fprintf(fgnu,", \"Be_rat%s.dat\" u 10:11:12 w e %d, \"Be_rat%s.dat\" u 10:11 w p %d, \"Be_rat%s.dat\" u  9:11 w l %d\\\n",gev,1,gev,6,gev,1);
         fprintf(fgnu,", \"Be_rat%s.dat\" u 14:15:16 w e %d, \"Be_rat%s.dat\" u 14:15 w p %d, \"Be_rat%s.dat\" u 13:15 w l %d",gev,1,gev,5,gev,1);
      }
      if(datakey == 5) // B/C ratio
      {
	 fprintf(fgnu,"\\\n");
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u  1: 2: 3 w e %d, \"BC_rat1%s.dat\" u  1: 2 w p %d\\\n",gev,1,gev, 1);
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u  4: 5: 6 w e %d, \"BC_rat1%s.dat\" u  4: 5 w p %d\\\n",gev,1,gev, 2);
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u  7: 8: 9 w e %d, \"BC_rat1%s.dat\" u  7: 8 w p %d\\\n",gev,1,gev, 3);
         fprintf(fgnu,", \"BC_rat1%s.dat\" u 10:11:12 w e %d, \"BC_rat1%s.dat\" u 10:11 w p %d\\\n",gev,1,gev, 4);
         fprintf(fgnu,", \"BC_rat1%s.dat\" u 13:14:15 w e %d, \"BC_rat1%s.dat\" u 13:14 w p %d\\\n",gev,1,gev, 5);
         fprintf(fgnu,", \"BC_rat1%s.dat\" u 16:17:18 w e %d, \"BC_rat1%s.dat\" u 16:17 w p %d\\\n",gev,1,gev, 6);
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u 19:20:21 w e %d, \"BC_rat1%s.dat\" u 19:20 w p %d\\\n",gev,1,gev, 7);
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u 22:23:24 w e %d, \"BC_rat1%s.dat\" u 22:23 w p %d\\\n",gev,1,gev, 8);
//         fprintf(fgnu,", \"BC_rat1%s.dat\" u 25:26:27 w e %d, \"BC_rat1%s.dat\" u 25:26 w p %d\\\n",gev,1,gev,10);
         fprintf(fgnu,", \"BC_rat1%s.dat\" u 28:29:30 w e %d, \"BC_rat1%s.dat\" u 28:29 w p %d\\\n",gev,1,gev,9);
         fprintf(fgnu,", \"BC_rat2%s.dat\" u  2: 3: 4 w e %d, \"BC_rat2%s.dat\" u  2: 3 w p %d, \"BC_rat2%s.dat\" u 1:3 w l %d\\\n",gev,1,gev,2,gev,1);
         fprintf(fgnu,", \"BC_rat2%s.dat\" u  5: 6: 7 w e %d, \"BC_rat2%s.dat\" u  5: 6 w p %d\\\n",gev,1,gev,7);
         fprintf(fgnu,", \"BC_rat2%s.dat\" u  8: 9:10 w e %d, \"BC_rat2%s.dat\" u  8: 9 w p %d\\\n",gev,1,gev,8);
      }
      if(datakey == 13) // Al26/Al27
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"Al_rat%s.dat\" u  2: 3: 4 w e %d, \"Al_rat%s.dat\" u  2: 3 w p %d, \"Al_rat%s.dat\" u  1: 3 w l %d\\\n",gev,1,gev,3,gev,1);
         fprintf(fgnu,", \"Al_rat%s.dat\" u  6: 7: 8 w e %d, \"Al_rat%s.dat\" u  6: 7 w p %d, \"Al_rat%s.dat\" u  5: 7 w l %d",gev,1,gev,4,gev,1);
      }
      if(datakey == 17) // Cl36/Cl
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"Cl_rat%s.dat\" u  2: 3: 4 w e %d, \"Cl_rat%s.dat\" u  2: 3 w p %d, \"Cl_rat%s.dat\" u  1: 3 w l %d\\\n",gev,1,gev,3,gev,1);
         fprintf(fgnu,", \"Cl_rat%s.dat\" u  6: 7: 8 w e %d, \"Cl_rat%s.dat\" u  6: 7 w p %d, \"Cl_rat%s.dat\" u  5: 7 w l %d",gev,1,gev,4,gev,1);
      }
      if(datakey == 25) // Mn54/Mn
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"Mn_rat%s.dat\" u  2: 3: 4 w e %d, \"Mn_rat%s.dat\" u  2: 3 w p %d, \"Mn_rat%s.dat\" u  1: 3 w l %d\\\n",gev,1,gev,3,gev,1);
         fprintf(fgnu,", \"Mn_rat%s.dat\" u  6: 7: 8 w e %d, \"Mn_rat%s.dat\" u  6: 7 w p %d, \"Mn_rat%s.dat\" u  5: 7 w l %d\\\n",gev,1,gev,4,gev,1);
         fprintf(fgnu,", \"Mn_rat%s.dat\" u  9:10:11 w e %d, \"Mn_rat%s.dat\" u  9:10 w p %d",gev,1,gev,5);
      }
      if(datakey == 26) // Al26/Si28
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"misc_rat%s.dat\" u 14:15:16 w e %d, \"misc_rat%s.dat\" u 14:15 w p %d, \"misc_rat%s.dat\" u 13:15 w l %d",gev,1,gev,3,gev,1);
      }
      if(datakey == 36) // Cl36/Fe56
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"misc_rat%s.dat\" u  2: 3: 4 w e %d, \"misc_rat%s.dat\" u  2: 3 w p %d, \"misc_rat%s.dat\" u  1: 3 w l %d",gev,1,gev,3,gev,1);
      }
      if(datakey == 54) // Mn54/Fe56
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"misc_rat%s.dat\" u  6: 7: 8 w e %d, \"misc_rat%s.dat\" u  6: 7 w p %d, \"misc_rat%s.dat\" u  5: 7 w l %d\\\n",gev,1,gev,3,gev,1);
         fprintf(fgnu,", \"misc_rat%s.dat\" u 10:11:12 w e %d, \"misc_rat%s.dat\" u 10:11 w p %d, \"misc_rat%s.dat\" u  9:11 w l %d",gev,1,gev,4,gev,1);
      }
      if(datakey == 56) // subFe/Fe 
      {
	 fprintf(fgnu,"\\\n");
         fprintf(fgnu,", \"subFe_Fe%s.dat\" u  1: 2: 3 w e %d, \"subFe_Fe%s.dat\" u  1: 2 w p %d\\\n",gev,1,gev,4);
         fprintf(fgnu,", \"subFe_Fe%s.dat\" u  4: 5: 6 w e %d, \"subFe_Fe%s.dat\" u  4: 5 w p %d\\\n",gev,1,gev,5);
         fprintf(fgnu,", \"subFe_Fe%s.dat\" u  7: 8: 9 w e %d, \"subFe_Fe%s.dat\" u  7: 8 w p %d\\\n",gev,1,gev,6);
         fprintf(fgnu,", \"subFe_Fe%s.dat\" u 10:11:12 w e %d, \"subFe_Fe%s.dat\" u 10:11 w p %d\\\n",gev,1,gev,7);
      }
      fprintf(fgnu,"\n");
      fclose(fin);
   }

// PLOT SPECTRA OF ELECTRONS AND POSITRONS

   if(key==3)
   {
      llimit=100; ulimit = 5.e6;
      if(z[0]==1) ulimit = 3.e5;
      if(!spec_type) ulimit =1.e5;
      if(spec_units)  { llimit /=1000.; ulimit /=1000.;}
      yscale = 8.e-12;
      if(spec_type || spec_units) yscale *= 1.e7;
      if(spec_type) yscale *= 100.;
      fprintf(fgnu,"set logscale xy 10\n");
      fprintf(fgnu,"set size 1,0.7\n");
      fprintf(fgnu,"set title \"PARTICLE SPECTRA: ");
      if(z[0]==-1) fprintf(fgnu,"electrons\n");
      if(z[0]== 1) fprintf(fgnu,"positrons\n");
      if(spec_type)                // E^2 Flux
	 if(spec_units)            // GeV/m^2
         { 
            fprintf(fgnu,"set xlabel \"Kinetic energy, GeV\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"E^2 Flux, GeV m^-2 s^-1 sr^-1\" 0.,0.\n");
         } else                    // MeV/cm^2
         { 
            fprintf(fgnu,"set xlabel \"Kinetic energy, MeV\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"E^2 Flux, MeV cm^-2 s^-1 sr^-1\" 0.,0.\n");
         }
      else                         // Flux
	 if(spec_units)            // 1/(m^2 GeV)
         {
            fprintf(fgnu,"set xlabel \"Kinetic energy, GeV\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"Flux, m^-2 s^-1 sr^-1 GeV^-1\" 0.,0.\n");
         } else                    // 1/(cm^2 MeV)
         {
            fprintf(fgnu,"set xlabel \"Kinetic energy, MeV\" 0.,0.\n");
            fprintf(fgnu,"set ylabel \"Flux, cm^-2 s^-1 sr^-1 MeV^-1\" 0.,0.\n");
         }
//      fprintf(fgnu,"set xrange [*: *] noreverse nowriteback  # (currently [-10:10] )\n");
//      fprintf(fgnu,"set yrange [*:*] noreverse nowriteback  # (currently [-10:10] )\n");
      fprintf(fgnu,"set zlabel \"\" 0.000000,0.000000\n");
//      fprintf(fgnu,"set zrange [*:*] noreverse nowriteback  # (currently [-10:10] )\n");
//      fprintf(fgnu,"set locale \"C\"\n");
      fprintf(fgnu,"set label \"Phi = %.0f MV\" at %.1g,%.1g\n",Phi, 10*llimit, 5*yscale);

      if(datakey) fprintf(fgnu,"set nokey\n");
      fin = fopen(argv[1],"r");
      fprintf(fgnu,"plot [%.1g:%.1g][%.1g:*] ",llimit, ulimit, yscale ); 
      for(i=0; i<argc-1; i++)
      {
	 if(i!=0) fprintf(fgnu,", ");
         fprintf(fgnu,"\"%s\" w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:3 w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:4 w l %d,",argv[i+1],i+1);
         fprintf(fgnu,"\"%s\" u 1:5 w l %d ",argv[i+1],i+1);
      }
      if(datakey ==-1)  // electrons
      {
	 if(spec_type)
	 {
            fprintf(fgnu,", \"electrons_boezio00%s.dat\" w l %d, \"electrons_boezio00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"electrons_barwick98%s.dat\" w l %d, \"electrons_barwick98%s.dat\" u 3:4 w p %d",gev,1,gev,6);
//            fprintf(fgnu,", \"electrons_taira93%s.dat\" w l %d",gev,-1);
            fprintf(fgnu,", \"electrons_kobayashi99%s.dat\" w l %d, \"electrons_kobayashi99%s.dat\" u 3:4 w p %d",gev,1,gev,7);
            fprintf(fgnu,", \"electrons_grimani02%s.dat\" w l %d, \"electrons_grimani02%s.dat\" u 3:4 w p %d",gev,1,gev,8);
            fprintf(fgnu,", \"electrons_duvernois01%s.dat\" w l %d, \"electrons_duvernois01%s.dat\" u 3:4 w p %d",gev,1,gev,9);
            fprintf(fgnu,", \"electrons_alcaraz00%s.dat\" w l %d, \"electrons_alcaraz00%s.dat\" u 3:4 w p %d",gev,1,gev,10);
         } else
         {
            fprintf(fgnu,", \"electrons0_boezio00.dat\" w l %d, \"electrons0_boezio00.dat\" u 3:4 w p %d",1,5);
            fprintf(fgnu,", \"electrons0_barwick98.dat\" w l %d, \"electrons0_barwick98.dat\" u 3:4 w p %d",1,6);
            fprintf(fgnu,", \"electrons.kobayashi99\" u 3:4:5 w e %d, \"electrons.kobayashi99\" u 3:4:5 w e %d",1,7);
            fprintf(fgnu,", \"electrons0_grimani02%s.dat\" w l %d, \"electrons0_grimani02%s.dat\" u 3:4 w p %d",gev,1,gev,8);
            fprintf(fgnu,", \"electrons0_duvernois01%s.dat\" w l %d, \"electrons0_duvernois01%s.dat\" u 3:4 w p %d",gev,1,gev,9);
            fprintf(fgnu,", \"electrons0_alcaraz00%s.dat\" w l %d, \"electrons0_alcaraz00%s.dat\" u 3:4 w p %d",gev,1,gev,10);
         }
      }
      if(datakey == 1)  // positrons
      {
	 if(spec_type)
	 {
            fprintf(fgnu,", \"positrons_boezio00%s.dat\" w l %d, \"positrons_boezio00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"positrons_barwick98%s.dat\" w l %d, \"positrons_barwick98%s.dat\" u 3:4 w p %d",gev,1,gev,6);
            fprintf(fgnu,", \"positrons_grimani02%s.dat\" w l %d, \"positrons_grimani02%s.dat\" u 3:4 w p %d",gev,1,gev,7);
            fprintf(fgnu,", \"positrons_duvernois01%s.dat\" w l %d, \"positrons_duvernois01%s.dat\" u 3:4 w p %d",gev,1,gev,8);
            fprintf(fgnu,", \"positrons_alcaraz00%s.dat\" w l %d, \"positrons_alcaraz00%s.dat\" u 3:4 w p %d",gev,1,gev,9);
         } else
         {
            fprintf(fgnu,", \"positrons0_boezio00%s.dat\" w l %d, \"positrons0_boezio00%s.dat\" u 3:4 w p %d",gev,1,gev,5);
            fprintf(fgnu,", \"positrons0_barwick98%s.dat\" w l %d, \"positrons0_barwick98%s.dat\" u 3:4 w p %d",gev,1,gev,6);
            fprintf(fgnu,", \"positrons0_grimani02%s.dat\" w l %d, \"positrons0_grimani02%s.dat\" u 3:4 w p %d",gev,1,gev,7);
            fprintf(fgnu,", \"positrons0_duvernois01%s.dat\" w l %d, \"positrons0_duvernois01%s.dat\" u 3:4 w p %d",gev,1,gev,8);
            fprintf(fgnu,", \"positrons0_alcaraz00%s.dat\" w l %d, \"positrons0_alcaraz00%s.dat\" u 3:4 w p %d",gev,1,gev,9);
         }
      }
      if(datakey == 0)  // e+/etot
      {
      }
      fprintf(fgnu,"\n");
      fclose(fin);
   }

   fclose(fgnu);
   return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void readheader ( char* filename)
{
   fitsfile *fptr;         // pointer to the FITS file, defined in fitsio.h 
   int status=0, nkeys, keypos, hdutype, ii, jj;
   char card[FLEN_CARD];   // standard string lengths defined in fitsioc.h 

   if ( fits_open_file(&fptr, filename, READONLY, &status) ) printerror( status );

// attempt to move to next HDU, until we get an EOF error 
    for (ii = 1; !(fits_movabs_hdu(fptr, ii, &hdutype, &status) ); ii++) 
    {
// get no. of keywords 
       if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) ) printerror( status );
       printf("Header listing for HDU #%d:\n", ii);
       for (jj = 1; jj <= nkeys; jj++)  
       {
          if ( fits_read_record(fptr, jj, card, &status) ) printerror( status );
          printf("%s\n", card);            // print the keyword card 
        }
        printf("END\n\n");                 // terminate listing with END 
   }
   if (status == END_OF_FILE) status = 0;  // got the expected EOF error; reset = 0  
   else printerror( status );              // got an unexpected error                
   if ( fits_close_file(fptr, &status) ) printerror( status );
   return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void printerror( int status)
{
   if (status)
   {
      fits_report_error(stderr, status); // print error report 
      exit( status );    // terminate the program, returning error status 
   }
   return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

float add_nuclei(int option, long ZZ, long AA, float energy, float scale)
{
   static int i=0, Z[99],A[99];
   static float abund[99];
   int j;
   float spectrum=0., Jacobian, amu=931.494, Pa, Po,  rho,rho0; // MeV/nucleon

   if(option != 0)
   {
      for(j=0; j<i; j++) 
	 if(100*ZZ+AA == 100*Z[j]+A[j]) break; // printf("%4d%4d  %.3e\n",Z[j],A[j],abund[j]);
      if(j==i) return 0.;
      Pa = sqrt(pow(energy+amu,2) -pow(amu,2)); // momentum per nucleon
      Po = sqrt(pow(scale +amu,2) -pow(amu,2)); // momentum scale
   }

   switch(option)
   {
   case 0: // read abundancies from a file
      FILE *fin;
      fin = fopen("add-abund.dat","r");
      while(fscanf(fin,"%d%d%e",Z+i,A+i,abund+i)!=EOF) i++;
//      for(j=0;j<i;j++) printf("%4d%4d  %.3e\n",Z[j],A[j],abund[j]);
      fclose(fin);
      break;

   case 1:
      spectrum = (energy/scale >100.) ? 0.: abund[j]*exp(-energy/scale)/energy;
//      cout<<" >> "<<energy<<" "<<spectrum<<endl;
      break;

   case 2:
      Jacobian = (energy+amu)/Pa;
      spectrum = (energy/scale >100.) ? 0.: abund[j]*Jacobian/Pa*exp(-Pa/Po);
      break;

   case 3:
      rho = AA*1./ZZ*sqrt(pow(energy+amu,2) -pow(amu,2)); // rigidity
      rho0= 2.      *sqrt(pow(scale +amu,2) -pow(amu,2)); // scale rigidity for A/Z=2.
      spectrum = (energy/scale >100.) ? 0.: abund[j]/rho*exp(-rho/rho0);
      break;

   case 4:
      rho = AA*1./ZZ*sqrt(pow(energy+amu,2) -pow(amu,2)); // rigidity
      rho0= AA*1./ZZ*sqrt(pow(scale +amu,2) -pow(amu,2)); // scale rigidity for A/Z=2.
      Jacobian = AA*1./ZZ*(energy+amu)/Pa;
      spectrum = (energy/scale >100.) ? 0.: abund[j]*Jacobian/rho*exp(-rho/rho0);
      break;

   default:
      return 0.;
   }
//cout<<" energy,spectrum = "<<energy<<" "<<spectrum<<endl;
   return spectrum;
}

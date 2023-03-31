
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_pi0_decay_emiss.cc *                      galprop package * 2/16/2005
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]

emissivity (cm^-3 s^-1 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,p_beam }  ) n(E)E dlog(E)]

pp_meson has Egamma in GeV, beam momentum in GeV
The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
BUT UNITS OF density/momentum = flux/(KE/nucleon)..... CHECK CAREFULLY, ALSO factor A
cross section from pp_meson in barns GeV^-1
factor= 1.0e-24* *1.0e-3 log(Ekin_factor)
*/
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

int Galprop::gen_pi0_decay_emiss()   // generate pi0-decay emissivity
{
   cout<<" >>>> gen_pi0_decay_emiss"<<endl;
   cout<<"generating pi0-decay emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   double factor,cs_p_HI,cs_p_He,cs_He_HI,cs_He_He, cs[2][4], y, P=-1.,Ekin=1.,Etot,beta,gamma,rig, yp,fp,yh,fh;
   double x,EL,Ep_min,Mp1=1.e3*Mp,Mr=Mpi0/Mp,Pth0=780.; //MeV/amu,Pth0-threshold momentum
   int stat=0, NA1, NA2, iprotons, iHelium, i, key,test;
   Distribution protons;

// identify CR protons
   if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid,                 gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   protons=0.;
   for(i=0, iprotons=-1; i<n_species; i++)  
      if(101==100*gcr[i].Z+gcr[i].A)
      {
         iprotons=i;
 	 protons+=gcr[iprotons].cr_density;
         cout<<"  CR protons found as species #"<<iprotons<<endl;
      }
   if(iprotons==-1) { cout<<"  CR protons not found!"<<endl; protons.delete_array(); return 1; }

// identify CR Helium
   for(i=0, iHelium =-1; i<n_species; i++) if(204==100*gcr[i].Z+gcr[i].A) iHelium =i;
   if(iHelium ==-1) cout<<"  CR Helium  not found!"<<endl;
   else cout<<"  CR Helium  found as species #"<<iHelium <<endl;

   int key1=0; // pi0-decay
   galaxy.pi0_decay_emiss=0.;
   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
   {
      for(i=0,key=0;i<4;i++)  cs[0][i] = cs[1][i] = 0.;
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {  // energy conservation check
	 if(gcr[iprotons].p[ip+1]<=Pth0) if(iHelium !=-1) { if(gcr[iHelium].p[ip+1]/gcr[iHelium].A<=Pth0) continue; }
	                                 else continue;

         if(galdef.integration_mode==1) // ### old integration ###
	 {         //pp_meson_cc( Esec, Pp1, NA1, NA2, key1 );
            NA1=1;    NA2=1; //           beam+target     p+HI
            cs_p_HI = (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                      blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
            NA1=1;    NA2=4; //           beam+target     p+He
            cs_p_He= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1): //IMOS20050216
                                             blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
            if(iHelium !=-1)
	    { 
               NA1=4;    NA2=1; //           beam+target     He+HI
               cs_He_HI= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1): //IMOS20050216
                                                 blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
               NA1=4;    NA2=4; //           beam+target     He+He
               cs_He_He= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1): //IMOS20050216
                                                 blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1);
	       if(galdef.verbose==-478) cout<<"  cs_p_HI, cs_p_He, cs_He_HI, cs_He_He ="<<cs_p_HI<<"  "<<cs_p_He<<"  "<<cs_He_HI<<"  "<<cs_He_He<<endl;
            }
            if(                cs_p_HI==0.) continue;

         } else // ### new integration - analytical ###
	 {
	    if(ip>=gcr[iprotons].n_pgrid-1) break; // don't consider the last point since already included
	    for(i=0;i<4;i++)  cs[0][i] = cs[1][i]; // reassign old points

// min energy of pi0-meson contributing to E_gamma
            EL = galaxy.E_gamma[iEgamma] +pow(Mpi0/2.*1.e3,2)/galaxy.E_gamma[iEgamma];
            EL = (EL< Mpi0*1.e3) ? Mpi0*1.e3 :EL;

            x=EL-Mp1/2.*(1.+Mr*Mr);
// min total energy of proton contributing to E_gamma
            if(fabs(x)<10.)  // avoiding singularity in denominator
	       Ep_min= 0.5*( Mp1 +Mp1*Mr*Mr/2. +Mp1*( 9.*(1.+Mr*Mr)+5.*pow(1.-Mr*Mr,2) )/6./(1.-Mr*Mr)
                            +x 
		            +x*( pow(1.-Mr*Mr,2) +20.*(1.+Mr*Mr)+9. )/6./(1.-Mr*Mr)
		            -x*pow( 36.*(1.+Mr*Mr) +20.*pow(1.-Mr*Mr,2),2 )/864./pow(1.-Mr*Mr,3) );
            else
	       Ep_min= Mp1/2./( 1. +Mr*Mr -2.*EL/Mp1) 
                          *(2. -Mr*Mr +EL/Mp1*Mr*Mr -2.*pow(EL/Mp1,2)
                           -sqrt( (pow(EL/Mp1,2) -Mr*Mr) *(4.*pow(EL/Mp1,2) +16.*EL/Mp1 +Mr*Mr*(-4.*EL/Mp1 -8. +Mr*Mr)) ) );
// min kinetic energy of proton contributing to E_gamma
            Ekin=Ep_min-Mp1;
            if(Ekin>gcr[iprotons].Ekin[ip+1]) continue;

//************************************** TEST
if(galdef.verbose==-217)
if(iEgamma==galaxy.n_E_gammagrid/2)
   cout<<" Egam= "<<galaxy.E_gamma[iEgamma]<<" Ep_kin= "<<gcr[iprotons].Ekin[ip]<<" cs0,1= "<<cs[0][0]<<"  "<<cs[1][0]<<" EL,Ekin= "<<EL<<"  "<<Ekin<<endl;
//if(cs[0][0]!=0.) break;
//***************************************/

// calculate new points in Etot
            NA1=1;    NA2=1; //           beam+target     p+HI
            cs[1][0]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                      blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1);
            NA1=1;    NA2=4; //           beam+target     p+He
            cs[1][1]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                      blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1);
            if(iHelium !=-1)
	    { 
               NA1=4;    NA2=1; //           beam+target     He+HI
               cs[1][2]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip+1]*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                         blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1);
               NA1=4;    NA2=4; //           beam+target     He+He
               cs[1][3]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, gcr[iHelium ].p[ip+1]*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                         blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, gcr[iprotons].p[ip+1]*1.e-3, NA1, NA2, key1);
            }
	    if(galdef.verbose==-478) cout<<"  cs[0][0], cs[0][1], cs[0][2], cs[0][3] ="<<cs[0][0]<<"  "<<cs[0][1]<<"  "<<cs[0][2]<<"  "<<cs[0][3]<<endl;
// fool proof
            if(                cs[1][0]==0.) continue;
// lower integration limit falls between the grid points 
            if(cs[0][0]==0. && cs[1][0]!=0.) // test p+p reaction
	    {
	       key=1;
	       Ekin =0.02*(gcr[iprotons].Ekin[ip+1] +49.*Ekin); //calc.lower int.limit
	       P =-1.;
               kinematic(gcr[iprotons].Z,gcr[iprotons].A,"nucleus",P,Ekin,Etot,beta,gamma,rig,test);
               NA1=1;    NA2=1; //           beam+target     p+HI
               cs[0][0]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                         blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, NA1, NA2, key1);
               NA1=1;    NA2=4; //           beam+target     p+He
               cs[0][1]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, NA1, NA2, key1): //IMOS20050216
	                                         blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3, NA1, NA2, key1);
               if(iHelium !=-1)
	       { 
                  NA1=4;    NA2=1; //           beam+target     He+HI
                  if(iHelium !=-1) 
		    cs[0][2]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3*gcr[iHelium].A, NA1, NA2, key1): //IMOS20050216
		                                      blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3*gcr[iHelium].A, NA1, NA2, key1);
                  NA1=4;    NA2=4; //           beam+target     He+He
                  if(iHelium !=-1) 
		    cs[0][3]= (galdef.pi0_decay==1) ? pp_meson_cc   (galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3*gcr[iHelium].A, NA1, NA2, key1): //IMOS20050216
		                                      blattnig_gamma(galaxy.E_gamma[iEgamma]*1.e-3, P*1.e-3*gcr[iHelium].A, NA1, NA2, key1);
	       }
	       if(galdef.verbose==-478) cout<<"  cs[0][0], cs[0][1], cs[0][2], cs[0][3] ="<<cs[0][0]<<"  "<<cs[0][1]<<"  "<<cs[0][2]<<"  "<<cs[0][3]<<endl;
//************************************** TEST
if(galdef.verbose==-217)
cout<<" lower integration limit falls between the grid points "<<cs[0][0]<<" "<<cs[1][0]<<" e= "<<gcr[iprotons].Ekin[ip]<<" g= "<<galaxy.E_gamma[iEgamma]<<endl;
//***************************************/
            }
         }

         if(galaxy.n_spatial_dimensions==2)
            for(int ir=0; ir<gcr[iprotons].n_rgrid-1; ir++)
               for(int iz=1; iz<gcr[iprotons].n_zgrid-1; iz++)
               { 
//************************************** TEST
if(galdef.verbose==-217)
{
if(galdef.integration_mode==1 && ip>0) Ekin=0.5*(gcr[iprotons].Ekin[ip-1]+gcr[iprotons].Ekin[ip  ]);
cs_p_HI=cs_p_He=cs_He_HI=cs_He_He=1.;
for(i=0;i<4;i++)  cs[0][i] = cs[1][i] = 1.;
protons.d2[ir][iz].s[ip  ] = pow(gcr[iprotons].Ekin[ip  ],-3);
protons.d2[ir][iz].s[ip+1] = pow(gcr[iprotons].Ekin[ip+1],-3);
if(iHelium !=-1) gcr[iHelium ].cr_density.d2[ir][iz].s[ip  ] =pow(gcr[iHelium ].Ekin[ip  ],-3);
if(iHelium !=-1) gcr[iHelium ].cr_density.d2[ir][iz].s[ip+1] =pow(gcr[iHelium ].Ekin[ip+1],-3);
galdef.He_H_ratio=0.;
}
//***************************************/

                  if(galdef.integration_mode==1) // ### old integration ###
	          {
                     galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]+=       
                        (cs_p_HI +cs_p_He *galdef.He_H_ratio) *protons.                 d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip];

                     if(iHelium !=-1) galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]+=       
                        (cs_He_HI+cs_He_He*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium ].A;

		     if(galdef.verbose==-701)//selectable debug AWS20040303
		     {
		      if(iEgamma==10)
	              {
		      cout<<"ir iz E_gamma pi0_decay_emiss "
                      <<ir<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma];
                      cout<<" gcr[iprotons].Ekin[ip]="<<gcr[iprotons].Ekin[ip];
		      cout<<" cs_p_HI cs_p_He cs_He_HI+cs_He_He "<<cs_p_HI<<" "<< cs_p_He<<" "<<  cs_He_HI<<" "<<  cs_He_He;
		      cout<<" protons= "<<protons.d2[ir][iz].s[ip]<<" Helium= "<<gcr[iHelium ].cr_density.d2[ir][iz].s[ip]<<endl;
		      }
		     }
                     continue;
                  }


                 if(galdef.integration_mode==0)// AWS20040303
                 {
// ### new integration - analytical ###
		  if(key==1)  // case: lower integration limit falls between the grid points
		  { // interpolate proton spectrum (power-law)
                     yp=log(     protons.d2[ir][iz].s[ip  ]/         protons.d2[ir][iz].s[ip+1])     // yp= power-law index
                       /log(gcr[iprotons].       Ekin[ip  ]/    gcr[iprotons].       Ekin[ip+1]);
                     fp=         protons.d2[ir][iz].s[ip  ]*pow(gcr[iprotons].       Ekin[ip  ],-yp);// normalization
                     fp*= pow(Ekin,yp);                                                              // proton flux @ Ekin

                    // interpolate Helium spectrum (power-law)
                     if(iHelium !=-1)
		     {
                        yh=log(gcr[iHelium ].cr_density.d2[ir][iz].s[ip  ]/    gcr[iHelium ].cr_density.d2[ir][iz].s[ip+1])
                          /log(gcr[iHelium ].                   Ekin[ip  ]/    gcr[iHelium].                    Ekin[ip+1]);
                        fh=    gcr[iHelium ].cr_density.d2[ir][iz].s[ip  ]*pow(gcr[iHelium].                    Ekin[ip  ],-yh);
                        fh*= pow(Ekin,yh);                                           // He flux @ Ekin
                     }
		  } else
		  {  // case: integration between the grid points
		     Ekin=gcr[iprotons].                   Ekin[ip  ]; // Ekin -kin.energy/nucleon, same for p & He
                     fp =      protons.d2[ir][iz].            s[ip  ];
                     if(iHelium !=-1) fh = gcr[iHelium ].cr_density.d2[ir][iz].s[ip  ];
                  }

// protons: derive y (= power-law index of the total expression)
		  y =log((cs[0][0] +cs[0][1]*galdef.He_H_ratio) *fp
	               /((cs[1][0] +cs[1][1]*galdef.He_H_ratio) *protons.d2[ir][iz].s[ip+1]))
                    /log(Ekin                              /gcr[iprotons].       Ekin[ip+1]);
// integrate analytically
                  galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma] +=(cs[1][0] +cs[1][1]*galdef.He_H_ratio) 
                              *protons.  d2[ir][iz].s[ip+1]
                         *gcr[iprotons].         Ekin[ip+1]/(y+1.)
                 *(1.-pow(                       Ekin/
                          gcr[iprotons].         Ekin[ip+1],y+1.));

// Helium: derive y (= power-law index of the total expression)
                  if(iHelium !=-1)
		  { 
		     y =log((cs[0][2] +cs[0][3]*galdef.He_H_ratio) *fh
	                  /((cs[1][2] +cs[1][3]*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d2[ir][iz].s[ip+1]))
                       /log(Ekin                                   /gcr[iHelium ].                   Ekin[ip+1]);
// integrate analytically
                     galaxy.pi0_decay_emiss.        d2[ir][iz].s[iEgamma] +=(cs[1][2] +cs[1][3]*galdef.He_H_ratio) 
                          *gcr[iHelium ].cr_density.d2[ir][iz].s[ip+1] *gcr[iHelium ].A
                          *gcr[iHelium ].                   Ekin[ip+1]/(y+1.)
                  *(1.-pow(                                 Ekin/
                           gcr[iHelium ].                   Ekin[ip+1],y+1.));
                  }
//************************************** TEST
//if((ir==0 && iz==1 && iEgamma==23) || (ir==0 && iz==1 && iEgamma==13))
//cout<<" protons @ ir,iz,iEgamma= >>> "<<ip<<" "<<protons.d2[ir][iz].s[ip]<<" "<<protons.d2[ir][iz].s[ip+1]<<" "<<" Egam= "<<galaxy.E_gamma[iEgamma]<<" Ep_kin= "<<gcr[iprotons].Ekin[ip]<<" Ep_kin+1= "<<gcr[iprotons].Ekin[ip+1]<<" cs= "<<cs[1][0]<<" "<<cs[1][1]<<" y= "<<y<<endl;
//cout<<"ir iz iEgamma pi0_decay_emiss "<<ir<<" "<<iz<<" "<<iEgamma<<" "<<galaxy.pi0_decay_emiss.d2[ir][iz].s[iEgamma]<<endl;
//***************************************/

		    } // integration_mode==0

               }//iz //ir //particle.n_spatial_dimensions==2



//************************************** TEST
if(galdef.verbose==-217)
if(iEgamma==galaxy.n_E_gammagrid/2)
  cout<<" Ekin= "<<Ekin<<"  "<<gcr[iprotons].Ekin[ip]<<" integral= "<<galaxy.pi0_decay_emiss.d2[gcr[iprotons].n_rgrid/2][gcr[iprotons].n_zgrid/2].s[iEgamma]<<"  >>>"<<pow(Ekin,-2)/2.<<"  " <<pow(Ekin,-2)*5./2.<<endl;
//***************************************/






         if(galaxy.n_spatial_dimensions==3)
            for(int ix=1; ix<gcr[iprotons].n_xgrid-1; ix++)
               for(int iy=1; iy<gcr[iprotons].n_ygrid-1; iy++)
                  for(int iz=1; iz<gcr[iprotons].n_zgrid-1; iz++)
                  {
                     if(galdef.integration_mode==1) // ### old integration ###
	             {
                        galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]+=           
                           (cs_p_HI +cs_p_He *galdef.He_H_ratio) *gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip];
                        if(iHelium !=-1) galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]+=            
                           (cs_He_HI+cs_He_He*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium ].A;
             
//cout<<"ix iy  iz  E_gamma_pi0_decay_emiss "<<ix<<" "<<iy<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma]<<endl; 
                        continue;
                     }
// ### new integration - analytical ###
		     if(key==1)  // case: lower integration limit falls between the grid points
		     { // interpolate proton spectrum (power-law)
                        yp=log(     protons.d3[ix][iy][iz].s[ip  ]/         protons.d3[ix][iy][iz].s[ip+1])     // power-law index
                          /log(gcr[iprotons].           Ekin[ip  ]/    gcr[iprotons].           Ekin[ip+1]);
                        fp=         protons.d3[ix][iy][iz].s[ip  ]*pow(gcr[iprotons].           Ekin[ip  ],-yp);// normalization
                        fp*= pow(Ekin,yp);                                                                      // prot.flux @ Ekin

                    // interpolate Helium spectrum (power-law)
                        if(iHelium !=-1)
		        {
                           yh=log(gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip  ]/    gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip+1])
                             /log(gcr[iHelium ].                       Ekin[ip  ]/    gcr[iHelium].                        Ekin[ip+1]);
                           fh=    gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip  ]*pow(gcr[iHelium].                        Ekin[ip  ],-yh);
                           fh*= pow(Ekin,yh);                                           // He flux @ Ekin
                        }
		     } else
		     {  // case: integration between the grid points
		        Ekin=gcr[iprotons].                   Ekin[ip  ]; // Ekin -kin.energy/nucleon, same for p & He
                        fp =      protons.d3[ix][iy][iz].        s[ip  ];
                        if(iHelium !=-1) fh = gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip  ];
                     }
// protons: derive y (= power-law index of the total expression)
		     y =log((cs[0][0] +cs[0][1]*galdef.He_H_ratio) *fp
	                  /((cs[1][0] +cs[1][1]*galdef.He_H_ratio) *protons.d3[ix][iy][iz].s[ip+1]))
                       /log(Ekin                              /gcr[iprotons].           Ekin[ip+1]);
// integrate analytically
                     galaxy.pi0_decay_emiss.d3[ix][iy][iz].s[iEgamma] +=(cs[1][0] +cs[1][1]*galdef.He_H_ratio) 
                                 *protons.  d3[ix][iy][iz].s[ip+1]
                            *gcr[iprotons].             Ekin[ip+1]/(y+1.)
                        *(1.-pow(                       Ekin/
                             gcr[iprotons].             Ekin[ip+1],y+1.));

// Helium: derive y (= power-law index of the total expression)
                     if(iHelium !=-1)
		     { 
		        y =log((cs[0][2] +cs[0][3]*galdef.He_H_ratio) *fh
	                     /((cs[1][2] +cs[1][3]*galdef.He_H_ratio) *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip+1]))
                          /log(Ekin                                   /gcr[iHelium ].                       Ekin[ip+1]);
// integrate analytically
                        galaxy.pi0_decay_emiss.        d3[ix][iy][iz].s[iEgamma] +=(cs[1][2] +cs[1][3]*galdef.He_H_ratio) 
                             *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip+1] *gcr[iHelium ].A
                             *gcr[iHelium ].                       Ekin[ip+1]/(y+1.)
                         *(1.-pow(                                 Ekin/
                              gcr[iHelium ].                       Ekin[ip+1],y+1.));
                     }
                  }//iz //iy //ix //particle.n_spatial_dimensions==3

//cout<<"gen_pi0_decay_emiss: iEgamma ip "<<iEgamma<<" "<<ip<<endl;
	 if(key==1) key=2;
      } //ip
//************************************** TEST
if(galdef.verbose==-217)
if(iEgamma==galaxy.n_E_gammagrid/2) 
{
factor= log(galdef.Ekin_factor);
galaxy.pi0_decay_emiss        *= factor;
if(galdef.integration_mode==1) cout<<"* integral= "<<galaxy.pi0_decay_emiss.d2[gcr[iprotons].n_rgrid/2][gcr[iprotons].n_zgrid/2].s[iEgamma]<<endl;
//exit(1);
}
//***************************************/
   } //iEgamma

   factor= (galdef.integration_mode==1) ? 1.0e-24*1.0e-3* log(galdef.Ekin_factor) : 1.0e-24*1.0e-3;
   galaxy.pi0_decay_emiss *= factor;
   protons.delete_array();                  // IMOS20030217

if(galdef.verbose==-702) // selectable debug AWS20041214
 { cout<<"   pi0-decay emissivity "<<endl; galaxy.pi0_decay_emiss.print(); }
   cout<<" <<<< gen_pi0_decay_emiss"<<endl;

if(galdef.verbose==-478) exit(0); // selectable debug IMOS20050216
   return stat;
}

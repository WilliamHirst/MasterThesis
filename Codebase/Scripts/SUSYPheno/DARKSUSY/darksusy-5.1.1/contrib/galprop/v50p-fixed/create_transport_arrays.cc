
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_transport_arrays.cc *                  galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"
//* JE FIX: Following line(s) added for Snow Leopard
#include<cstring>
#include<cstdlib>

static int key=-1;

int Galprop::create_transport_arrays(Particle &particle)
{
   cout<<" >>>> create_transport_arrays "<<particle.name<<endl;
   int stat=0, A1,Z2,A2,K_electron;                                               // IMOS20010816
   int galdef_network_par=0;         // imos network, use everywhere                 IMOS20010816
   double t_half[2];                                                              // IMOS20010816
   double fragment_p,fragment_He,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann; // IMOS20000607
   
// ASSIGNING PRIMARY SOURCE FUNCTION

   cout<<"      assigning primary source function"<<endl;

   double spec_shape;
   double g_0=0., rigid_br0=0.,                                   // IMOS20031012
          rigid_br=galdef.nuc_rigid_br,                           // IMOS20000607
          g_1=galdef.nuc_g_1,
          g_2=galdef.nuc_g_2;

   if(strcmp(particle.name,"primary_electrons")==0)
   {
      g_0=galdef.electron_g_0;                                    // IMOS20031012
      rigid_br0=galdef.electron_rigid_br0;
      g_1=galdef.electron_g_1;
      rigid_br=galdef.electron_rigid_br;
      g_2=galdef.electron_g_2;
   }
   cout<<particle.name<<" g_0="<<g_0<<"  rigid_br0= "<<rigid_br0  // IMOS20031012
                      <<" g_1="<<g_1<<"  rigid_br= " <<rigid_br <<" g_2="<<g_2<<endl;

   particle.primary_source_function = 0.; // IMOS20020418 whole 2D/3D particle.primary_source_function loops are changed

   if(galdef.n_spatial_dimensions==2)
   {
      for(int ip=0; ip<particle.n_pgrid; ip++)
      {
         if(strcmp(galdef.inj_spectrum_type,"Etot")==0) spec_shape=pow(particle.Etot[ip],-g_2); // IMOS20000613
         else                                                                                   // IMOS20000615
	 {
	    if(particle.rigidity[ip]< rigid_br0)                                   // IMOS20031012
               spec_shape =pow(particle.rigidity[ip]/rigid_br0,-g_0) *pow(rigid_br0/rigid_br,-g_1);
            if(rigid_br0<= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_1);
            if(rigid_br <= particle.rigidity[ip])
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_2);
         }
//         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape*=particle.beta[ip];      // IMOS20000615
         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape/=sqrt(1.+pow(2.e3/particle.rigidity[ip],2)); // IMOS20011210

         int ir=0, iz=particle.n_zgrid/2;
         if(galdef.source_specification==1)
	 {
            particle.primary_source_function.d2[ir][iz].s[ip]=particle.primary_abundance*spec_shape;
            continue;
         }
         for(ir=0; ir<particle.n_rgrid; ir++)
         {
            if(galdef.source_specification==2)
	    {
               particle.primary_source_function.d2[ir][iz].s[ip]=particle.primary_abundance*spec_shape;
               continue;
            }
            for(iz=0; iz<particle.n_zgrid; iz++)
            {
               if(galdef.source_specification==0)
               {
                  particle.primary_source_function.d2[ir]               [iz].s[ip]= 
                     source_distribution( particle. r[ir],0.0,particle.z[iz])*particle.primary_abundance*spec_shape;

                  if(galdef.verbose>=20) cout<<"r z source_distribution  "<<particle.r[ir]<<" "<<particle.z[iz]
                     <<" "<<source_distribution(     particle.r[ir],   0.0,particle.z[iz])<<endl;
               }
            }
         }
      }
   }

   if(galdef.n_spatial_dimensions==3)
   {
      for(int ip=0; ip<particle.n_pgrid; ip++)
      {
         if(strcmp(galdef.inj_spectrum_type,"Etot")==0) spec_shape=pow(particle.Etot[ip],-g_2); // IMOS20000613
         else                                                                                   // IMOS20000615
         {
	    if(particle.rigidity[ip]< rigid_br0)                                   // IMOS20031012
               spec_shape =pow(particle.rigidity[ip]/rigid_br0,-g_0) *pow(rigid_br0/rigid_br,-g_1);
            if(rigid_br0<= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_1);
            if(rigid_br <= particle.rigidity[ip])
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_2);
         }
//         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape*=particle.beta[ip];      // IMOS20000615
         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape/=sqrt(1.+pow(2.e3/particle.rigidity[ip],2)); // IMOS20011

	 int ix=particle.n_xgrid/2, iy=particle.n_ygrid/2, iz=particle.n_zgrid/2;
         if(galdef.use_symmetry==1) ix=iy=iz=0;
         if(galdef.source_specification==1)
            particle.primary_source_function.d3[ix][iy][iz].s[ip]=particle.primary_abundance;
         
         for(ix=0; ix<particle.n_xgrid; ix++)
         {
            for(iy=0; iy<particle.n_ygrid; iy++)
            {
               if(galdef.source_specification==2)
                  particle.primary_source_function.d3[ix][iy][iz].s[ip]=particle.primary_abundance;

               for(iz=0; iz<particle.n_zgrid; iz++)
               {
                  if(galdef.source_specification==0)
                     particle.primary_source_function.d3[ix]           [iy]           [iz].s[ip]= 
                        source_distribution(particle.  x[ix],particle.y[iy],particle.z[iz])*particle.primary_abundance;

// source spectral index dispersion
                  double spec_dg_ratio;                             //AWS20010411
                  if(strcmp(particle.name,"primary_electrons")==0)  //AWS20010411
                  { 
	             spec_dg_ratio=
	                pow(particle.rigidity[ip]/galdef.SNR_electron_dgpivot,1.*galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]);

                     if(galdef.verbose==-501)// selectable debug
	             cout<<"SNR_electron_dg="<<galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]
                         <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
		  }//electrons

                  if(strcmp(particle.name,"primary_electrons")!=0)  //AWS20010411
		  {
	             spec_dg_ratio=
	                pow(particle.rigidity[ip]/galdef.SNR_nuc_dgpivot,     1.*galaxy.SNR_nuc_dg.     d3[ix][iy][iz].s[0]);

                     if(galdef.verbose==-501)// selectable debug
                     cout<<"SNR_nuc_dg="<<galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0] 
                         <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
		  }//nucleons

                  particle.primary_source_function.d3[ix][iy][iz].s[ip]*=spec_shape*spec_dg_ratio;
               } //iz
            } //iy
         } //ix
      } //ip
   } // 3D

// CASE: PRIMARY NUCLEI                                                         AWS20000601.1

   if(strcmp(particle.species,"nucleus")==0)                                 // IMOS20000601.1
   {
      particle.primary_source_function*= pow(particle.A, g_2-1);             // IMOS20000613.10
      if(strcmp(galdef.inj_spectrum_type,"Etot")!=0)                         // IMOS20000613.9
         particle.primary_source_function*= pow(fabs(1.*particle.Z),-g_2);   // AWS20000601.2
   }
  
   if(strcmp(particle.name,"primary_electrons")!=0)
   {
      particle.primary_source_function *= galdef.source_normalization;
      cout<<" >>>>>>>>>>>>>>>>>>> norm "<<galdef.source_normalization<<" >>>>>>>>>>>>>>>>>>>"<<endl;
   }

// CASE: PRIMARY ELECTRONS                                                      IMOS20031016

   if(strcmp(particle.name,"primary_electrons")==0)
   {
      particle.primary_source_function *= galdef.electron_source_normalization;
      cout<<" >>>>>>>>>> electron_norm "<<galdef.electron_source_normalization<<" >>>>>>>>>>>>>>>>>>>"<<endl;
   }
  

// ASSIGNING DIFFUSION COEFFICIENT

   cout<<"      assigning diffusion coefficient"<<endl;

// compute beta at break rigidity so that formula will give galdef.D0_xx at this point
   double Ekin_br=-1.;                             // so that p will be used in kinematic
   double p_br=galdef.D_rigid_br*fabs(1.*particle.Z);
   double Etot_br, beta_br, gamma_br, rigidity_br; // output of kinematic
   char species[10];
////////////////////////V//IMOS20030214 all region
   int iprotons=-1;

   if(galdef.diff_reacc > 5)
   {
// identify CR protons
      for(int i=0; i<n_species; i++)  
         if(101==100*gcr[i].Z+gcr[i].A)
         {
            iprotons=i;
            cout<<"  CR protons found as species #"<<iprotons<<endl;
         }
      if(iprotons==-1) { cout<<"CR protons not found!"<<endl; return 1; }
      cout<<"create_transport_arrays>> "<<particle.Z*100+particle.A<<" "<<particle.p[0]<<endl;

// Zero approximation proton spectrum for calculation of damping
      if(galdef.network_iterations==1 && key<0)
      {
         if(galdef.n_spatial_dimensions==2)
            for(int ir=0; ir<particle.n_rgrid; ir++)
               for(int iz=0; iz<particle.n_zgrid; iz++)
                  for(int ip=0; ip<particle.n_pgrid; ip++)
		  {
		    gcr[iprotons].cr_density.d2[ir]    [iz].s[ip] = 
                           (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
                          *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
//cout<<"create_transport_arrays>> "<<gcr[iprotons].cr_density.d2[ir]    [iz].s[ip]<<endl;
                  }

         if(galdef.n_spatial_dimensions==3)
            for(int ix=0; ix<particle.n_xgrid; ix++)
               for(int iy=0; iy<particle.n_ygrid; iy++)
                  for(int iz=0; iz<particle.n_zgrid; iz++)
                     for(int ip=0; ip<particle.n_pgrid; ip++)
                        gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] = 
                           (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
		          *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
      }
      key=1;
   }
//   cout<<"create_transport_arrays>> "<<iprotons<<" "<<gcr[iprotons].cr_density.d2[9][9].s[10]<< " "<<protons.d2[9][9].s[10]<<endl;
//   for(int ip=0; ip<particle.n_pgrid; ip++) cout<<" "<<protons.d2[9][9].s[ip];
//   cout<<endl;
////////////////////////^//IMOS20030214

   strcpy(species,"nucleus");
   if(particle.A==0) strcpy(species,"electron");

// IMOS20000810.2
   if(kinematic(particle.Z,particle.A,species,p_br,Ekin_br,Etot_br,beta_br,gamma_br,rigidity_br,0)) exit(1);

   cout<<" beta at break rigidity="<<beta_br<<"   rigidity_br="<<rigidity_br
      <<"?= galdef.D_rigid_br="<<galdef.D_rigid_br<<endl;

   beta_br=1.0;  // to simulate fortran implementation

   if(galdef.n_spatial_dimensions==2)
     for(int ir=0; ir<particle.n_rgrid; ir++)
       for(int iz=0; iz<particle.n_zgrid; iz++)
	 for(int ip=particle.n_pgrid-1; ip>=0; ip--) // IMOS20060330 changed to reverse order (Wave-particle interactions)
	   {
	     D_xx(particle,iprotons,ir, 0, 0,iz,ip); //array assigned in D_xx IMOS20030129
	     particle.Dpp.d2[ir][iz].s[ip] =D_pp(particle.p[ip],galdef.D_g_1,galdef.v_Alfven,particle.Dxx.d2[ir][iz].s[ip]);
	   }
   if(galdef.n_spatial_dimensions==3)
     for(int ix=0; ix<particle.n_xgrid; ix++)
       for(int iy=0; iy<particle.n_ygrid; iy++)
	 for(int iz=0; iz<particle.n_zgrid; iz++)
	   for(int ip=particle.n_pgrid-1; ip>=0; ip--) // IMOS20060330 changed to reverse order
	     {
	       D_xx(particle,iprotons, 0,ix,iy,iz,ip); //array assigned in D_xx IMOS20030129
	       particle.Dpp.d3[ix][iy][iz].s[ip] =D_pp(particle.p[ip], galdef.D_g_1, galdef.v_Alfven, particle.Dxx.d3[ix][iy][iz].s[ip]);
	     }  // p
   
   if(galdef.verbose>=2)
   {
      cout<< "spatial   diffusion coefficient Dxx  for species "<<particle.name<<endl;
      particle.Dxx.print();
      cout<< "momentum diffusion coefficient Dpp  for species "<<particle.name<<endl;
      particle.Dpp.print();
   }

// ASSIGNING FRAGMENTATION RATE

   cout<<   "======== assigning fragmentation rate ======== "<<endl;
   int ZH=1, ZHe=2;                                                    //  IMOS20010816
   double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;          //  IMOS20010816

   particle.fragment=0.0;
   if(particle.A!=0)
   {
     double CSratio,CStot_ratio;
     
     for(int ip=0; ip<particle.n_pgrid; ip++)
       {                                                             // IMOS20000607 whole segment
	 A1 = 1;                                                    // nucleus
	 Z2 = particle.Z;  A2 = particle.A;                         // - " -
	 if(101==100*Z2+A2) { Z2 = 2;  A2 = 4; }                    // protons
	 if(-99==100*Z2+A2) { A1 =-1;  Z2 = 2;  A2 = 4; }           // antiprotons
	 nucleon_cs(galdef.total_cross_section,particle.Ekin[ip]*1.e-3,A1,Z2,A2, // AWS20010620
		    &PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);
	 He_to_H_CS(particle.Ekin[ip]*1.e-3,particle.Z,particle.A,999,999,&CSratio,&CStot_ratio);
	 
	 fragment_p = PA_inel;                                      // nuclei
	 fragment_He= PA_inel*CStot_ratio;                          // -"- 

// ELECTRON ATTACHMENT/STRIPPING CROSS SECTION                                  IMOS20010816
	 
	 if(galdef.K_capture)
	   {
	     for(K_electron=0;K_electron<=galdef.K_capture;K_electron++) 
	       nucdata(galdef_network_par,particle.Z,particle.A,K_electron,particle.Z,particle.A,&Z2,&A2,&t_half[K_electron]);
	   
	     if(t_half[0] != t_half[1])
	       {
		 Kcapture_cs(particle.Ekin[ip],particle.Z,ZH, &attach_H ,&strip_H );
		 Kcapture_cs(particle.Ekin[ip],particle.Z,ZHe,&attach_He,&strip_He);
		 if(particle.K_electron)
		   {
                     fragment_p += strip_H ;
                     fragment_He+= strip_He;
		   } else
		     {
		       fragment_p += attach_H ;
		       fragment_He+= attach_He;
		     }
		 
		 if(galdef.verbose==-502)// selectable debug
		   cout<<"create_transport_arrays: Ekin Z,A,K_electron,attach_H,strip_H: "
		       <<particle.Ekin[ip]<<" "<<particle.Z<<" "<<particle.A<<" "
		       <<particle.K_electron<<" "<<strip_H<<" "<<attach_H<<endl;
               }
	   }
	 if(101==100*particle.Z+particle.A)                         // protons
	   { 
	     fragment_p = PP_inel;
	     fragment_He= PA_inel;
	   }
	 if(-99==100*particle.Z+particle.A)                         // antiprotons
	   { 
	     fragment_p = aPP_non+aPP_ann;
	     fragment_He= aPA_non+aPA_ann;
	   }

	 if(galdef.n_spatial_dimensions==2)
	   for(int ir=0; ir<particle.n_rgrid; ir++)                   // IMOS20010816
	     for(int iz=0; iz<particle.n_zgrid; iz++)
	       particle.fragment.d2[ir] [iz].s[ip]= particle.beta[ip]*C
		 *(galaxy.n_HI.d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0]+galaxy.n_HII.d2[ir] [iz].s[0]) 
		 *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;
     	 
	 if(galdef.n_spatial_dimensions==3)
	   for(int ix=0; ix<particle.n_xgrid; ix++)                // IMOS20010816
	     for(int iy=0; iy<particle.n_ygrid; iy++)
	       for(int iz=0; iz<particle.n_zgrid; iz++)
		 particle.fragment.d3[ix][iy][iz].s[ip]= particle.beta[ip]*C
		   *(galaxy.n_HI.d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
		   *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;
       }  //  p
   }  //  A!=0
     
   if(galdef.verbose>=2)
     {
       cout<< "fragmentation for species "<<particle.name<<endl;
       particle.fragment.print();
     }

// ASSIGNING MOMENTUM LOSS RATE

   cout<<   "======== assigning momentum loss rate ======== "<<endl;

   if(galdef.n_spatial_dimensions==2)
   {
      for(int ir=0; ir<particle.n_rgrid; ir++)
      {
         for(int iz=0; iz<particle.n_zgrid; iz++)
         {
            for(int ip=0; ip<particle.n_pgrid; ip++)
            {
	       double aion,coul;                                      // NUCLEONS

               if(particle.A!=0) particle.dpdt.d2[ir] [iz].s[ip]= 
                   nucleon_loss(particle.Z,particle.A,particle.Etot[ip],
                   galaxy.n_HI .d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0],
                   galaxy.n_HII.d2[ir] [iz].s[0], galdef.He_H_ratio, 
                   &aion,  &coul) / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1
      
               if(particle.A==0)
               {
                  double uevcm3=0., bevcm3, brem1,brem2,sync,cmptn;   // ELECTRONS
		                                                                       // IMOS200008016
                  bevcm3=pow(galaxy.B_field.d2[ir][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3
                  particle.dpdt.  d2[ir] [iz].s[ip]= electron_loss( particle.Etot[ip],
                     galaxy.n_HI .d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0],
                     galaxy.n_HII.d2[ir] [iz].s[0], galdef.He_H_ratio, uevcm3, bevcm3,
                     &aion, &coul,&brem1,&brem2,&sync,&cmptn)
	             / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1
               }  //  A==0
// cout<<" dpdt="<<particle.dpdt.d2[ix]    [iz].s[ip]<<" aion="<<aion<<endl;
            }  //  p
         }  //  z
      }  //  r
   }

   if(galdef.n_spatial_dimensions==3)
   {
      for(int ix=0; ix<particle.n_xgrid; ix++)
      {
         for(int iy=0; iy<particle.n_ygrid; iy++)
         {
            for(int iz=0; iz<particle.n_zgrid; iz++)
            {
               for(int ip=0; ip<particle.n_pgrid; ip++)
               {
                  double aion,coul;                                        // NUCLEONS

                  if(particle.A!=0) particle.dpdt.d3[ix][iy][iz].s[ip]=
                     nucleon_loss(particle.Z,particle.A,particle.Etot[ip],
                     galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
                     galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio, 
                     &aion,&coul) / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1

                  if(particle.A==0)
                  {
		     double uevcm3=0., bevcm3=0.,brem1,brem2,sync,cmptn;   // ELECTRONS
		                                                                              // IMOS200008016
                     bevcm3=pow(galaxy.B_field.d3[ix][iy][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3
                     particle.  dpdt.d3[ix][iy][iz].s[ip]= electron_loss(particle.Etot[ip],
                        galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
                        galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio, uevcm3, bevcm3,
                        &aion,&coul,&brem1,&brem2,&sync,&cmptn)
                 	/ particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1
        	  }  //  A==0
//  cout<<" dpdt="<<particle.dpdt.d3[ix][iy][iz].s[ip]<<" p="<<particle.p[ip] <<endl;
               }  //  p
            }  //  z
         }  //  y
      }  //  x
   }

// IF ELECTRON ADD KLEIN_NISHINA LOSSES

   if(particle.A==0) e_KN_loss(particle);  // MeV s-1

   if(galdef.verbose>=2)
   {
      cout<< "dpdt for species "<<particle.name<<endl;
      particle.dpdt.print();
   }

// ASSIGNING DECAY RATE

   if(particle.t_half!=0.0)
   {
      cout<<   "======== assigning decay rate ======== "<<endl;

      if(galdef.n_spatial_dimensions==2)
      {
         for(int ir=0; ir<particle.n_rgrid; ir++)
         {
            for(int iz=0; iz<particle.n_zgrid; iz++)
            {
               for(int ip=0; ip<particle.n_pgrid; ip++)
                  particle.decay.d2[ir][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
            }  // z
         }  //  r
      }

      if(galdef.n_spatial_dimensions==3)
      {
         for(int ix=0; ix<particle.n_xgrid; ix++)
         {
            for(int iy=0; iy<particle.n_ygrid; iy++)
            {
               for(int iz=0; iz<particle.n_zgrid; iz++)
               {
                  for(int ip=0; ip<particle.n_pgrid; ip++)
                     particle.decay.d3[ix][iy][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
               }  //  z
            }  //  y
         }  //  x
      }

      if(galdef.verbose>=1)
      {
         cout<< "decay for species "<<particle.name<<endl;
         particle.decay.print();
      }
   }  //  t_half!=0.0

   if(galdef.verbose>=1)
   {
      cout<< "primary source function for species "<<particle.name<<endl;
      particle.primary_source_function.print();
   }

//particle.print();
   cout<<"============== completed creation of transport arrays for "<<particle.name<<endl<<endl;
   cout<<" <<<< create_transport_arrays"<<endl;
   return stat;
}





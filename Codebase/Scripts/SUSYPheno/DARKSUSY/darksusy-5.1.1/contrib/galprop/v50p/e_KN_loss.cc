
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * e_KN_loss.cc *                                galprop package * 8/10/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// energy loss of electron (MeV s^-1) using Klein_Nishina cross-section
//
// ISRF is kept in [nu U_nu] = Hz eV cm-3 Hz-1, transformation to the number
// density is obtained by 
//    nu*d(density)/d(nu) (Hz cm-3 Hz-1) = ISRF * eV_to_erg / (h_planck * nu).
//    energy loss = integral d(nu) d(density)/d(nu) dp/dt 
//                = integral d(log nu) nu*d(density)/d(nu) dp/dt
// d(log nu) is constant in this ISRF, dp/dt is given by e_loss_compton_cc (MeV s^-1).
//
// ISRF_w = (h_planck*nu)* erg_to_eV * 1.0d-6 /m_electron = target photon energy in mc2 units
// factor=(eV_to_erg / h_planck) * LOG(nu(2)/nu(1)) 

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

int Galprop::e_KN_loss(Particle &particle)
{
   cout<<"e_KN_loss"<<endl;
   cout<<"particle name="<<particle.name<<endl;

   int stat=0;
   float *ISRF_w;
   float *ISRF_over_nu;
   float *e_loss_compton_matrix;

   int ir,ix,iy,iz,ip;
   float sum;
   int i_comp;

   cout<<"e_KN_loss: generating ISRF_w and ISRF_over_nu"<<endl;
   ISRF_w      =new float[galaxy.ISRF[0].n_pgrid];            // an array of photon energies
   ISRF_over_nu=new float[galaxy.ISRF[0].n_pgrid];

   e_loss_compton_matrix=new float[particle.n_pgrid*galaxy.ISRF[0].n_pgrid];
   for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)          // target photon energy in mc^2 units
      ISRF_w[inu]=h_planck*galaxy.nu_ISRF[inu]* erg_to_eV * 1.0e-6 / m_electron; 
 
   int i_mat=0;
   for(ip=0; ip<particle.n_pgrid; ip++)
   {
      for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
      {
         e_loss_compton_matrix[i_mat]=e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip]);
// cout<<" e_KN_loss matrix "<<e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip])<<" "
// <<e_loss_compton_matrix[i_mat]<<endl;
         i_mat++;
      }
   }

   double factor=(eV_to_erg / h_planck)  *  log(galaxy.nu_ISRF[1]/ galaxy.nu_ISRF[0] );

   if(particle.n_spatial_dimensions==2)
   {
      cout<<"generating KN losses for particle.n_spatial_dimensions="<<particle.n_spatial_dimensions<<endl;

      for(ir=0; ir<particle.n_rgrid; ir++)
      {
         for(iz=0; iz<particle.n_zgrid; iz++)
         {
            for(i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)
            {
// cout<<" "<<ir<<" "<<iz<<" "<<i_comp<<endl;
               for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
               {
                  ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d2[ir][iz].s[inu]/galaxy.nu_ISRF[inu];
//cout<<inu<<" "<<ISRF_over_nu[inu];
               }
               i_mat=0;
               for(ip=0; ip<particle.n_pgrid; ip++)
               {
                  sum=0.0;
                  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
                  {
// sum+=ISRF_over_nu[inu]*e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip]);
// cout<<" e_KN_loss "<<e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip])<<" "<<e_loss_compton_matrix[i_mat]<<endl;
                     sum+=ISRF_over_nu[inu]*e_loss_compton_matrix[i_mat];
                     i_mat++;
	          }
//cout<<" "<<ir<<" "<<iz<<" "<<i_comp<<" "<<particle.gamma[ip]<<" "<<sum*factor<<endl;
                  particle.dpdt.d2[ir][iz].s[ip]+=sum*factor;
               }  //  ip
            }  //  ISRF_components
         }  //  iz
      }  //  ir
   }  //  particle.n_spatial_dimensions==2

   if(particle.n_spatial_dimensions==3)
   {
      cout<<"generating KN losses for particle.n_spatial_dimensions="<<particle.n_spatial_dimensions<<endl;

      for(ix=0; ix<particle.n_xgrid; ix++)
      {
         for(iy=0; iy<particle.n_ygrid; iy++)
         {
            for(iz=0; iz<particle.n_zgrid; iz++)
            {
// cout<<" "<<ix<<" "<<iy<<" "<<iz<<endl;

               for(i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)
               {
//cout<<"e_KN_loss ix iy iz i_comp  "<<ix<<" "<<iy<<" "<<iz<<" "<<i_comp<<endl;

                  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
                  {
                     ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d3[ix][iy][iz].s[inu]/galaxy.nu_ISRF[inu];
//cout<<inu<<" "<<ISRF_over_nu[inu];
                  }
                  i_mat=0;
                  for(ip=0; ip<particle.n_pgrid; ip++)
                  {
                     sum=0.0;
                     for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
                     {
//sum+=ISRF_over_nu[inu]*e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip]);
//cout<<" e_KN_loss "<<e_loss_compton_cc(ISRF_w[inu],particle.gamma[ip])<<" "<<e_loss_compton_matrix[i_mat]<<endl;
                        sum+=ISRF_over_nu[inu]*e_loss_compton_matrix[i_mat];
                        i_mat++;
	             }
//cout<<" "<<ir<<" "<<iz<<" "<<i_comp<<" "<<particle.gamma[ip]<<" "<<sum*factor<<endl;
                     particle.dpdt.d3[ix][iy][iz].s[ip]+=sum*factor;
                  }  //  ip
               }  //  ISRF_components
            }  //  iz
         }  //  iy
      }  //  ix
   }  //  particle.n_spatial_dimensions==3

   if(galdef.verbose>=2)particle.dpdt.print();
   cout<<" <<<< e_KN_loss"<<endl;
   return stat;
}

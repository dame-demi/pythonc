/*******************************************************************************
* CONVERGENT SCIENCE CONFIDENTIAL                                              *
* All rights reserved.                                                         *
* All information contained herein is the property of Convergent Science.      *
* The intellectual and technical concepts contained herein are                 *
* proprietary to Convergent Science.                                           *
* Dissemination of this information or reproduction of this material           *
* is strictly forbidden unless prior written permission is obtained from       *
* Convergent Science.                                                          *
*******************************************************************************/
#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "user_header.h"

/**********************************************************************/
/*                                                                    */
/* Name: user_spray_inject_custom                                     */
/*                                                                    */
/* Description: user_spray_inject_custom allows the user to customize */
/* certain drop properties (velocity components, spatial distribution */
/* components, drop radius)                                           */
/*                                                                    */
/* Called when: user_spray_inject_custom_flag=1 in udf.in             */
/*                                                                    */
/* Variables passed in:   injector_index injector number              */
/*                        nozzle_index   nozzle number                */
/*                        velocity      injection velocity            */
/*                        csuba         area contraction coefficient  */
/*                        random1       random number                 */
/*                        random2       random number                 */
/*                        random_fill   random number                 */
/*                        random_angle  random number                 */
/*                                                                    */
/* Variables passed back: *vel_axi      drop velocity component       */
/*                        *vel_norm     drop velocity component       */
/*                        *vel_other    drop velocity component       */
/*                        *xdist        drop spatial component        */
/*                        *ydist        drop spatial component        */
/*                        *zdist        drop spatial component        */
/*                        *radius_set   drop radius                   */
/*                                                                    */
/* Subroutines called: none                                           */
/*                                                                    */
/* Variable use:                                                      */
/*                                                                    */
/*   velocity components (mm=0,1,2)                                   */
/*     uu[mm]=vel_axi*injector[n]->nozzle[m]->axi_vec[mm]+            */
/*            vel_norm*injector[n]->nozzle[m]->norm_vec[mm]+          */
/*            vel_other*injector[n]->nozzle[m]->other_vec[mm]         */
/*                                                                    */
/*   spatial components                                               */
/*     xx[0]=injector[n]->nozzle[m]->x_noz+xdist;                     */
/*     xx[1]=injector[n]->nozzle[m]->y_noz+ydist;                     */
/*     xx[2]=injector[n]->nozzle[m]->z_noz+zdist;                     */
/*                                                                    */
/* Note: Random numbers are passed in for quality control checks.     */
/*       Use following syntax if additional random numbers are        */
/*       needed (where randomx is the random number):                 */
/*       random_real(&randomx);                                       */
/*                                                                    */
/**********************************************************************/


CONVERGE_UDF(spray_inject_custom,
             IN(VALUE(CONVERGE_index_t, passed_injector_index),
                VALUE(CONVERGE_index_t, passed_nozzle_index),
                VALUE(CONVERGE_precision_t, passed_velocity),
                VALUE(CONVERGE_precision_t, passed_csuba),
                VALUE(CONVERGE_precision_t, passed_random1),
                VALUE(CONVERGE_precision_t, passed_random2),
                VALUE(CONVERGE_precision_t, passed_random_fill),
                VALUE(CONVERGE_precision_t, passed_random_angle)),
             OUT(REF(CONVERGE_precision_t, passed_vel_axi),
                 REF(CONVERGE_precision_t, passed_vel_norm),
                 REF(CONVERGE_precision_t, passed_vel_other),
                 REF(CONVERGE_precision_t, passed_xdist),
                 REF(CONVERGE_precision_t, passed_ydist),
                 REF(CONVERGE_precision_t, passed_zdist),
                 REF(CONVERGE_precision_t, passed_smd_distribution),
                 REF(CONVERGE_precision_t, passed_radius_set)))
{

	int rank;
	CONVERGE_mpi_comm_rank(&rank);

   CONVERGE_precision_t spray_angle1, spray_angle2;
   CONVERGE_precision_t rad_local, angle_fill;
   CONVERGE_precision_t vel_local;
   CONVERGE_vec3_t norm_vec;
   CONVERGE_vec3_t other_vec;

   CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(passed_injector_index);
   CONVERGE_nozzle_t nozzle     = CONVERGE_injector_get_nozzle_with_id(injector, passed_nozzle_index);

   /**********************************************************************/
   /* this example implements a standard solid cone injection            */
   /**********************************************************************/

   // calculate the spray angles to be used in setting the velocity components
   spray_angle1 = CONVERGE_nozzle_get_cone_angle(nozzle);
   spray_angle1 = spray_angle1 * (fabs(passed_random1 - 0.5));
   spray_angle2 = 2.0 * M_PI * passed_random2;

   // velocity can be modified from passed in value through vel_local if desired
   vel_local = passed_velocity;

   // calculate the velocity components
   *passed_vel_axi   = (vel_local)*cos(spray_angle1);
   *passed_vel_norm  = (vel_local)*sin(spray_angle1) * cos(spray_angle2);
   *passed_vel_other = (vel_local)*sin(spray_angle1) * sin(spray_angle2);

   // calculate the position of the parcel relative to the nozzle center
   rad_local  = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_RADIUS_INJECT);
   angle_fill = 2.0 * M_PI * passed_random_angle;
   CONVERGE_nozzle_get_other_vec(nozzle, &other_vec);
   CONVERGE_nozzle_get_normal_vec(nozzle, &norm_vec);

   *passed_xdist = passed_random_fill * rad_local * (norm_vec[0] * cos(angle_fill) + other_vec[0] * sin(angle_fill));
   *passed_ydist = passed_random_fill * rad_local * (norm_vec[1] * cos(angle_fill) + other_vec[1] * sin(angle_fill));
   *passed_zdist = passed_random_fill * rad_local * (norm_vec[2] * cos(angle_fill) + other_vec[2] * sin(angle_fill));

//printf ("spray_inject_custom.c ambient_pres = %g ambient_temp = %g MW_Fuel = %g \n", ambient_pres, ambient_temp, MW_Fuel);
//printf ("spray_inject_custom.c den_liq = %g pvap_liq = %g  hvap_liq = %g fuel_temp = %g\n", fuel_den, fuel_pvap, fuel_hvap, fuel_temp);
   //printf ("spray_inject_custom.c sigma = %g vap_den = %g\n", sigma, fuel_den);

//Calculate reduced droplet diameter ---- Eqn references from Star CCM+ Manual on Flash Boiling

CONVERGE_precision_t theta_c, y, C_c0, r_D, T_s, rho_star, D_b, bub_freq, Nstar_nuc;
CONVERGE_precision_t F_rhostar, N_nuc, S_noz, V_b, vdot_vap, C_c, v_mean, v_vena, v_eff, vdot_tot;
CONVERGE_precision_t drop_radius_fact, g;

theta_c = 45.78; //bubble surface contact angle
y = 4.4; //eqn 2895 exponent
C_c0 = 0.611; //baseline contraction coefficient
r_D = 0.07; //roundness ratio at nozzle inlet
g = 9.81;
T_s = 340; //t_sat
 
 rho_star = (fuel_den - vap_den) / vap_den;
 D_b = 2.64E-5 * theta_c * pow((sigma/(g*(fuel_den - vap_den))),0.5) * pow(rho_star,0.9); //eqn 2893    
//Rewrite bub_freq and Nstar Nuc
 bub_freq = (1.18/D_b)*pow(((sigma*g*(fuel_den - vap_den))/pow(fuel_den,2)),0.25); //eqn 2892
 Nstar_nuc = pow((D_b*(fuel_temp - T_s)*vap_den*fuel_hvap/(2*sigma*T_s)),y); //eqn 2895
 
 F_rhostar = 2.157E-07 * pow(rho_star,-3.12) * pow((1 + 0.0049*rho_star),4.13);
 N_nuc = (1/pow(D_b,2)) * Nstar_nuc * F_rhostar; //eqn 2894
 S_noz = noz_dia * noz_len * M_PI; //inner surface area of nozzle
 V_b = pow(D_b,3) * 0.33333 * M_PI; //bubble volume at departure
 vdot_vap = bub_freq * V_b * S_noz * N_nuc; //eqn 2891
//vdot_vap = 2.5E-05;

//rewrite C_c
 C_c = pow((1/(pow(C_c0,2))-11.4*r_D),-0.5); //eqn 2889
  
  v_mean = mdot / (noz_area * fuel_den);
  v_vena = v_mean / C_c; //eqn 2888
  v_eff = (noz_area/mdot) * (fuel_pvap - ambient_pres) + v_vena; //eqn 2887
  vdot_tot = v_vena * noz_area; //calculate total vol flowrate in nozzle
  drop_radius_fact = (vdot_tot-vdot_vap)/vdot_tot;
//  drop_count = (1/drop_radius_fact)^3;


   //printf("spray_inject_custom.c original radius = %d\n", *passed_radius_set);
   // calculate the drop radius
   *passed_radius_set = *passed_radius_set*drop_radius_fact;
   printf("spray_inject_custom.c reduction factor = %g\n", drop_radius_fact);
   //printf("spray_inject_custom.c rho star = %g D_b = %g bub_freq = %g Nstar_nuc = %g\n", rho_star, D_b, bub_freq, Nstar_nuc);

   //*passed_radius_set = *passed_radius_set;
   // reset drop SMD if desired
   *passed_smd_distribution = *passed_smd_distribution;

   return;
}

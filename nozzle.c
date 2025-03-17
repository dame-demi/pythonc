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
/* Name: nozzle                                                       */
/*                                                                    */
/* Description: user_nozzle allows the user to customize the nozzle   */
/* injection calculation                                              */
/*                                                                    */
/* Called when: user_nozzle_flag=1 in udf.in                          */
/*                                                                    */
/* Variables passed in:   passed_injector_index injector index        */
/*                        passed_nozzle_index   nozzle number         */
/*                                                                    */
/* Variables passed back: *vel_eff     effective injection velocity   */
/*                        *diam_eff    effective nozzle diameter      */
/*                        *angle_eff   effective spray angle (degrees)*/
/*                                                                    */
/* Subroutines called: none                                           */
/*                                                                    */
/**********************************************************************/

CONVERGE_UDF(nozzle,
             IN(VALUE(CONVERGE_index_t, passed_injector_index),
                VALUE(CONVERGE_index_t, passed_nozzle_index),
                VALUE(CONVERGE_species_t, passed_species)),
             OUT(REF(CONVERGE_precision_t, passed_velocity_eff),
                 REF(CONVERGE_precision_t, passed_diameter_eff),
                 REF(CONVERGE_precision_t, passed_cone_angle_eff)))
{
   printf ("This is the start of nozzle.c \n");
	int rank;
	CONVERGE_mpi_comm_rank(&rank);

   CONVERGE_index_t num_parcel_species = CONVERGE_species_num_parcel(passed_species);

   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(passed_species, &parcel_species_it);
   CONVERGE_table_t *rho_table  = NULL;
   CONVERGE_table_t *visc_table = NULL;
   CONVERGE_table_t *surf_table = NULL;
   CONVERGE_table_t *pvap_table = NULL;
   CONVERGE_table_t *hvap_table = NULL;
   CONVERGE_table_t *cond_table = NULL;

   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_DENSITY_ID, &rho_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_VISC_ID, &visc_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_SURFACE_TENSION_ID, &surf_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_PVAP_ID, &pvap_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_HVAP_ID, &hvap_table);
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_COND_ID, &cond_table);

   CONVERGE_precision_t den_liq  = 0.0;
   CONVERGE_precision_t visc_liq = 0.0;
   CONVERGE_precision_t surf_liq = 0.0;
   CONVERGE_precision_t pvap_liq = 0.0;
   CONVERGE_precision_t hvap_liq = 0.0;
   CONVERGE_precision_t cond_liq = 0.0;

   CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(passed_injector_index);
   CONVERGE_nozzle_t nozzle     = CONVERGE_injector_get_nozzle_with_id(injector, passed_nozzle_index);
   CONVERGE_precision_t temp    = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_TEMP);
   CONVERGE_precision_t *mfrac  = malloc(sizeof(CONVERGE_precision_t));
   CONVERGE_injector_get_mass_fraction(injector, mfrac);

   for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
   {
      den_liq  = den_liq + mfrac[isp] * CONVERGE_table_lookup(rho_table[isp], temp);
      visc_liq = visc_liq + mfrac[isp] * CONVERGE_table_lookup(visc_table[isp], temp);
      surf_liq = surf_liq + mfrac[isp] * CONVERGE_table_lookup(surf_table[isp], temp);
      pvap_liq = pvap_liq + mfrac[isp] * CONVERGE_table_lookup(pvap_table[isp], temp);
      hvap_liq = hvap_liq + mfrac[isp] * CONVERGE_table_lookup(hvap_table[isp], temp);
      cond_liq = cond_liq + mfrac[isp] * CONVERGE_table_lookup(cond_table[isp], temp);
   }

   *passed_diameter_eff = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_DIAMETER);

   // following assumes that c_v=1.0

   CONVERGE_precision_t area_coeff = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_DISCHARGE_COEFF);
   CONVERGE_precision_t injected_mass = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_MASS); // line by James
   CONVERGE_precision_t injector_duration = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_DURATION); // line by James
   mdot = injected_mass/injector_duration; // line by James

   if(area_coeff > 1.0)
      area_coeff = 1.0;
   *passed_velocity_eff = (*passed_velocity_eff) / area_coeff;

   CONVERGE_precision_t area_noz = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_AREA);
   CONVERGE_precision_t noz_len = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_LENGTH); // line by James

   *passed_diameter_eff = 2.0 * sqrt((area_noz * area_coeff) / PI);

   printf ("This is from nozzle.c \n");
   printf ("nozzle.c rank = %d ambient_pres = %g ambient_temp = %g MW_Fuel = %g \n", rank, ambient_pres, ambient_temp, MW_Fuel);
   printf ("nozzle.c den_liq = %g pvap_liq = %g  hvap_liq = %g surf_liq = %g\n", den_liq, pvap_liq, hvap_liq, surf_liq);
   //CONVERGE_precision_t cone_angle = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_CONE_ANGLE);

   // Calculate increased cone angle due to flashing - From Price model. Equation references from Star CCM+ manual

   CONVERGE_precision_t a_0, v_l, big_theta, k_b, R_p, beta, jam_a, jam_b, jam_c, u, bulk_mod, c_vel;

   jam_a = -3.208; //coeff a in eqn 2890
   jam_b = 366.61; //coeff b in eqn 2890
   jam_c = -10324; //coeff c in eqn 2890
   u = 1.660539066E-27; //atomic mass unit
   v_l = 2.977E-28; //liquid molecular volume
   k_b = 1.380649E-23; //boltzmann constant
   bulk_mod = 1.9E9; //bulk modulus of liquid fuel

   a_0 = pow((36*M_PI), 0.333333) * pow(v_l, 0.666666); //Calculate molecular area
   big_theta = (a_0 * surf_liq) / (k_b * temp);
   R_p = pvap_liq / ambient_pres;
   beta = log10(pow(R_p, 2)*pow(big_theta, 3) / pow(MW_Fuel*u, 2));
   c_vel = sqrt(bulk_mod / den_liq); //calculate speed of sound in liquid fuel

//Rename variables to avoid problems when broadcasting
   fuel_den = den_liq;
   fuel_pvap = pvap_liq;
   fuel_hvap = hvap_liq;
   fuel_temp = temp;
   fuel_cond = cond_liq;
   fuel_visc = visc_liq;
   noz_dia = *passed_diameter_eff;
   noz_area = area_noz;
   sigma = surf_liq;
   sound_speed = c_vel;

//Write variables to header file for use in other routines
   CONVERGE_mpi_bcast(&fuel_den,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&fuel_hvap,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&fuel_cond,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&fuel_pvap,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&fuel_temp,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&fuel_visc,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&noz_dia,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&noz_area,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&noz_len,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&mdot,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&sigma,1,CONVERGE_DOUBLE,0);
   CONVERGE_mpi_bcast(&sound_speed,1,CONVERGE_DOUBLE,0);

//Pass cone angle to main routine
   *passed_cone_angle_eff = jam_a*pow(beta, 2) + jam_b*beta + jam_c; //calculate spray cone angle (eqn 2890 StarCCM+)

   //printf("Theta = %f\n", *passed_cone_angle_eff);
   unload_species_tables(parcel_species_it, &rho_table);
   unload_species_tables(parcel_species_it, &visc_table);
   unload_species_tables(parcel_species_it, &surf_table);
   unload_species_tables(parcel_species_it, &pvap_table);
   unload_species_tables(parcel_species_it, &hvap_table);
   unload_species_tables(parcel_species_it, &cond_table);

   CONVERGE_iterator_destroy(&parcel_species_it);

   return;
}

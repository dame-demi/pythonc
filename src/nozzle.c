#include "lagrangian/env.h"

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "user_header.h"

CONVERGE_UDF(nozzle,
             IN(VALUE(CONVERGE_index_t, passed_injector_index),
                VALUE(CONVERGE_index_t, passed_nozzle_index),
                VALUE(CONVERGE_species_t, passed_species)),
             OUT(REF(CONVERGE_precision_t, passed_velocity_eff),
                 REF(CONVERGE_precision_t, passed_diameter_eff),
                 REF(CONVERGE_precision_t, passed_cone_angle_eff)))
{
    int rank;
    CONVERGE_mpi_comm_rank(&rank);
    // printf("Rank: %d\n", rank);

    CONVERGE_index_t num_parcel_species = CONVERGE_species_num_parcel(passed_species);
    // printf("Number of parcel species: %d\n", num_parcel_species);

    CONVERGE_iterator_t parcel_species_it;
    CONVERGE_species_parcel_iterator_create(passed_species, &parcel_species_it);
    // printf("Parcel species iterator created.\n");

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
    // printf("Species tables loaded.\n");

    CONVERGE_precision_t den_liq  = 0.0;
    CONVERGE_precision_t visc_liq = 0.0;
    CONVERGE_precision_t surf_liq = 0.0;
    CONVERGE_precision_t pvap_liq = 0.0;
    CONVERGE_precision_t hvap_liq = 0.0;
    CONVERGE_precision_t cond_liq = 0.0;

    CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(passed_injector_index);
    CONVERGE_nozzle_t nozzle     = CONVERGE_injector_get_nozzle_with_id(injector, passed_nozzle_index);
    CONVERGE_precision_t temp    = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_TEMP);
    // printf("Injector and nozzle parameters retrieved.\n");

    CONVERGE_precision_t *mfrac  = (CONVERGE_precision_t *)malloc(sizeof(CONVERGE_precision_t) * num_parcel_species);
    if (mfrac == NULL) {
       // printf("Memory allocation for mfrac failed.\n");
        return;
    }
    CONVERGE_injector_get_mass_fraction(injector, mfrac);
    // printf("Mass fractions retrieved.\n");

    for (CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++) {
        den_liq  = den_liq + mfrac[isp] * CONVERGE_table_lookup(rho_table[isp], temp);
        visc_liq = visc_liq + mfrac[isp] * CONVERGE_table_lookup(visc_table[isp], temp);
        surf_liq = surf_liq + mfrac[isp] * CONVERGE_table_lookup(surf_table[isp], temp);
        pvap_liq = pvap_liq + mfrac[isp] * CONVERGE_table_lookup(pvap_table[isp], temp);
        hvap_liq = hvap_liq + mfrac[isp] * CONVERGE_table_lookup(hvap_table[isp], temp);
        cond_liq = cond_liq + mfrac[isp] * CONVERGE_table_lookup(cond_table[isp], temp);
    }
    // printf("Liquid properties calculated.\n");

    *passed_diameter_eff = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_DIAMETER);
    // printf("Effective diameter: %f\n", *passed_diameter_eff);

    CONVERGE_precision_t area_coeff = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_DISCHARGE_COEFF);
    CONVERGE_precision_t injected_mass = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_MASS);
    CONVERGE_precision_t injector_duration = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_INJECT_DURATION);
    CONVERGE_precision_t mdot = injected_mass / injector_duration;
    if (area_coeff > 1.0)
        area_coeff = 1.0;
    *passed_velocity_eff = (*passed_velocity_eff) / area_coeff;
    // printf("Effective velocity: %f\n", *passed_velocity_eff);

    CONVERGE_precision_t area_noz = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_AREA);
    CONVERGE_precision_t noz_len = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_LENGTH);
    *passed_diameter_eff = 2.0 * sqrt((area_noz * area_coeff) / M_PI);
    // printf("Nozzle area: %.12f, Nozzle length: %f\n", area_noz, noz_len);

    CONVERGE_precision_t v_l, big_theta, k_b, R_p, beta, jam_a, jam_b, jam_c, u;
    double a_0;

    //jam_a = -3.208;
    //jam_b = 366.61;
    //jam_c = -10324;
    u = 1.660539066E-27;
    v_l = 2.977E-28;
    k_b = 1.380649E-23;

    // a_0 = pow((36 * M_PI), 0.333333) * pow(v_l, 0.666666);
    a_0 = pow((36 * M_PI), 0.333333) * pow(v_l, 0.666666);
    // printf("v_l: %.12e\n", v_l);
    // printf("Intermediate value: %.12e\n", pow((36 * M_PI), 0.333333));
   // printf("a_0: %.12e\n", a_0);
    big_theta = (a_0 * surf_liq) / (k_b * temp);
    // printf("big_theta: %.12e\n", big_theta);
    R_p = pvap_liq / ambient_pres;
    // printf("R_p: %.12e\n", R_p);
    beta = log10(pow(R_p, 2) * pow(big_theta, 3) / pow(MW_Fuel * u, 2));
    // printf("beta: %.12e\n", beta);
    // printf("Cone angle parameters calculated.\n");
    // printf("Molar weight of the fuel: %.12e\n", MW_Fuel);

    fuel_den = den_liq;
    fuel_pvap = pvap_liq;
    fuel_hvap = hvap_liq;
    fuel_temp = temp;
    fuel_cond = cond_liq; // Not used in this function, but would be broadcast for later use
    noz_dia = *passed_diameter_eff;
    noz_area = area_noz;
    sigma = surf_liq;
    sound_speed = 1050; // default sound speed for IC3H18, will be overwritten if needed

    CONVERGE_mpi_bcast(&fuel_den, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&fuel_hvap, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&fuel_pvap, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&fuel_temp, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&fuel_cond, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&noz_dia, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&noz_area, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&noz_len, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&mdot, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&sigma, 1, CONVERGE_DOUBLE, 0);
    CONVERGE_mpi_bcast(&sound_speed, 1, CONVERGE_DOUBLE, 0);
   // printf("MPI broadcast completed.\n");

    CONVERGE_precision_t cone_angle = CONVERGE_nozzle_get_parameter_precision(nozzle, NOZZLE_CONE_ANGLE);

    *passed_cone_angle_eff = cone_angle * (180.0 / PI);

    //*passed_cone_angle_eff = jam_a * pow(beta, 2) + jam_b * beta + jam_c;
   // printf("Effective cone angle: %f\n", *passed_cone_angle_eff);

    unload_species_tables(parcel_species_it, &rho_table);
    unload_species_tables(parcel_species_it, &visc_table);
    unload_species_tables(parcel_species_it, &surf_table);
    unload_species_tables(parcel_species_it, &pvap_table);
    unload_species_tables(parcel_species_it, &hvap_table);
    unload_species_tables(parcel_species_it, &cond_table);
    // printf("Species tables unloaded.\n");

    CONVERGE_iterator_destroy(&parcel_species_it);
    // printf("Parcel species iterator destroyed.\n");

    free(mfrac);
   // printf("Memory freed.\n");

    return;
}
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

int USER_LAG_VAR;
   int USER_LAG_VARi;
   int USER_LAG_VARv3;
   int USER_LAG_VARv3b;

   // Lagrangian Data Field IDs
   int LAGRANGIAN_REL_VEL;
   int LAGRANGIAN_UPRIME;
   int LAGRANGIAN_UU;
   int LAGRANGIAN_ISBIG;
   int LAGRANGIAN_FROM_INJECTOR;
   int LAGRANGIAN_V_NU;
   int LAGRANGIAN_V_SH;
   int LAGRANGIAN_TEMP;
   int LAGRANGIAN_TEMP_TM1;
   int LAGRANGIAN_TEMP_STARM1;
   int LAGRANGIAN_REY_NUM;
   int LAGRANGIAN_REL_VEL_MAG;
   int LAGRANGIAN_RADIUS;
   int LAGRANGIAN_RADIUS_TM1;
   int LAGRANGIAN_PARENT;
   int LAGRANGIAN_DENSITY;
   int LAGRANGIAN_DENSITY_TM1;
   int LAGRANGIAN_MFRAC;
   int LAGRANGIAN_MFRAC_TM1;
   int LAGRANGIAN_NUM_DROP;
   int LAGRANGIAN_SURF_TEN;
   int LAGRANGIAN_VISCOSITY;
   int LAGRANGIAN_DISTORT;
   int LAGRANGIAN_DISTORT_DOT;
   int LAGRANGIAN_SURF_TEMP;
   int LAGRANGIAN_TBREAK_KH;
   int LAGRANGIAN_SHED_NUM_DROP;
   int LAGRANGIAN_SHED_MASS;
   int LAGRANGIAN_SACTIVE;
   int LAGRANGIAN_SACTIVE_TM1;
   int LAGRANGIAN_DM_DT;
   int LAGRANGIAN_DRDT;
   int LAGRANGIAN_WALL_HEAT_EXCHANGE;
   int LAGRANGIAN_L_RR;
   int LAGRANGIAN_L_RC;
   int LAGRANGIAN_L_TEMP1;
   int LAGRANGIAN_L_TEMP2;
   int LAGRANGIAN_JUST_HIT;
   int LAGRANGIAN_JUST_HIT_LEIDEN;
   int LAGRANGIAN_ISTHICK;
   int LAGRANGIAN_FILM_SHED;
   int LAGRANGIAN_XX;
   int LAGRANGIAN_ON_TRIANGLE;
   int LAGRANGIAN_FROM_NOZZLE;
   int LAGRANGIAN_TBREAK_RT;
   int LAGRANGIAN_SURF_TEMP_TM1;
   int LAGRANGIAN_NUM_DROP_TM1;
   int LAGRANGIAN_FILM_FLAG;
   int LAGRANGIAN_FILM_ENERGY;
   int LAGRANGIAN_T_TURB;
   int LAGRANGIAN_T_TURB_ACCUM;
   int LAGRANGIAN_FILM_THICKNESS;
   int LAGRANGIAN_AREA_IN_FILM;
   int LAGRANGIAN_FILM_THICKNESS_TM1;
   int LAGRANGIAN_AREA_REDUCTION;
   int LAGRANGIAN_TKE0;
   int LAGRANGIAN_EPS0;
   int LAGRANGIAN_LIFETIME;
   int LAGRANGIAN_UU_TM1;
   int LAGRANGIAN_FORCE_COEFFICIENT;
   int LAGRANGIAN_DROP_GAS_SRC;
   int LAGRANGIAN_GAS_DENSITY;
   int LAGRANGIAN_FILM_ACCUM_BIT_FLAG;
   int LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG;

// Nozzle Parameter IDs
int NOZZLE_AXIAL_VEC;
int NOZZLE_NORMAL_VEC;
int NOZZLE_OTHER_VEC;
int NOZZLE_AREA;
int NOZZLE_DIAMETER;
int NOZZLE_LENGTH;
int NOZZLE_SMD_DISTRIBUTION;
int NOZZLE_AMP_DISTORT;
int NOZZLE_CONE_ANGLE;
int NOZZLE_THICKNESS;
int NOZZLE_RADIAL_DISTANCE;
int NOZZLE_AXIAL_DISTANCE;
int NOZZLE_POSITION_VEC;
int NOZZLE_XX_VEC;
int NOZZLE_YY_VEC;
int NOZZLE_ZZ_VEC;
int NOZZLE_POSITION;
int NOZZLE_X;
int NOZZLE_Y;
int NOZZLE_Z;
int NOZZLE_MAX_DISTANCE;
int NOZZLE_PENETRATION;
int NOZZLE_VAPOR_PENETRATION;
int NOZZLE_ECN_PENETRATION;
int NOZZLE_TOTAL_MASS;
int NOZZLE_RADIUS_INJECT;
int NOZZLE_VOL_COUNT;
int NOZZLE_CONE_INPUT_FLAG;
int NOZZLE_DIAMETER_INPUT_FLAG;
int NOZZLE_RADIUS_INJECT_INPUT_FLAG;
int NOZZLE_NOZ_THETA;
int NOZZLE_NOZ_ANGLE_XY;
int NOZZLE_NOZ_ANGLE_XZ;

// Injector Parameter IDs
int INJECTOR_MFRAC;
int INJECTOR_START_INJECT_IS_FILE;
int INJECTOR_DUR_INJECT_IS_FILE;
int INJECTOR_MASS_INJECT_IS_FILE;
int INJECTOR_NOZZLE_INIT_FLAG;
int INJECTOR_SPRAY_INJECT_BC_FLAG;
int INJECTOR_ELSA_FLAG;
int INJECTOR_TEMP_FLAG;
int INJECTOR_TKE_FLAG;
int INJECTOR_EPS_FLAG;
int INJECTOR_INIT_CELL_TURB_FLAG;
int INJECTOR_KH_FLAG;
int INJECTOR_KHACT_NOZZLE_FLOW_FLAG;
int INJECTOR_KH_NEW_PARCEL_FLAG;
int INJECTOR_KH_NO_ENLARGE_FLAG;
int INJECTOR_INJECT_DISTRIBUTION_FLAG;
int INJECTOR_RT_FLAG;
int INJECTOR_RT_DISTRIBUTION_FLAG;
int INJECTOR_TAB_FLAG;
int INJECTOR_TAB_DISTRIBUTION_FLAG;
int INJECTOR_LISA_FLAG;
int INJECTOR_DISCHARGE_COEFF_INPUT_FLAG;
int INJECTOR_LISA_DISTRIBUTION_FLAG;
int INJECTOR_RATE_SHAPE_INPUT_FLAG;
int INJECTOR_DISCHARGE_COEFF_FLAG;
int INJECTOR_INJECT_INDEX;
int INJECTOR_NUM_PARCELS_PER_NOZZLE;
int INJECTOR_NUM_NOZZLES;
int INJECTOR_STREAM_INDEX;
int INJECTOR_NUM_INJECT_PARCELS;
int INJECTOR_CONE_FLAG;
int INJECTOR_TEMPORAL_TYPE;
int INJECTOR_VELOCITY_COEFF_INPUT_FLAG;
int INJECTOR_BREAKUP_FLAG;
int INJECTOR_POLAR_CPY_NUM;
int INJECTOR_Q_RR;
int INJECTOR_GAMMA_RR_X_DIST;
int INJECTOR_LISA_LENGTH_CONST;
int INJECTOR_LISA_SIZE_CONST;
int INJECTOR_LISA_INJECTION_PRES;
int INJECTOR_LISA_KV;
int INJECTOR_CYCLIC_PERIOD;
int INJECTOR_INJECT_START_TIME;
int INJECTOR_INJECT_DURATION;
int INJECTOR_END_INJECT;
int INJECTOR_INJECT_MASS;
int INJECTOR_INJECT_TEMP;
int INJECTOR_AREA;
int INJECTOR_KH_NEW_PARCEL_CUTOFF;
int INJECTOR_KH_SHED_FACTOR;
int INJECTOR_KH_BALPHA;
int INJECTOR_KH_CONST1;
int INJECTOR_KH_CONST2;
int INJECTOR_RT_LENGTH_CONST;
int INJECTOR_RT_CONST3;
int INJECTOR_RT_CONST2;
int INJECTOR_KHACT_TURB_KC;
int INJECTOR_KHACT_TURB_KE;
int INJECTOR_KHACT_TURB_S;
int INJECTOR_KHACT_C_TCAV;
int INJECTOR_TIME_OFFSET;
int INJECTOR_MASS_PER_PARCEL;
int INJECTOR_VELOCITY;
int INJECTOR_VELOCITY_TM1;
int INJECTOR_VELOCITY_OUT_OLD;
int INJECTOR_VELOCITY_OUT_NEW;
int INJECTOR_VELOCITY_COEFF;
int INJECTOR_MASS_TM1;
int INJECTOR_MASS_TM2;
int INJECTOR_DISCHARGE_COEFF;
int INJECTOR_RATE_SHAPE;
int INJECTOR_ANGLE_XY_INJ;
int INJECTOR_ANGLE_XZ_INJ;
int INJECTOR_RHO_LIQUID;
int INJECTOR_SWIRL_FRAC;
int INJECTOR_TOT_INJECTED_MASS;
int INJECTOR_TOT_INJECTED_MASS_TM;
int INJECTOR_AZIMUTH_ANGLE_START;
int INJECTOR_AZIMUTH_ANGLE_END;
int INJECTOR_SCALE_TEMP;
int INJECTOR_OFFSET_TEMP;
int INJECTOR_INJECT_TKE;
int INJECTOR_SCALE_TKE;
int INJECTOR_OFFSET_TKE;
int INJECTOR_INJECT_EPS;
int INJECTOR_SCALE_EPS;
int INJECTOR_OFFSET_EPS;
int INJECTOR_SWIRLER_MASS_FLOW_RATE;
int INJECTOR_SWIRLER_MEAN_ANGLE;
int INJECTOR_SWIRLER_RADIUS;
int INJECTOR_VOF_SPRAY_MASS_PER_PARCEL;
int INJECTOR_VOF_SPRAY_LIQ_FRAC_THRESHOLD;
int INJECTOR_X_CEN;
int INJECTOR_AXI_VEC;
int INJECTOR_NORM_VEC;
int INJECTOR_OTHER_VEC;

/** Load the spray environment
 */
CONVERGE_ONLOAD(spray_env, IN(CONVERGE_VOID))
{
   // Register a simple double data parcel field
   CONVERGE_variable_register("user_lag_var", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);

   // Register a simple int data parcel field
   CONVERGE_variable_register("user_lag_var_i", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);

   // User defined component names, overrides automatic nameing for CONVERGE_VEC3
   const char *user_lag_var_v3_comp_names[] = {"user_lag_var0", "user_lag_var1", "user_lag_var2"};
   CONVERGE_variable_register(
      "user_lag_var_v3",
      CONVERGE_DOUBLE,
      DEFAULT_PARCEL_VARIABLE_SETTINGS,
      // Double wtih dimension 3 may also be simply CONVERGE_VEC3.
      // The purpose of specifying the dimension manually is to allow for component names to be specified separately
      "dimension",
      3,
      // Pass the component names after the "component_names" parameter
      "component_names",
      user_lag_var_v3_comp_names,
      END_ARG_LIST);
   CONVERGE_variable_register(
      "user_lag_var_v3b",
      // CONVERGE_VEC3 will automattically append _1/_2/_3 to the end of the variable name
      CONVERGE_VEC3,
      DEFAULT_PARCEL_VARIABLE_SETTINGS,
      // The purpose of specifying the dimension manually here is to demonstrate it is permitted to mix dimension with any CONVERGE_APIType
      "dimension",
      3,
      END_ARG_LIST);
   // Get dynamic IDs to Lagrangian Cloud fields
   USER_LAG_VAR    = CONVERGE_lagrangian_field_id("user_lag_var");
   USER_LAG_VARi   = CONVERGE_lagrangian_field_id("user_lag_var_i");
   USER_LAG_VARv3  = CONVERGE_lagrangian_field_id("user_lag_var_v3");
   USER_LAG_VARv3b = CONVERGE_lagrangian_field_id("user_lag_var_v3b");

   LAGRANGIAN_FROM_INJECTOR = CONVERGE_lagrangian_field_id("LAGRANGIAN_FROM_INJECTOR");
   LAGRANGIAN_ON_TRIANGLE   = CONVERGE_lagrangian_field_id("LAGRANGIAN_ON_TRIANGLE");
   LAGRANGIAN_FROM_NOZZLE   = CONVERGE_lagrangian_field_id("LAGRANGIAN_FROM_NOZZLE");
   LAGRANGIAN_FILM_FLAG     = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_FLAG");

   LAGRANGIAN_JUST_HIT        = CONVERGE_lagrangian_field_id("LAGRANGIAN_JUST_HIT");
   LAGRANGIAN_JUST_HIT_LEIDEN = CONVERGE_lagrangian_field_id("LAGRANGIAN_JUST_HIT_LEIDEN");
   LAGRANGIAN_ISTHICK         = CONVERGE_lagrangian_field_id("LAGRANGIAN_ISTHICK");

   LAGRANGIAN_REL_VEL = CONVERGE_lagrangian_field_id("LAGRANGIAN_REL_VEL");
   LAGRANGIAN_UPRIME  = CONVERGE_lagrangian_field_id("LAGRANGIAN_UPRIME");
   LAGRANGIAN_UU      = CONVERGE_lagrangian_field_id("LAGRANGIAN_UU");
   LAGRANGIAN_UU_TM1  = CONVERGE_lagrangian_field_id("LAGRANGIAN_UU_TM1");
   LAGRANGIAN_XX      = CONVERGE_lagrangian_field_id("LAGRANGIAN_XX");

   LAGRANGIAN_V_NU               = CONVERGE_lagrangian_field_id("LAGRANGIAN_V_NU");
   LAGRANGIAN_V_SH               = CONVERGE_lagrangian_field_id("LAGRANGIAN_V_SH");
   LAGRANGIAN_TEMP               = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP");
   LAGRANGIAN_TEMP_TM1           = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP_TM1");
   LAGRANGIAN_TEMP_STARM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP_STARM1");
   LAGRANGIAN_REY_NUM            = CONVERGE_lagrangian_field_id("LAGRANGIAN_REY_NUM");
   LAGRANGIAN_REL_VEL_MAG        = CONVERGE_lagrangian_field_id("LAGRANGIAN_REL_VEL_MAG");
   LAGRANGIAN_RADIUS             = CONVERGE_lagrangian_field_id("LAGRANGIAN_RADIUS");
   LAGRANGIAN_RADIUS_TM1         = CONVERGE_lagrangian_field_id("LAGRANGIAN_RADIUS_TM1");
   LAGRANGIAN_PARENT             = CONVERGE_lagrangian_field_id("LAGRANGIAN_PARENT");
   LAGRANGIAN_DENSITY            = CONVERGE_lagrangian_field_id("LAGRANGIAN_DENSITY");
   LAGRANGIAN_DENSITY_TM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_DENSITY_TM1");
   LAGRANGIAN_GAS_DENSITY        = CONVERGE_lagrangian_field_id("LAGRANGIAN_GAS_DENSITY");
   LAGRANGIAN_MFRAC              = CONVERGE_lagrangian_field_id("LAGRANGIAN_MFRAC");
   LAGRANGIAN_MFRAC_TM1          = CONVERGE_lagrangian_field_id("LAGRANGIAN_MFRAC_TM1");
   LAGRANGIAN_NUM_DROP           = CONVERGE_lagrangian_field_id("LAGRANGIAN_NUM_DROP");
   LAGRANGIAN_SURF_TEMP          = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEMP");
   LAGRANGIAN_TBREAK_KH          = CONVERGE_lagrangian_field_id("LAGRANGIAN_TBREAK_KH");
   LAGRANGIAN_SHED_NUM_DROP      = CONVERGE_lagrangian_field_id("LAGRANGIAN_SHED_NUM_DROP");
   LAGRANGIAN_SHED_MASS          = CONVERGE_lagrangian_field_id("LAGRANGIAN_SHED_MASS");
   LAGRANGIAN_SACTIVE            = CONVERGE_lagrangian_field_id("LAGRANGIAN_SACTIVE");
   LAGRANGIAN_SACTIVE_TM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_SACTIVE_TM1");
   LAGRANGIAN_SURF_TEN           = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEN");
   LAGRANGIAN_VISCOSITY          = CONVERGE_lagrangian_field_id("LAGRANGIAN_VISCOSITY");
   LAGRANGIAN_DISTORT            = CONVERGE_lagrangian_field_id("LAGRANGIAN_DISTORT");
   LAGRANGIAN_DISTORT_DOT        = CONVERGE_lagrangian_field_id("LAGRANGIAN_DISTORT_DOT");
   LAGRANGIAN_DM_DT              = CONVERGE_lagrangian_field_id("LAGRANGIAN_DM_DT");
   LAGRANGIAN_DRDT               = CONVERGE_lagrangian_field_id("LAGRANGIAN_DRDT");
   LAGRANGIAN_WALL_HEAT_EXCHANGE = CONVERGE_lagrangian_field_id("LAGRANGIAN_WALL_HEAT_EXCHANGE");
   LAGRANGIAN_L_RR               = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_RR");
   LAGRANGIAN_L_RC               = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_RC");
   LAGRANGIAN_L_TEMP1            = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_TEMP1");
   LAGRANGIAN_L_TEMP2            = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_TEMP2");
   LAGRANGIAN_FILM_SHED          = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_SHED");
   LAGRANGIAN_TBREAK_RT          = CONVERGE_lagrangian_field_id("LAGRANGIAN_TBREAK_RT");
   LAGRANGIAN_SURF_TEMP_TM1      = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEMP_TM1");
   LAGRANGIAN_NUM_DROP_TM1       = CONVERGE_lagrangian_field_id("LAGRANGIAN_NUM_DROP_TM1");
   LAGRANGIAN_FILM_ENERGY        = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ENERGY");
   LAGRANGIAN_T_TURB             = CONVERGE_lagrangian_field_id("LAGRANGIAN_T_TURB");
   LAGRANGIAN_T_TURB_ACCUM       = CONVERGE_lagrangian_field_id("LAGRANGIAN_T_TURB_ACCUM");
   LAGRANGIAN_FILM_THICKNESS     = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_THICKNESS");
   LAGRANGIAN_AREA_REDUCTION     = CONVERGE_lagrangian_field_id("LAGRANGIAN_AREA_REDUCTION");
   LAGRANGIAN_TKE0               = CONVERGE_lagrangian_field_id("LAGRANGIAN_TKE0");
   LAGRANGIAN_EPS0               = CONVERGE_lagrangian_field_id("LAGRANGIAN_EPS0");
   LAGRANGIAN_LIFETIME           = CONVERGE_lagrangian_field_id("LAGRANGIAN_LIFETIME");
   LAGRANGIAN_FORCE_COEFFICIENT  = CONVERGE_lagrangian_field_id("LAGRANGIAN_FORCE_COEFFICIENT");
   LAGRANGIAN_DROP_GAS_SRC       = CONVERGE_lagrangian_field_id("LAGRANGIAN_DROP_GAS_SRC");

   LAGRANGIAN_FILM_THICKNESS_TM1 = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_THICKNESS_TM1");
   LAGRANGIAN_AREA_IN_FILM       = CONVERGE_lagrangian_field_id("LAGRANGIAN_AREA_IN_FILM");

   LAGRANGIAN_FILM_ACCUM_BIT_FLAG = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ACCUM_BIT_FLAG");
   LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG");

   // Get Dynamic IDs for Injectors parameters
   INJECTOR_MFRAC                        = CONVERGE_get_parameter_id("injector.mfrac");
   INJECTOR_START_INJECT_IS_FILE         = CONVERGE_get_parameter_id("injector.start_inject_is_file");
   INJECTOR_DUR_INJECT_IS_FILE           = CONVERGE_get_parameter_id("injector.dur_inject_is_file");
   INJECTOR_MASS_INJECT_IS_FILE          = CONVERGE_get_parameter_id("injector.mass_inject_is_file");
   INJECTOR_NOZZLE_INIT_FLAG             = CONVERGE_get_parameter_id("injector.nozzle_init_flag");
   INJECTOR_SPRAY_INJECT_BC_FLAG         = CONVERGE_get_parameter_id("injector.spray_inject_bc_flag");
   INJECTOR_ELSA_FLAG                    = CONVERGE_get_parameter_id("injector.elsa_flag");
   INJECTOR_TEMP_FLAG                    = CONVERGE_get_parameter_id("injector.temp_flag");
   INJECTOR_TKE_FLAG                     = CONVERGE_get_parameter_id("injector.tke_flag");
   INJECTOR_EPS_FLAG                     = CONVERGE_get_parameter_id("injector.eps_flag");
   INJECTOR_INIT_CELL_TURB_FLAG          = CONVERGE_get_parameter_id("injector.init_cell_turb_flag");
   INJECTOR_KH_FLAG                      = CONVERGE_get_parameter_id("injector.kh_flag");
   INJECTOR_KHACT_NOZZLE_FLOW_FLAG       = CONVERGE_get_parameter_id("injector.khact_nozzle_flow_flag");
   INJECTOR_KH_NEW_PARCEL_FLAG           = CONVERGE_get_parameter_id("injector.kh_new_parcel_flag");
   INJECTOR_KH_NO_ENLARGE_FLAG           = CONVERGE_get_parameter_id("injector.kh_no_enlarge_flag");
   INJECTOR_INJECT_DISTRIBUTION_FLAG     = CONVERGE_get_parameter_id("injector.inject_distribution_flag");
   INJECTOR_RT_FLAG                      = CONVERGE_get_parameter_id("injector.rt_flag");
   INJECTOR_RT_DISTRIBUTION_FLAG         = CONVERGE_get_parameter_id("injector.rt_distribution_flag");
   INJECTOR_TAB_FLAG                     = CONVERGE_get_parameter_id("injector.tab_flag");
   INJECTOR_TAB_DISTRIBUTION_FLAG        = CONVERGE_get_parameter_id("injector.tab_distribution_flag");
   INJECTOR_LISA_FLAG                    = CONVERGE_get_parameter_id("injector.lisa_flag");
   INJECTOR_DISCHARGE_COEFF_INPUT_FLAG   = CONVERGE_get_parameter_id("injector.discharge_coeff_input_flag");
   INJECTOR_LISA_DISTRIBUTION_FLAG       = CONVERGE_get_parameter_id("injector.lisa_distribution_flag");
   INJECTOR_RATE_SHAPE_INPUT_FLAG        = CONVERGE_get_parameter_id("injector.rate_shape_input_flag");
   INJECTOR_DISCHARGE_COEFF_FLAG         = CONVERGE_get_parameter_id("injector.discharge_coeff_flag");
   INJECTOR_INJECT_INDEX                 = CONVERGE_get_parameter_id("injector.inject_index");
   INJECTOR_NUM_PARCELS_PER_NOZZLE       = CONVERGE_get_parameter_id("injector.num_parcels_per_nozzle");
   INJECTOR_NUM_NOZZLES                  = CONVERGE_get_parameter_id("injector.num_nozzles");
   INJECTOR_STREAM_INDEX                 = CONVERGE_get_parameter_id("injector.stream_index");
   INJECTOR_NUM_INJECT_PARCELS           = CONVERGE_get_parameter_id("injector.num_inject_parcels");
   INJECTOR_CONE_FLAG                    = CONVERGE_get_parameter_id("injector.cone_flag");
   INJECTOR_TEMPORAL_TYPE                = CONVERGE_get_parameter_id("injector.temporal_type");
   INJECTOR_VELOCITY_COEFF_INPUT_FLAG    = CONVERGE_get_parameter_id("injector.velocity_coeff_input_flag");
   INJECTOR_BREAKUP_FLAG                 = CONVERGE_get_parameter_id("injector.breakup_flag");
   INJECTOR_POLAR_CPY_NUM                = CONVERGE_get_parameter_id("injector.polar_cpy_num");
   INJECTOR_Q_RR                         = CONVERGE_get_parameter_id("injector.q_rr");
   INJECTOR_GAMMA_RR_X_DIST              = CONVERGE_get_parameter_id("injector.gamma_rr_x_dist");
   INJECTOR_LISA_LENGTH_CONST            = CONVERGE_get_parameter_id("injector.lisa_length_const");
   INJECTOR_LISA_SIZE_CONST              = CONVERGE_get_parameter_id("injector.lisa_size_const");
   INJECTOR_LISA_INJECTION_PRES          = CONVERGE_get_parameter_id("injector.lisa_injection_pressure");
   INJECTOR_LISA_KV                      = CONVERGE_get_parameter_id("injector.lisa_kv");
   INJECTOR_CYCLIC_PERIOD                = CONVERGE_get_parameter_id("injector.cyclic_period");
   INJECTOR_INJECT_START_TIME            = CONVERGE_get_parameter_id("injector.inject_start_time");
   INJECTOR_INJECT_DURATION              = CONVERGE_get_parameter_id("injector.inject_duration");
   INJECTOR_END_INJECT                   = CONVERGE_get_parameter_id("injector.end_inject");
   INJECTOR_INJECT_MASS                  = CONVERGE_get_parameter_id("injector.inject_mass");
   INJECTOR_INJECT_TEMP                  = CONVERGE_get_parameter_id("injector.inject_temperature");
   INJECTOR_AREA                         = CONVERGE_get_parameter_id("injector.tot_area");
   INJECTOR_KH_NEW_PARCEL_CUTOFF         = CONVERGE_get_parameter_id("injector.kh_new_parcel_cutoff");
   INJECTOR_KH_SHED_FACTOR               = CONVERGE_get_parameter_id("injector.kh_shed_factor");
   INJECTOR_KH_BALPHA                    = CONVERGE_get_parameter_id("injector.kh_balpha");
   INJECTOR_KH_CONST1                    = CONVERGE_get_parameter_id("injector.kh_const1");
   INJECTOR_KH_CONST2                    = CONVERGE_get_parameter_id("injector.kh_const2");
   INJECTOR_RT_LENGTH_CONST              = CONVERGE_get_parameter_id("injector.rt_length_const");
   INJECTOR_RT_CONST3                    = CONVERGE_get_parameter_id("injector.rt_const3");
   INJECTOR_RT_CONST2                    = CONVERGE_get_parameter_id("injector.rt_const2");
   INJECTOR_KHACT_TURB_KC                = CONVERGE_get_parameter_id("injector.khact_turb_kc");
   INJECTOR_KHACT_TURB_KE                = CONVERGE_get_parameter_id("injector.khact_turb_ke");
   INJECTOR_KHACT_TURB_S                 = CONVERGE_get_parameter_id("injector.khact_turb_s");
   INJECTOR_KHACT_C_TCAV                 = CONVERGE_get_parameter_id("injector.khact_c_tcav");
   INJECTOR_TIME_OFFSET                  = CONVERGE_get_parameter_id("injector.time_offset");
   INJECTOR_MASS_PER_PARCEL              = CONVERGE_get_parameter_id("injector.mass_per_parcel");
   INJECTOR_VELOCITY                     = CONVERGE_get_parameter_id("injector.velocity");
   INJECTOR_VELOCITY_TM1                 = CONVERGE_get_parameter_id("injector.velocity_tm1");
   INJECTOR_VELOCITY_OUT_OLD             = CONVERGE_get_parameter_id("injector.velocity_out_old");
   INJECTOR_VELOCITY_OUT_NEW             = CONVERGE_get_parameter_id("injector.velocity_out_new");
   INJECTOR_VELOCITY_COEFF               = CONVERGE_get_parameter_id("injector.velocity_coeff");
   INJECTOR_MASS_TM1                     = CONVERGE_get_parameter_id("injector.mass_tm1");
   INJECTOR_MASS_TM2                     = CONVERGE_get_parameter_id("injector.mass_tm2");
   INJECTOR_DISCHARGE_COEFF              = CONVERGE_get_parameter_id("injector.discharge_coeff");
   INJECTOR_RATE_SHAPE                   = CONVERGE_get_parameter_id("injector.rate_shape");
   INJECTOR_ANGLE_XY_INJ                 = CONVERGE_get_parameter_id("injector.angle_xy_inj");
   INJECTOR_ANGLE_XZ_INJ                 = CONVERGE_get_parameter_id("injector.angle_xz_inj");
   INJECTOR_RHO_LIQUID                   = CONVERGE_get_parameter_id("injector.rho_liquid");
   INJECTOR_SWIRL_FRAC                   = CONVERGE_get_parameter_id("injector.swirl_frac");
   INJECTOR_TOT_INJECTED_MASS            = CONVERGE_get_parameter_id("injector.tot_inj_mass");
   INJECTOR_TOT_INJECTED_MASS_TM         = CONVERGE_get_parameter_id("injector.tot_inj_mass_tm");
   INJECTOR_AZIMUTH_ANGLE_START          = CONVERGE_get_parameter_id("injector.azimuth_angle_start");
   INJECTOR_AZIMUTH_ANGLE_END            = CONVERGE_get_parameter_id("injector.azimuth_angle_end");
   INJECTOR_SCALE_TEMP                   = CONVERGE_get_parameter_id("injector.scale_temperature");
   INJECTOR_OFFSET_TEMP                  = CONVERGE_get_parameter_id("injector.offset_temperature");
   INJECTOR_INJECT_TKE                   = CONVERGE_get_parameter_id("injector.inject_tke");
   INJECTOR_SCALE_TKE                    = CONVERGE_get_parameter_id("injector.scale_tke");
   INJECTOR_OFFSET_TKE                   = CONVERGE_get_parameter_id("injector.offset_tke");
   INJECTOR_INJECT_EPS                   = CONVERGE_get_parameter_id("injector.inject_epsilon");
   INJECTOR_SCALE_EPS                    = CONVERGE_get_parameter_id("injector.scale_epsilon");
   INJECTOR_OFFSET_EPS                   = CONVERGE_get_parameter_id("injector.offset_epsilon");
   INJECTOR_SWIRLER_MASS_FLOW_RATE       = CONVERGE_get_parameter_id("injector.swirler_mass_flow_rate");
   INJECTOR_SWIRLER_MEAN_ANGLE           = CONVERGE_get_parameter_id("injector.swirler_mean_angle");
   INJECTOR_SWIRLER_RADIUS               = CONVERGE_get_parameter_id("injector.swirler_radius");
   INJECTOR_VOF_SPRAY_MASS_PER_PARCEL    = CONVERGE_get_parameter_id("injector.vof_spray_mass_per_parcel");
   INJECTOR_VOF_SPRAY_LIQ_FRAC_THRESHOLD = CONVERGE_get_parameter_id("injector.vof_spray_liq_frac_threshold");
   INJECTOR_X_CEN                        = CONVERGE_get_parameter_id("injector.x_cen_inj");
   INJECTOR_AXI_VEC                      = CONVERGE_get_parameter_id("injector.axi_vec");
   INJECTOR_NORM_VEC                     = CONVERGE_get_parameter_id("injector.norm_vec");
   INJECTOR_OTHER_VEC                    = CONVERGE_get_parameter_id("injector.other_vec");

   // Get Dynamic IDs for Nozzle parameters
   NOZZLE_AXIAL_VEC                = CONVERGE_get_parameter_id("nozzle.axial_vec");
   NOZZLE_NORMAL_VEC               = CONVERGE_get_parameter_id("nozzle.normal_vec");
   NOZZLE_OTHER_VEC                = CONVERGE_get_parameter_id("nozzle.other_vec");
   NOZZLE_AREA                     = CONVERGE_get_parameter_id("nozzle.area_noz");
   NOZZLE_DIAMETER                 = CONVERGE_get_parameter_id("nozzle.diameter");
   NOZZLE_LENGTH                   = CONVERGE_get_parameter_id("nozzle.length");
   NOZZLE_SMD_DISTRIBUTION         = CONVERGE_get_parameter_id("nozzle.smd_distribution");
   NOZZLE_AMP_DISTORT              = CONVERGE_get_parameter_id("nozzle.amp_distort");
   NOZZLE_CONE_ANGLE               = CONVERGE_get_parameter_id("nozzle.cone_angle");
   NOZZLE_THICKNESS                = CONVERGE_get_parameter_id("nozzle.thickness");
   NOZZLE_RADIAL_DISTANCE          = CONVERGE_get_parameter_id("nozzle.radial_distance");
   NOZZLE_AXIAL_DISTANCE           = CONVERGE_get_parameter_id("nozzle.axial_distance");
   NOZZLE_POSITION_VEC             = CONVERGE_get_parameter_id("nozzle.position_vec");
   NOZZLE_XX_VEC                   = CONVERGE_get_parameter_id("nozzle.xx_vec");
   NOZZLE_YY_VEC                   = CONVERGE_get_parameter_id("nozzle.yy_vec");
   NOZZLE_ZZ_VEC                   = CONVERGE_get_parameter_id("nozzle.zz_vec");
   NOZZLE_POSITION                 = CONVERGE_get_parameter_id("nozzle.position");
   NOZZLE_X                        = CONVERGE_get_parameter_id("nozzle.x");
   NOZZLE_Y                        = CONVERGE_get_parameter_id("nozzle.y");
   NOZZLE_Z                        = CONVERGE_get_parameter_id("nozzle.z");
   NOZZLE_MAX_DISTANCE             = CONVERGE_get_parameter_id("nozzle.max_distance");
   NOZZLE_PENETRATION              = CONVERGE_get_parameter_id("nozzle.penetration");
   NOZZLE_VAPOR_PENETRATION        = CONVERGE_get_parameter_id("nozzle.vapor_penetration");
   NOZZLE_ECN_PENETRATION          = CONVERGE_get_parameter_id("nozzle.ecn_penetration");
   NOZZLE_TOTAL_MASS               = CONVERGE_get_parameter_id("nozzle.total_mass");
   NOZZLE_RADIUS_INJECT            = CONVERGE_get_parameter_id("nozzle.radius_inject");
   NOZZLE_VOL_COUNT                = CONVERGE_get_parameter_id("nozzle.vol_count");
   NOZZLE_CONE_INPUT_FLAG          = CONVERGE_get_parameter_id("nozzle.cone_noz_input_flag");
   NOZZLE_DIAMETER_INPUT_FLAG      = CONVERGE_get_parameter_id("nozzle.noz_diameter_input_flag");
   NOZZLE_RADIUS_INJECT_INPUT_FLAG = CONVERGE_get_parameter_id("nozzle.radius_inject_input_flag");
   NOZZLE_NOZ_THETA                = CONVERGE_get_parameter_id("nozzle.noz_theta");
   NOZZLE_NOZ_ANGLE_XY             = CONVERGE_get_parameter_id("nozzle.noz_angle_xy");
   NOZZLE_NOZ_ANGLE_XZ             = CONVERGE_get_parameter_id("nozzle.noz_angle_xz");
}

/** Load all of the cloud data from a CONVERGE_cloud_t to a wrapper structure
 */
void load_user_cloud(struct ParcelCloud *parcel_cloud_loc, CONVERGE_cloud_t c)
{
   // Load user defined field data
   parcel_cloud_loc->user_temp_starm1 = (double *)CONVERGE_cloud_get_field_data(c, USER_LAG_VAR);
   parcel_cloud_loc->user_lag_var_i   = (int *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARi);
   parcel_cloud_loc->user_lag_var_v3  = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3);
   parcel_cloud_loc->user_lag_var_v3b = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3b);

   parcel_cloud_loc->from_injector = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FROM_INJECTOR);
   parcel_cloud_loc->on_triangle   = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_ON_TRIANGLE);
   parcel_cloud_loc->from_nozzle   = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FROM_NOZZLE);
   parcel_cloud_loc->film_flag     = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_FLAG);

   parcel_cloud_loc->just_hit        = (short *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_JUST_HIT);
   parcel_cloud_loc->just_hit_leiden = (short *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_JUST_HIT_LEIDEN);
   parcel_cloud_loc->is_thick        = (char *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_ISTHICK);

   parcel_cloud_loc->rel_vel      = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REL_VEL);
   parcel_cloud_loc->uprime       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UPRIME);
   parcel_cloud_loc->uu           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UU);
   parcel_cloud_loc->uu_tm1       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UU_TM1);
   parcel_cloud_loc->xx           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_XX);
   parcel_cloud_loc->drop_gas_src = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DROP_GAS_SRC);

   parcel_cloud_loc->v_nu          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_V_NU);
   parcel_cloud_loc->v_sh          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_V_SH);
   parcel_cloud_loc->temp          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP);
   parcel_cloud_loc->temp_tm1      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP_TM1);
   parcel_cloud_loc->temp_starm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP_STARM1);
   parcel_cloud_loc->rey_num       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REY_NUM);
   parcel_cloud_loc->rel_vel_mag   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REL_VEL_MAG);
   parcel_cloud_loc->radius        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_RADIUS);
   parcel_cloud_loc->radius_tm1    = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_RADIUS_TM1);
   parcel_cloud_loc->parent        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_PARENT);
   parcel_cloud_loc->density       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DENSITY);
   parcel_cloud_loc->density_tm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DENSITY_TM1);
   parcel_cloud_loc->gas_density   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_GAS_DENSITY);
   parcel_cloud_loc->mfrac         = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_MFRAC);
   parcel_cloud_loc->mfrac_tm1     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_MFRAC_TM1);
   parcel_cloud_loc->num_drop      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_NUM_DROP);
   parcel_cloud_loc->surf_temp     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEMP);
   parcel_cloud_loc->tbreak_kh     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TBREAK_KH);
   parcel_cloud_loc->shed_num_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SHED_NUM_DROP);
   parcel_cloud_loc->shed_mass     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SHED_MASS);
   parcel_cloud_loc->sactive       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SACTIVE);
   parcel_cloud_loc->sactive_tm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SACTIVE_TM1);
   parcel_cloud_loc->surf_ten      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEN);
   parcel_cloud_loc->viscosity     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_VISCOSITY);
   parcel_cloud_loc->distort       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DISTORT);
   parcel_cloud_loc->distort_dot   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DISTORT_DOT);
   parcel_cloud_loc->dm_dt         = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DM_DT);
   parcel_cloud_loc->drdt          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DRDT);
   parcel_cloud_loc->wall_heat_exchange =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_WALL_HEAT_EXCHANGE);
   parcel_cloud_loc->l_rr          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_RR);
   parcel_cloud_loc->l_rc          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_RC);
   parcel_cloud_loc->l_temp1       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_TEMP1);
   parcel_cloud_loc->l_temp2       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_TEMP2);
   parcel_cloud_loc->film_shed     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_SHED);
   parcel_cloud_loc->tbreak_rt     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TBREAK_RT);
   parcel_cloud_loc->surf_temp_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEMP_TM1);
   parcel_cloud_loc->num_drop_tm1  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_NUM_DROP_TM1);
   parcel_cloud_loc->film_energy   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ENERGY);
   parcel_cloud_loc->t_turb        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_T_TURB);
   parcel_cloud_loc->t_turb_accum  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_T_TURB_ACCUM);
   parcel_cloud_loc->film_thickness =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_THICKNESS);
   parcel_cloud_loc->area_reduction =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_AREA_REDUCTION);
   parcel_cloud_loc->tke0     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TKE0);
   parcel_cloud_loc->eps0     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_EPS0);
   parcel_cloud_loc->lifetime = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_LIFETIME);
   parcel_cloud_loc->force_coefficient = 
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FORCE_COEFFICIENT);
   parcel_cloud_loc->film_accum_bit_flag = 
      (unsigned long long*)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ACCUM_BIT_FLAG);
   parcel_cloud_loc->film_accum_plus_bit_flag = 
      (unsigned int*)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG);

   if(CONVERGE_cloud_type(c) == LAGRANGIAN_FILM)
   {
      parcel_cloud_loc->film_thickness_tm1 =
         (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_THICKNESS_TM1);
      parcel_cloud_loc->area_in_film =
         (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_AREA_IN_FILM);
   }
   else
   {
      parcel_cloud_loc->film_thickness_tm1 = NULL;
      parcel_cloud_loc->area_in_film       = NULL;
   }
}

/* load urea information*/

void load_urea_parameters(struct UreaInfo *urea_info)
{
   urea_info->isp_nh3_urea             = CONVERGE_get_int("urea_in.isp_nh3_urea");
   urea_info->isp_hnco_urea            = CONVERGE_get_int("urea_in.isp_hnco_urea");
   urea_info->isp_co2_urea             = CONVERGE_get_int("urea_in.isp_co2_urea");
   urea_info->isp_cn_urea              = CONVERGE_get_int("urea_in.isp_cn_urea");
   urea_info->isp_nh2_urea             = CONVERGE_get_int("urea_in.isp_nh2_urea");
   urea_info->isp_water_parcel         = CONVERGE_get_int("urea_in.isp_water_parcel");
   urea_info->isp_urea_parcel          = CONVERGE_get_int("urea_in.isp_urea_parcel");
   urea_info->isp_urea_parcel_sl       = CONVERGE_get_int("urea_in.isp_urea_parcel_sl");
   urea_info->urea_a                   = CONVERGE_get_double("urea_in.urea_a");
   urea_info->urea_ea                  = CONVERGE_get_double("urea_in.urea_ea");
   urea_info->gamma_scr                = CONVERGE_get_double("urea_in.gamma_scr");
   urea_info->urea_hdcmp               = CONVERGE_get_double("urea_in.urea_hdcmp");
}

void load_urea_dd_parameters(struct UreaInfo *urea_info)
{
   urea_info->num_react_parcel_species = CONVERGE_get_int("urea_in.num_react_parcel_species");
   urea_info->num_react_gas_species    = CONVERGE_get_int("urea_in.num_react_gas_species");
   urea_info->num_reactions            = CONVERGE_get_int("urea_in.num_reactions");
   urea_info->isp_nh3                  = CONVERGE_get_int("urea_in.isp_nh3");//gas
   urea_info->isp_hnco                 = CONVERGE_get_int("urea_in.isp_hnco");//gas
   urea_info->isp_h2o_aq               = CONVERGE_get_int("urea_in.isp_h2o_aq");//aqueous
   urea_info->isp_urea_aq              = CONVERGE_get_int("urea_in.isp_urea_aq");//aqueous
   urea_info->isp_urea_sl              = CONVERGE_get_int("urea_in.isp_urea_sl");//solid
   urea_info->scaling_factor_a         = CONVERGE_get_double("urea_in.scaling_factor_a");
   urea_info->scaling_factor_e         = CONVERGE_get_double("urea_in.scaling_factor_e");
   urea_info->mw_h2o                   = CONVERGE_get_double("urea_in.mw_h2o");
   urea_info->mw_urea                  = CONVERGE_get_double("urea_in.mw_urea");

   /*CONVERGE_urea_get_isp_lists(&(urea_info->isp_parcel_species_list),
                              &(urea_info->isp_gas_species_list),
                              &(urea_info->isp_react_species_list));

   CONVERGE_urea_get_double_lists(&(urea_info->mass_dot_gases),
                                 &(urea_info->mass_parcel_reac_sp),
                                 &(urea_info->mass_dot_parcel_reac_sp),
                                 &(urea_info->c_isp),
                                 &(urea_info->k_cr),
                                 &(urea_info->mw),
                                 &(urea_info->a_react),
                                 &(urea_info->e_react)); */
}

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
#include "/home/tony/Downloads/instdir/include/cvode/cvode.h"            /* prototypes for CVODE fcts., consts.  */
#include "/home/tony/Downloads/instdir/include/cvode/cvode_diag.h"            /* prototypes for CVODE fcts., consts.  */
#include "/home/tony/Downloads/instdir/include/nvector/nvector_serial.h" /* access to serial N_Vector            */
#include "/home/tony/Downloads/instdir/include/sundials/sundials_types.h" /* access to dense SUNMatrix            */
#include "/home/tony/Downloads/instdir/include/sunlinsol/sunlinsol_dense.h" /* access to dense SUNLinearSolver      */
#include "/home/tony/Downloads/instdir/include/sunmatrix/sunmatrix_dense.h" /* access to dense SUNMatrix            */
#include "/home/tony/Downloads/instdir/include/nvector/nvector_serial.h" /* access to dense SUNNonlinearSolver */
#include "/home/tony/Downloads/instdir/include/sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to dense SUNMatrix            */
#include "/home/tony/Downloads/instdir/include/sunnonlinsol/sunnonlinsol_newton.h" /* access to dense SUNMatrix */
#include "user_header.h"
/**********************************************************************/
/*                                                                    */
/* Name: user_kh                                                      */
/*                                                                    */
/* Description: Kelvin-Helmholtz Breakup Model                        */
/*                                                                    */
/* Called when: user_kh_flag=1 in udf.in                              */
/*                                                                    */
/**********************************************************************/

static void init_tables(CONVERGE_species_t species);
static void destroy_tables(CONVERGE_species_t species);
static CONVERGE_table_t *pvap_table = NULL;

/* Helpful macros */
#ifndef SQR
#define SQR(A) ((A) * (A))
#endif
#ifndef CUB
#define CUB(A) ((A) * (A) * (A))
#endif

/* Shared Problem Constants */
#define ATOL SUN_RCONST(1.0e-10)
#define RTOL SUN_RCONST(1.0e-8)
#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define THREE  SUN_RCONST(3.0)
#define FOUR   SUN_RCONST(4.0)
#define EIGHT  SUN_RCONST(8.0)
#define PI     SUN_RCONST(3.141592653589793)

/* Problem Parameters */
CONVERGE_precision_t P1_T0, P1_T1, P1_DTOUT, P1_TOL_FACTOR;
CONVERGE_precision_t P1_RBInit, P1_RBDotInit, P1_CL, P1_Pv, P1_pAmbient, P1_sigma, P1_Pr, P1_Pr0, P1_mu, P1_kappa, P1_rho;

/* ODE Variables */
#define P1_NEQ  2
#define P1_NOUT 2

/* Function Prototypes */
static int f1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Problem1(sunrealtype* radius_array, CONVERGE_precision_t* bubble_growth_rate);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu, sunrealtype hu,
                        sunrealtype* radius_array, int index);
static void PrintFinalStats(void* cvode_mem, sunrealtype ero);

/* Check function return value */
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;
    if (opt == 0 && returnvalue == NULL)
    {
        fprintf(stderr, "SUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname);
        return 1;
    }
    if (opt == 1)
    {
        retval = (int*)returnvalue;
        if (*retval < 0)
        {
            fprintf(stderr, "SUNDIALS_ERROR: %s() failed with retval = %d\n", funcname, *retval);
            return 1;
        }
    }
    return 0;
}

/* Dimensional RHS function for simplified Rayleigh-Plesset equation */
static int f1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    sunrealtype R = NV_Ith_S(y, 0); // Dimensional radius (m)
    sunrealtype R_dot = NV_Ith_S(y, 1); // Dimensional velocity (m/s)

    /* Check for valid inputs */
    if (R <= 1e-15 || isnan(R) || P1_rho <= 0.0 || isnan(P1_rho))
    {
        fprintf(stderr, "Error in f1: t = %.12e, R = %.12e, R_dot = %.12e, P1_rho = %.12e, "
                        "P1_Pv = %.12e, P1_Pr0 = %.12e, P1_sigma = %.12e, P1_Pr = %.12e\n",
                t, R, R_dot, P1_rho, P1_Pv, P1_Pr0, P1_sigma, P1_Pr);
        return -1;
    }

    /* Simplified Rayleigh-Plesset equation (no viscosity terms) */
    sunrealtype term1 = P1_Pv / P1_rho; // Vapor pressure term
    sunrealtype term2 = (P1_Pr0 + 2.0 * P1_sigma / P1_RBInit) * pow(P1_RBInit / R, 3) / P1_rho; // Initial pressure term
    sunrealtype term3 = -2.0 * P1_sigma / (R * P1_rho); // Surface tension term
    sunrealtype term4 = -P1_Pr / P1_rho; // Ambient pressure term
    sunrealtype pressure_term = term1 + term2 + term3 + term4;
    sunrealtype accel = pressure_term - 1.5 * R_dot * R_dot / R;

    /* Log terms for debugging */
    if (t == P1_T0)
    {
        fprintf(stderr, "f1 terms at t = %.12e: R = %.12e, R_dot = %.12e\n", t, R, R_dot);
        fprintf(stderr, "  term1 = %.12e, term2 = %.12e, term3 = %.12e, term4 = %.12e\n",
                term1, term2, term3, term4);
        fprintf(stderr, "  pressure_term = %.12e, accel = %.12e\n", pressure_term, accel);
    }

    /* Dimensional derivatives */
    NV_Ith_S(ydot, 0) = R_dot; // dR/dt = R_dot
    NV_Ith_S(ydot, 1) = accel; // dR_dot/dt = acceleration

    if (t == P1_T0)
    {
        fprintf(stderr, "Initial derivatives: ydot[0] = %.12e, ydot[1] = %.12e\n",
                NV_Ith_S(ydot, 0), NV_Ith_S(ydot, 1));
    }

    return 0;
}

/* Dimensional Jacobian function for simplified Rayleigh-Plesset equation */
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sunrealtype R = NV_Ith_S(y, 0); // Dimensional radius (m)
    sunrealtype R_dot = NV_Ith_S(y, 1); // Dimensional velocity (m/s)

    /* Check for valid inputs */
    if (R <= 1e-15 || isnan(R) || P1_rho <= 0.0 || isnan(P1_rho))
    {
        fprintf(stderr, "Error in Jac: t = %.12e, R = %.12e, R_dot = %.12e, P1_rho = %.12e\n",
                t, R, R_dot, P1_rho);
        return -1;
    }

    /* Derivatives for simplified RP equation */
    sunrealtype A = P1_Pr0 + 2.0 * P1_sigma / P1_RBInit;
    sunrealtype B = -2.0 * P1_sigma;
    sunrealtype dp_dR = -3.0 * A * pow(P1_RBInit, 3) / pow(R, 4) / P1_rho - B / (R * R * P1_rho);
    sunrealtype J10 = (dp_dR * R - (P1_Pv/P1_rho + A * pow(P1_RBInit, 3) / pow(R, 3) / P1_rho +
                                    B / (R * P1_rho) - P1_Pr/P1_rho)) / (R * R) +
                      1.5 * R_dot * R_dot / (R * R);
    sunrealtype J11 = -3.0 * R_dot / R;

    /* Jacobian matrix */
    SM_ELEMENT_D(J, 0, 0) = 0.0;
    SM_ELEMENT_D(J, 0, 1) = 1.0;
    SM_ELEMENT_D(J, 1, 0) = J10;
    SM_ELEMENT_D(J, 1, 1) = J11;

    return 0;
}

/* Main SUNDIALS solver function */
int Problem1(sunrealtype* radius_array, CONVERGE_precision_t* bubble_growth_rate)
{
    void* cvode_mem = NULL;
    N_Vector y = NULL;
    SUNMatrix A = NULL;
    SUNLinearSolver LS = NULL;
    SUNContext ctx = NULL;
    sunrealtype t;
    int flag;

    /* Create SUNContext */
    flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
    if (check_retval(&flag, "SUNContext_Create", 1)) {
        return -1;
    }

    /* Create state vector */
    y = N_VNew_Serial(2, ctx); // 2D system (R, R_dot)
    if (check_retval((void*)y, "N_VNew_Serial", 0)) {
        SUNContext_Free(&ctx);
        return -1;
    }

    /* Initialize state vector */
    NV_Ith_S(y, 0) = radius_array[0]; // Initial radius (m)
    NV_Ith_S(y, 1) = radius_array[1]; // Initial growth rate (m/s)

    /* Create CVODE solver */
    cvode_mem = CVodeCreate(CV_BDF, ctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) {
        N_VDestroy(y);
        SUNContext_Free(&ctx);
        return -1;
    }

    /* Initialize CVODE */
    flag = CVodeInit(cvode_mem, f1, P1_T0, y);
    if (check_retval(&flag, "CVodeInit", 1)) {
        N_VDestroy(y);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    /* Set tolerances */
    flag = CVodeSStolerances(cvode_mem, P1_TOL_FACTOR, P1_TOL_FACTOR);
    if (check_retval(&flag, "CVodeSStolerances", 1)) {
        N_VDestroy(y);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    /* Set user data */
    flag = CVodeSetUserData(cvode_mem, NULL);
    if (check_retval(&flag, "CVodeSetUserData", 1)) {
        N_VDestroy(y);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    /* Set Jacobian */
    A = SUNDenseMatrix(2, 2, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) {
        N_VDestroy(y);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return -1;
    }

    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) {
        N_VDestroy(y);
        SUNMatDestroy(A);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return -1;
    }

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&flag, "CVodeSetLinearSolver", 1)) {
        N_VDestroy(y);
        SUNMatDestroy(A);
        SUNLinSolFree(LS);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    flag = CVodeSetJacFn(cvode_mem, Jac);
    if (check_retval(&flag, "CVodeSetJacFn", 1)) {
        N_VDestroy(y);
        SUNMatDestroy(A);
        SUNLinSolFree(LS);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    /* Integrate to P1_T1 */
    flag = CVode(cvode_mem, P1_T1, y, &t, CV_NORMAL);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error: CVode failed with flag = %d\n", flag);
        N_VDestroy(y);
        SUNMatDestroy(A);
        SUNLinSolFree(LS);
        CVodeFree(&cvode_mem);
        SUNContext_Free(&ctx);
        return flag;
    }

    /* Store results */
    radius_array[0] = NV_Ith_S(y, 0); // Final radius
    radius_array[1] = NV_Ith_S(y, 1); // Final growth rate
    *bubble_growth_rate = radius_array[1];

    /* Output integration stats */
    PrintFinalStats(cvode_mem, 0.0);

    /* Clean up */
    N_VDestroy(y);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
    CVodeFree(&cvode_mem);
    SUNContext_Free(&ctx);
    return 0; // Success
}

/* Output function */
static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu, sunrealtype hu,
                        sunrealtype* radius_array, int index)
{
    radius_array[index] = y0;
    /* fprintf(stderr, "t = %.12e, radius = %.12e, growth_rate = %.12e, order = %d, step = %.12e\n",
            t, y0, y1, qu, hu); */
}

/* Final statistics output */
static void PrintFinalStats(void* cvode_mem, sunrealtype ero)
{
    long int nst, nfe, ncfn, netf, nje, nfeLS;
    int flag;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    if (check_retval(&flag, "CVodeGetNumSteps", 1)) return;

    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    if (check_retval(&flag, "CVodeGetNumRhsEvals", 1)) return;

    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    if (check_retval(&flag, "CVodeGetNumErrTestFails", 1)) return;

    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    if (check_retval(&flag, "CVodeGetNumNonlinSolvConvFails", 1)) return;

    flag = CVodeGetNumJacEvals(cvode_mem, &nje);
    if (check_retval(&flag, "CVodeGetNumJacEvals", 1)) return;

    flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    if (check_retval(&flag, "CVodeGetNumLinRhsEvals", 1)) return;

    fprintf(stderr, "CVode stats: Steps = %ld, RHS evals = %ld, Error test fails = %ld, "
                    "Nonlinear conv fails = %ld, Jac evals = %ld, Lin RHS evals = %ld\n",
            nst, nfe, netf, ncfn, nje, nfeLS);
}

CONVERGE_UDF(spray_kh,
             IN(VALUE(CONVERGE_index_t, passed_parcel_idx),
                VALUE(CONVERGE_cloud_t, passed_spray_cloud),
                VALUE(CONVERGE_species_t, passed_species),
                VALUE(CONVERGE_precision_t, passed_density),
                VALUE(CONVERGE_precision_t, passed_pressure),
                VALUE(CONVERGE_precision_t *, passed_velocity),
                VALUE(CONVERGE_precision_t *, passed_species_mass_fraction)),
             OUT(CONVERGE_VOID))
{
   struct ParcelCloud parcel_cloud;
   load_user_cloud(&parcel_cloud, passed_spray_cloud);

   init_tables(passed_species);

   CONVERGE_index_t kh_new_parcel_flag;
   CONVERGE_index_t nnn, num_kh, ikh, temp_index;
   CONVERGE_precision_t drop_nu, rey_num, weber_term;
   CONVERGE_precision_t weber_gas, weber_liq, ohn, taylor;
   CONVERGE_precision_t wave_length_denom, frequency;
   CONVERGE_precision_t wave_length, growth_rate, breakup_time;
   CONVERGE_precision_t radius_eq1, radius_eq2, radius_old;
   CONVERGE_precision_t radius_equil, new_radius;
   CONVERGE_precision_t cnst2_parcel, balpha_parcel;
   CONVERGE_precision_t old_num_drop, new_parcel_num_drop;
   CONVERGE_precision_t new_parcel_mass, percent_lost;
   CONVERGE_precision_t newparcel_cutoff_parcel, num_drop_tmp;
   CONVERGE_precision_t shed_factor_parcel;
   CONVERGE_precision_t length_kh, length_cav, length_tur;
   CONVERGE_precision_t time_kh, time_cav, time_tur, time_collapse, time_thermal;
   CONVERGE_precision_t radius_cav, coef_area_reduction;
   CONVERGE_precision_t parcel_tke, parcel_eps, parcel_uprime, tke_term1;
   CONVERGE_precision_t time_burst, khact_c_tcav;
   CONVERGE_precision_t scale_ratio, scale_ratio_cav, scale_ratio_tur;
   CONVERGE_index_t kh_model_flag;
   CONVERGE_precision_t fr_prop, pvap, fr_temp, tot_pvap;
   CONVERGE_index_t it_prop, kh_no_enlarge_flag;
   CONVERGE_precision_t oldparent_radius;
   CONVERGE_precision_t liquid_mass_in_cell;
   CONVERGE_precision_t radius_bubble, Fp, rMean, Fs, Thermal_radius, new_radius_thermal, old_drop, old_radius;
   CONVERGE_precision_t deltat, ini_bubble_radius;

   CONVERGE_index_t spray_elsa_flag    = CONVERGE_get_int("lagrangian.elsa_flag");
   CONVERGE_index_t num_gas_species    = CONVERGE_species_num_gas(passed_species);
   CONVERGE_index_t num_fluid_species  = CONVERGE_species_num_fluid(passed_species);
   CONVERGE_index_t num_parcel_species = CONVERGE_species_num_parcel(passed_species);
   CONVERGE_precision_t dt             = CONVERGE_simulation_dt();
   CONVERGE_iterator_t psp_it;
   CONVERGE_species_parcel_iterator_create(passed_species, &psp_it);

   // set KH model constants cnst2_parcel and balpha_parcel based on the injector
   // that the parcel came from

   CONVERGE_injector_t injector = CONVERGE_get_injector_with_id(parcel_cloud.from_injector[passed_parcel_idx]);
   cnst2_parcel                 = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KH_CONST2);
   balpha_parcel                = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KH_BALPHA);
   khact_c_tcav                 = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KHACT_C_TCAV);
   coef_area_reduction          = parcel_cloud.area_reduction[passed_parcel_idx];

   ikh                = CONVERGE_injector_get_parameter_flag(injector, INJECTOR_KH_FLAG);
   kh_no_enlarge_flag = CONVERGE_injector_get_parameter_flag(injector, INJECTOR_KH_NO_ENLARGE_FLAG);

   // set the new parcel flag and the new parcel cutoff mass fraction

   kh_new_parcel_flag      = CONVERGE_injector_get_parameter_flag(injector, INJECTOR_KH_NEW_PARCEL_FLAG);
   newparcel_cutoff_parcel = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KH_NEW_PARCEL_CUTOFF);
   shed_factor_parcel      = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KH_SHED_FACTOR);

   // calculate liquid viscosity and liquid reynolds number

   drop_nu = parcel_cloud.viscosity[passed_parcel_idx] / parcel_cloud.density[passed_parcel_idx];
   rey_num = (parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.rel_vel_mag[passed_parcel_idx]) / drop_nu;

   // calculate gas and liquid phase weber numbers

   liquid_mass_in_cell = 0.0;
   if(spray_elsa_flag)
   {
      for(CONVERGE_index_t isp = num_gas_species; isp < num_fluid_species; isp++)
      {
         liquid_mass_in_cell += passed_species_mass_fraction[isp] * passed_density;
      }
   }

   weber_term = (parcel_cloud.rel_vel_mag[passed_parcel_idx] * parcel_cloud.rel_vel_mag[passed_parcel_idx] *
                 parcel_cloud.radius[passed_parcel_idx]) /
                parcel_cloud.surf_ten[passed_parcel_idx];
   weber_gas = (passed_density - liquid_mass_in_cell) * weber_term;
   weber_liq = parcel_cloud.density[passed_parcel_idx] * weber_term;

   // calculate ohnesorge and taylor numbers

   ohn    = CONVERGE_sqrt(weber_liq) / rey_num;
   taylor = ohn * CONVERGE_sqrt(weber_gas);

   // KH wavelength, Eq. 4 in KH model paper referenced in header

   wave_length_denom = CONVERGE_pow((1.0 + 0.87 * CONVERGE_pow(weber_gas, 1.67)), 0.6);
   wave_length       = parcel_cloud.radius[passed_parcel_idx] * (9.02 * (1.0 + 0.45 * CONVERGE_sqrt(ohn)) *
                                                           (1.0 + 0.4 * CONVERGE_pow(taylor, 0.7)) / wave_length_denom);

   // KH growth rate, Eq. 5 in KH model paper referenced in header

   growth_rate =
      (0.34 + 0.38 * weber_gas * CONVERGE_sqrt(weber_gas)) / ((1.0 + ohn) * (1.0 + 1.4 * CONVERGE_pow(taylor, 0.6)));
   frequency   = CONVERGE_sqrt(parcel_cloud.surf_ten[passed_parcel_idx] /
                             (parcel_cloud.density[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                              parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx]));
   growth_rate = growth_rate * frequency;

   // KH breakup time, Eq. 12 in KH model paper referenced in header

   radius_equil = balpha_parcel * wave_length;
   breakup_time = 3.726 * cnst2_parcel * parcel_cloud.radius[passed_parcel_idx] / (growth_rate * wave_length);
   new_radius =
      (parcel_cloud.radius[passed_parcel_idx] + (dt / breakup_time) * radius_equil) / (1.0 + (dt / breakup_time));

   kh_model_flag = 1;

   length_kh   = parcel_cloud.radius[passed_parcel_idx] - radius_equil;
   time_kh     = breakup_time;
   scale_ratio = length_kh / time_kh;
   // KH-ACT Model
   if(ikh == 2 && parcel_cloud.parent[passed_parcel_idx] == 1)
   {
      // clip the tke0 to make sure it is greater than 0
      parcel_cloud.tke0[passed_parcel_idx] = fmax(parcel_cloud.tke0[passed_parcel_idx], 1.0e-10);

      // get the parcel properties given the current temperature and current mass fractions
      temp_index = (CONVERGE_index_t)(parcel_cloud.temp[passed_parcel_idx] / 10.0);
      fr_temp    = parcel_cloud.temp[passed_parcel_idx] / 10.0 - (CONVERGE_precision_t)temp_index;
      tot_pvap   = 0.0;
      for(CONVERGE_index_t isp = 0, isp_global = CONVERGE_iterator_first(psp_it); isp < num_parcel_species;
          isp++, isp_global                    = CONVERGE_iterator_next(psp_it))
      {
         const CONVERGE_precision_t isp_tcrit = CONVERGE_species_tcrit(passed_species, isp_global);
         if(parcel_cloud.temp[passed_parcel_idx] < isp_tcrit)
         {
            it_prop = temp_index;
            fr_prop = fr_temp;
         }
         else
         {
            it_prop = (CONVERGE_index_t)(isp_tcrit / 10.0);
            fr_prop = isp_tcrit / 10.0 - (CONVERGE_precision_t)it_prop;
         }
         double temp1 = (it_prop + fr_prop) * 10.;
         pvap         = CONVERGE_table_lookup(pvap_table[isp], temp1);

         if(pvap < 1.0e-10)
         {
            pvap = 1.0e-10;
         }
         if(pvap > passed_pressure)
         {
            pvap = passed_pressure;
         }
         tot_pvap = tot_pvap + parcel_cloud.mfrac[num_parcel_species * passed_parcel_idx + isp] * pvap;
      }

      CONVERGE_precision_t turb_in_ceps2 = CONVERGE_get_double("turbulence.ceps2");
      tke_term1                          = parcel_cloud.tke0[passed_parcel_idx] + parcel_cloud.eps0[passed_parcel_idx] *
                                                            parcel_cloud.lifetime[passed_parcel_idx] *
                                                            (turb_in_ceps2 - 1.0);
      tke_term1 = (CONVERGE_pow(parcel_cloud.tke0[passed_parcel_idx], turb_in_ceps2) / tke_term1);

      parcel_tke = CONVERGE_pow(tke_term1, (1.0 / (turb_in_ceps2 - 1.0)));
      parcel_eps = parcel_cloud.eps0[passed_parcel_idx] *
                   CONVERGE_pow(parcel_tke / parcel_cloud.tke0[passed_parcel_idx], turb_in_ceps2);
      length_tur = CONVERGE_get_double("turbulence.keps_cmu") * CONVERGE_pow_3over2(parcel_tke) / parcel_eps;
      time_tur   = CONVERGE_get_double("turbulence.keps_cmu") * (parcel_tke / parcel_eps);

      CONVERGE_nozzle_t nozzle =
         CONVERGE_injector_get_nozzle_with_id(injector, parcel_cloud.from_nozzle[passed_parcel_idx]);
      CONVERGE_precision_t noz_diameter = CONVERGE_nozzle_get_diameter(nozzle);
      radius_cav                        = noz_diameter / 2.0 * CONVERGE_sqrt(1.0 - coef_area_reduction);
      time_collapse = 0.9145 * radius_cav * CONVERGE_sqrt(parcel_cloud.density[passed_parcel_idx] / tot_pvap);
      if(radius_cav == 0.0)
      {
         time_collapse = 1.0e5;
      }

      parcel_uprime = CONVERGE_sqrt(2.0 / 3.0 * parcel_tke);
      time_burst    = (noz_diameter / 2.0 - radius_cav) / parcel_uprime;
      length_cav    = radius_cav;
      time_cav      = fmin(time_burst, time_collapse);

      scale_ratio     = length_kh / time_kh;
      scale_ratio_cav = length_cav / time_cav;
      scale_ratio_tur = length_tur / time_tur;

      if(scale_ratio < scale_ratio_cav)
      {
         scale_ratio   = scale_ratio_cav;
         kh_model_flag = 0;
      }
      if(scale_ratio < scale_ratio_tur)
      {
         scale_ratio   = scale_ratio_tur;
         kh_model_flag = 0;
      }
   }


   if(kh_no_enlarge_flag == 1)   // user can use this flag to turn off the drop enlargement code below
   {
      parcel_cloud.tbreak_kh[passed_parcel_idx] = 0.0;
   }

   if(wave_length > (parcel_cloud.radius[passed_parcel_idx] / balpha_parcel))
   {
      if(parcel_cloud.tbreak_kh[passed_parcel_idx] > 0.0)
      {
         // KH model drop enlargement

         if(parcel_cloud.tbreak_kh[passed_parcel_idx] < breakup_time)
         {
            parcel_cloud.tbreak_kh[passed_parcel_idx] = parcel_cloud.tbreak_kh[passed_parcel_idx] + dt;
         }
         else
         {
            // breakup occurs

            parcel_cloud.tbreak_kh[passed_parcel_idx] = 0.0;

            // calculate equilibrium radius, Eq. 10b in KH model paper referenced in header

            radius_eq1 = CONVERGE_cbrt(0.75 * wave_length * parcel_cloud.radius[passed_parcel_idx] *
                                       parcel_cloud.radius[passed_parcel_idx]);
            radius_eq2 = CONVERGE_cbrt(1.5 * PI * parcel_cloud.radius[passed_parcel_idx] *
                                       parcel_cloud.radius[passed_parcel_idx] *
                                       parcel_cloud.rel_vel_mag[passed_parcel_idx] / growth_rate);
            radius_old = parcel_cloud.radius[passed_parcel_idx];

            // update parcel's radius (minimum of radius_eq1 and radius_eq2) and drop number
            // also zero the RT breakup time and shed_num_drop variable if breakup takes place

            parcel_cloud.radius[passed_parcel_idx]     = (radius_eq1 < radius_eq2) ? radius_eq1 : radius_eq2;
            parcel_cloud.radius_tm1[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
            parcel_cloud.num_drop[passed_parcel_idx] =
               parcel_cloud.num_drop[passed_parcel_idx] * radius_old * radius_old * radius_old /
               (parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                parcel_cloud.radius[passed_parcel_idx]);
            parcel_cloud.tbreak_rt[passed_parcel_idx]     = 0.0;
            parcel_cloud.shed_num_drop[passed_parcel_idx] = 0.0;
            parcel_cloud.shed_mass[passed_parcel_idx]     = 0.0;
         }
      }
   }
   // KH model drop breakup
   else
   {
      // calculate equilibrium/child drop radius, Eq. 10a in KH model paper referenced in header

      radius_equil = balpha_parcel * wave_length;

      // calculate the updated radius, based on Eq. 11 in KH model paper referenced in header

      if(kh_model_flag == 1)
      {
          // Initialize deltaT
          deltat = parcel_cloud.temp[passed_parcel_idx] - 349;
          CONVERGE_precision_t MIN_DELTAT = 1.0e-6;
          if (deltat <= 0) {
              deltat = MIN_DELTAT;
          }

          // Bubble initialization
          const CONVERGE_precision_t C1 = 1.0e-23;
          const CONVERGE_precision_t C2 = -0.3;
          ini_bubble_radius = fmax(CONVERGE_cbrt(C1 * exp(C2 / deltat)), 1.0e-7);
          P1_RBInit = fmin(ini_bubble_radius, parcel_cloud.radius[passed_parcel_idx] / 10.0);
          P1_RBDotInit = 1.0e-3;
          P1_T0 = 0.0;
          P1_T1 = 1.5e-2; // 15 ms
          P1_DTOUT = dt / 12;
          P1_TOL_FACTOR = 1.0e-6;
          P1_Pv = fuel_pvap;
          P1_pAmbient = ambient_pres;
          P1_sigma = parcel_cloud.surf_ten[passed_parcel_idx];
          P1_Pr = P1_pAmbient;
          P1_Pr0 = P1_pAmbient;
          P1_mu = parcel_cloud.viscosity[passed_parcel_idx]; // Not used
          P1_kappa = 0.0; // Not used
          P1_rho = parcel_cloud.density[passed_parcel_idx];
          sunrealtype y0_array[P1_NOUT];

          // Initialize ODE
          y0_array[0] = P1_RBInit;
          y0_array[1] = P1_RBDotInit;

          // Save old values
          old_radius = parcel_cloud.radius[passed_parcel_idx];
          old_drop = parcel_cloud.num_drop[passed_parcel_idx];

          CONVERGE_precision_t bubble_growth_rate;
          int ode_status = Problem1(y0_array, &bubble_growth_rate);
          if (ode_status != 0) {
              printf("Warning: ODE solver failed for parcel %ld, using initial bubble radius\n", passed_parcel_idx);
              radius_bubble = P1_RBInit;
              bubble_growth_rate = P1_RBDotInit;
          } else {
              radius_bubble = y0_array[0];
              bubble_growth_rate = y0_array[1];
          }
          printf("radius_bubble = %.12e, bubble_growth_rate = %.12e\n, ini_bubble_radius = %.12e\n, old_radius = %.12e\n", radius_bubble, bubble_growth_rate, ini_bubble_radius, old_radius);
          radius_bubble = fmin(radius_bubble, parcel_cloud.radius[passed_parcel_idx] * 0.99);
          if (radius_bubble <= 0.0) {
              printf("Warning: radius_bubble clamped to 1.0e-8 for parcel %ld\n", passed_parcel_idx);
              radius_bubble = 1.0e-8;
          }

          // Calculate forces and proceed with thermal breakup...

          Fp = 4 * PI * (P1_Pv * CONVERGE_pow(radius_bubble, 2) - P1_pAmbient * CONVERGE_pow(old_radius, 2)); // Pressure force
          printf("Fp = %f\n", Fp);
          rMean = (radius_bubble + old_radius) / 2; // Mean radius
          printf("rMean = %f\n", rMean);
          Fs = 2 * PI * rMean * P1_sigma; // Surface tension force
          // printf("Fs = %f\n", Fs);

         if (Fp > Fs) {
             CONVERGE_LAGRANGE_Spray_Parcel_make_small(passed_parcel_idx, passed_spray_cloud);

             // Save old values for mass conservation check
             CONVERGE_precision_t old_volume = CONVERGE_pow(old_radius, 3);
             CONVERGE_precision_t bubble_volume = CONVERGE_pow(radius_bubble, 3);

             // Calculate Thermal_radius
             Thermal_radius = CONVERGE_cbrt((old_volume - bubble_volume) / 2);
             if (Thermal_radius <= 0.0) {
                 printf("Warning: Thermal_radius clamped to 1.0e-6 for parcel %ld\n", passed_parcel_idx);
                 Thermal_radius = 1.0e-6;
             }
             printf("Thermal radius = %.12e\n", Thermal_radius);

             // Calculate old and new mass for verification
             CONVERGE_precision_t mass_old = parcel_cloud.num_drop[passed_parcel_idx] * (4.0 / 3.0) * PI *
                                             old_volume * parcel_cloud.density[passed_parcel_idx];

             // Update parent parcel's radius and drop number
             parcel_cloud.radius[passed_parcel_idx] = Thermal_radius;
             parcel_cloud.radius_tm1[passed_parcel_idx] = Thermal_radius;
             parcel_cloud.num_drop[passed_parcel_idx] =
                 parcel_cloud.num_drop[passed_parcel_idx] * old_volume / CONVERGE_pow(Thermal_radius, 3);

             // Verify mass conservation
             CONVERGE_precision_t mass_new = parcel_cloud.num_drop[passed_parcel_idx] * (4.0 / 3.0) * PI *
                                            CONVERGE_pow(Thermal_radius, 3) * parcel_cloud.density[passed_parcel_idx];
             if (fabs(mass_new - mass_old) / mass_old > 1.0e-6) {
                 printf("Warning: Mass conservation error for parcel %ld: old_mass = %.12e, new_mass = %.12e\n",
                        passed_parcel_idx, mass_old, mass_new);
             }
             printf("old_mass = %.12e, new_mass = %.12e\n",
                         mass_old, mass_new);

             /* since breakup has taken place, zero some of the parcel properties (resetting some of these properties
             can potentially have significant effects on the simulation) */

             // parcel_cloud.tbreak_rt[passed_parcel_idx] = 0.0;
             // parcel_cloud.distort[passed_parcel_idx] = 0.0;
             // parcel_cloud.distort_dot[passed_parcel_idx] = 0.0;
             parcel_cloud.shed_num_drop[passed_parcel_idx] = 0.0;
             parcel_cloud.shed_mass[passed_parcel_idx] = 0.0;
             // parcel_cloud.tbreak_kh[passed_parcel_idx] = 0.0;
             // parcel_cloud.tke0[passed_parcel_idx] = 0.0;
             // parcel_cloud.eps0[passed_parcel_idx] = 0.0;
             // parcel_cloud.lifetime[passed_parcel_idx] = 0.0;
             // parcel_cloud.area_reduction[passed_parcel_idx] = 0.0;
          }
          // kh model breakup happens here
          new_radius =
                  (parcel_cloud.radius[passed_parcel_idx] + (dt / breakup_time) * radius_equil) / (1.0 + (dt / breakup_time));
      }
      else
      // kh-act model activation which occurs if the kh model flag is 0
      {
         new_radius = parcel_cloud.radius[passed_parcel_idx] - khact_c_tcav * scale_ratio * dt;
         new_radius = fmax(new_radius, 1.0e-20);
      }

      // save parcel's drop number and calculate updated drop number
      old_num_drop = parcel_cloud.num_drop[passed_parcel_idx];
      parcel_cloud.num_drop[passed_parcel_idx] =
         old_num_drop * parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
         parcel_cloud.radius[passed_parcel_idx] / (new_radius * new_radius * new_radius);

      // if new (child) parcels are not allowed, update the parcel radius

      if(kh_new_parcel_flag == 0)
      {
         parcel_cloud.radius[passed_parcel_idx]     = new_radius;
         parcel_cloud.radius_tm1[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
      }
      // if new parcels are allowed, calculate how many child parcels are to be added
      else
      {
         if(parcel_cloud.shed_num_drop[passed_parcel_idx] == 0.0)
         {
            // initialize shed_num_drop
            parcel_cloud.shed_num_drop[passed_parcel_idx] = old_num_drop;
         }

         parcel_cloud.shed_mass[passed_parcel_idx] =
            parcel_cloud.shed_mass[passed_parcel_idx] +
            shed_factor_parcel * parcel_cloud.shed_num_drop[passed_parcel_idx] * PI * (4.0 / 3.0) *
               parcel_cloud.density[passed_parcel_idx] *
               (parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                   parcel_cloud.radius[passed_parcel_idx] -
                new_radius * new_radius * new_radius);

         CONVERGE_precision_t mass_per_parcel =
         CONVERGE_injector_get_parameter_precision(injector, INJECTOR_MASS_PER_PARCEL);
         percent_lost = parcel_cloud.shed_mass[passed_parcel_idx] / (mass_per_parcel * newparcel_cutoff_parcel);

         num_kh = (CONVERGE_index_t)(percent_lost);

         if(num_kh == 0)
         {
            // do not add child parcels yet, but update parent parcel radius (num_drop was updated above)
            parcel_cloud.radius[passed_parcel_idx]     = new_radius;
            parcel_cloud.radius_tm1[passed_parcel_idx] = parcel_cloud.radius[passed_parcel_idx];
         }
         else
         {
            // one or more parcels is predicted to shed off
            CONVERGE_LAGRANGE_Spray_Parcel_make_small(passed_parcel_idx, passed_spray_cloud);
            new_parcel_mass     = newparcel_cutoff_parcel * mass_per_parcel;
            new_parcel_num_drop = new_parcel_mass / (PI * (4.0 / 3.0) * parcel_cloud.density[passed_parcel_idx] *
                                                     radius_equil * radius_equil * radius_equil);

            for(nnn = 0; nnn < num_kh; nnn++)
            {
               // calculate the parent parcel's drop number

               num_drop_tmp = (old_num_drop * parcel_cloud.radius[passed_parcel_idx] *
                                  parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] -
                               ((CONVERGE_precision_t)(num_kh)) * new_parcel_num_drop * radius_equil * radius_equil *
                                  radius_equil) /
                              (new_radius * new_radius * new_radius);
               parcel_cloud.num_drop[passed_parcel_idx] = old_num_drop;

               // create new child parcel and initialize its properties

               CONVERGE_precision_t kh_const1 = CONVERGE_injector_get_parameter_precision(injector, INJECTOR_KH_CONST1);
               CONVERGE_spray_child_parcel(passed_velocity,
                                           growth_rate,
                                           wave_length,
                                           radius_equil,
                                           new_parcel_num_drop,
                                           kh_const1,
                                           passed_parcel_idx,
                                           passed_spray_cloud);
               // reload after adding parcels
               load_user_cloud(&parcel_cloud, passed_spray_cloud);

               // since breakup occurred, zero the parent parcel's RT breakup time
               parcel_cloud.num_drop[passed_parcel_idx] = fmax(num_drop_tmp, 0.0);
               parcel_cloud.tbreak_rt[passed_parcel_idx] = 0.0;
            }

            // update parent drop's radius and shed_mass

            oldparent_radius                              = parcel_cloud.radius[passed_parcel_idx];
            parcel_cloud.radius[passed_parcel_idx]        = new_radius;
            parcel_cloud.radius_tm1[passed_parcel_idx]    = parcel_cloud.radius[passed_parcel_idx];
            parcel_cloud.shed_num_drop[passed_parcel_idx] = 0.0;
            parcel_cloud.shed_mass[passed_parcel_idx] =
               parcel_cloud.shed_mass[passed_parcel_idx] - (CONVERGE_precision_t)(num_kh)*new_parcel_mass;

            //parent parcel's density has not changed, so updating sactive based on volume scaling
            parcel_cloud.sactive[passed_parcel_idx] =
               (parcel_cloud.num_drop[passed_parcel_idx]/old_num_drop) * CONVERGE_cube(new_radius / oldparent_radius) * parcel_cloud.sactive[passed_parcel_idx];
            parcel_cloud.sactive_tm1[passed_parcel_idx] = parcel_cloud.sactive[passed_parcel_idx];
         }
      }
   }

   CONVERGE_LAGRANGE_Parcel_discretize_big(CONVERGE_get_int("lagrangian.evap_layers_per_drop"),
                                           0,
                                           passed_parcel_idx,
                                           passed_spray_cloud);
   destroy_tables(passed_species);
   CONVERGE_iterator_destroy(&psp_it);
}

void init_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   // Get the parcel species tables
   load_species_tables(parcel_species_it, TEMPERATURE_TABLE_PVAP_ID, &pvap_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}

static void destroy_tables(CONVERGE_species_t species)
{
   CONVERGE_iterator_t parcel_species_it;
   CONVERGE_species_parcel_iterator_create(species, &parcel_species_it);

   // Destroy local parcel tables
   unload_species_tables(parcel_species_it, &pvap_table);

   CONVERGE_iterator_destroy(&parcel_species_it);
}

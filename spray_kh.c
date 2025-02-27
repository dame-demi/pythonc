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
#include <stdbool.h>
#include "/home/tony/Downloads/instdir/include/cvode/cvode.h"            /* prototypes for CVODE fcts., consts.  */
#include "/home/tony/Downloads/instdir/include/cvode/cvode_diag.h"            /* prototypes for CVODE fcts., consts.  */
#include "/home/tony/Downloads/instdir/include/nvector/nvector_serial.h" /* access to serial N_Vector            */
#include "/home/tony/Downloads/instdir/include/sundials/sundials_types.h" /* access to dense SUNMatrix            */
#include "/home/tony/Downloads/instdir/include/sunlinsol/sunlinsol_dense.h" /* access to dense SUNLinearSolver      */
#include "/home/tony/Downloads/instdir/include/sunmatrix/sunmatrix_dense.h" /* access to dense SUNMatrix            */
#include "/home/tony/Downloads/instdir/include/sunlinsol/sunlinsol_band.h" /* access to banded linear solver            */
#include "/home/tony/Downloads/instdir/include/sunmatrix/sunmatrix_band.h" /* access to banded SUNMatrix            */
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

/* helpful macros */

#ifndef SQR
#define SQR(A) ((A) * (A))
#endif

#ifndef CUB
#define CUB(A) ((A) * (A) * (A))
#endif

/* Shared Problem Constants */

#define ATOL SUN_RCONST(1.0e-6)
#define RTOL SUN_RCONST(1.0e-4)

#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define THREE  SUN_RCONST(3.0)
#define FOUR   SUN_RCONST(4.0)
#define EIGHT  SUN_RCONST(8.0)

// Declare the variables for the Sundials ODE
CONVERGE_precision_t P1_T0, P1_T1, P1_DTOUT, P1_TOL_FACTOR;
CONVERGE_precision_t P1_RBInit, P1_RBDotInit, P1_CL, P1_Pv, P1_pAmbient, P1_sigma, P1_Pr, P1_Pr0, P1_mu, P1_kappa, P1_rho;

// ODE Variables

#define P1_NEQ        2
#define P1_NOUT       10

/* Linear Solver Options */
/* Enum for solver options */
enum
{
   FUNC,
   DENSE_DQ,
   DIAG,
   BAND_USER,
   BAND_DQ
 };

/* Private Helper Functions */

static void Problem1(sunrealtype* y1_array);
static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu, sunrealtype hu, sunrealtype* y1_array, int index);
static int PrepareNextRun(SUNContext sunctx, void* cvode_mem, int lmm, int miter, N_Vector y, SUNMatrix* A, sunindextype mu, sunindextype ml, SUNLinearSolver* LS, SUNNonlinearSolver* NLS);
static void PrintErrOutput(sunrealtype tol_factor);
static void PrintFinalStats(void* cvode_mem, int miter, sunrealtype ero);
static void PrintErrInfo(int nerr);

/* Functions Called by the Solver */

static int f1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */

static int check_retval(void* returnvalue, const char* funcname, int opt);

static void Problem1(sunrealtype* y1_array)
{
    sunrealtype reltol = RTOL, abstol = ATOL, t, tout, ero, er;
    int miter, retval, temp_retval, iout, nerr = 0;
    N_Vector y;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    void* cvode_mem;
    sunbooleantype firstrun;
    int qu;
    sunrealtype hu;
    SUNContext sunctx;

    y         = NULL;
    A         = NULL;
    LS        = NULL;
    NLS       = NULL;
    cvode_mem = NULL;

    /* Create the SUNDIALS context */
    retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) { return; }

    y = N_VNew_Serial(P1_NEQ, sunctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) { return; }
    /*PrintIntro1();*/

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    NLS = NULL;
    LS  = NULL;
    A   = NULL;

    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return; }

    for (miter = DENSE_DQ; miter <= DENSE_DQ; miter++)
    {
        ero            = ZERO;
        NV_Ith_S(y, 0) = P1_RBInit;
        NV_Ith_S(y, 1) = P1_RBDotInit;

        firstrun = (miter == DENSE_DQ);

        retval = CVodeInit(cvode_mem, f1, P1_T0, y);
        if (check_retval(&retval, "CVodeInit", 1)) { return; }

        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        if (check_retval(&retval, "CVodeSStolerances", 1)) { return; }

        retval = CVodeReInit(cvode_mem, P1_T0, y);
        if (check_retval(&retval, "CVodeReInit", 1)) { return; }

        retval = PrepareNextRun(sunctx, cvode_mem, CV_BDF, miter, y, &A, 0, 0, &LS, &NLS);
        if (check_retval(&retval, "PrepareNextRun", 1)) { return; }

        /*PrintHeader1();*/

        for (iout = 1, tout = P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT)
        {
            retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
            check_retval(&retval, "CVode", 1);
            temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
            if (check_retval(&temp_retval, "CVodeGetLastOrder", 1)) { ++nerr; }
            temp_retval = CVodeGetLastStep(cvode_mem, &hu);
            if (check_retval(&temp_retval, "CVodeGetLastStep", 1)) { ++nerr; }
            PrintOutput1(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), qu, hu, y1_array, iout - 1);
            if (retval != CV_SUCCESS)
            {
                nerr++;
                break;
            }
            if (iout % 2 == 0)
            {
                er = fabs(NV_Ith_S(y, 0)) / abstol;
                if (er > ero) { ero = er; }
                if (er > P1_TOL_FACTOR)
                {
                    nerr++;
                    PrintErrOutput(P1_TOL_FACTOR);
                }
            }
        }

        /*PrintFinalStats(cvode_mem, miter, ero);*/
    }

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    N_VDestroy(y);
    SUNContext_Free(&sunctx);
}

static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu, sunrealtype hu, sunrealtype* y1_array, int index)
{
   // Save y1 into the array at the specified index
   y1_array[index] = y0;

#if defined(SUNDIALS_EXTENDED_PRECISION)
   printf("%10.5Lf    %12.5Le   %12.5Le   %2d    %6.4Le\n", t, y0, y1, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
   printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#else
   printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#endif

   return;
}

static int f1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
   sunrealtype y0, y1;

   y0 = NV_Ith_S(y, 0);
   y1 = NV_Ith_S(y, 1);

   NV_Ith_S(ydot, 0) = y1;
   NV_Ith_S(ydot, 1) = ((ONE + y1 / P1_CL) * (
                          P1_Pv + (P1_Pr0 + TWO * P1_sigma / P1_RBInit) * CUB(P1_RBInit/y0) - TWO * P1_sigma / y0 - FOUR
                          * P1_mu * y1 / y0 - FOUR * P1_kappa * y1 / SQR(y0) - P1_Pr) / P1_rho + y0 * (
                          ONE / (P1_CL * P1_rho)) * (
                          -THREE * CUB(P1_RBInit/y0) * (P1_Pr0 + TWO * P1_sigma / P1_RBInit) * (y1 / y0) + TWO * P1_sigma
                          * y1 / SQR(y0) + FOUR * P1_mu * SQR(y1) / SQR(y0) + EIGHT * P1_kappa * SQR(y1) / CUB(y0)) - (
                          ONE - y1 / (THREE * P1_CL)) * (THREE * SQR(y1) / TWO)) / (
                         (ONE - y1 / P1_CL) * y0 + y0 * (ONE / (P1_CL * P1_rho)) * (
                           FOUR * P1_mu / y0 + FOUR * P1_kappa / SQR(y0)));

   return (0);
}

static int PrepareNextRun(SUNContext sunctx, void* cvode_mem, int lmm,
                          int miter, N_Vector y, SUNMatrix* A, sunindextype mu,
                          sunindextype ml, SUNLinearSolver* LS,
                          SUNNonlinearSolver* NLS)
{
  int retval = CV_SUCCESS;

  if (*NLS) { SUNNonlinSolFree(*NLS); }
  if (*LS) { SUNLinSolFree(*LS); }
  if (*A) { SUNMatDestroy(*A); }

  printf("\n\n-------------------------------------------------------------");

  printf("\n\nLinear Multistep Method : ");
  if (lmm == CV_ADAMS) { printf("ADAMS\n"); }
  else { printf("BDF\n"); }

  printf("Iteration               : ");
  if (miter == FUNC)
  {
    printf("FIXEDPOINT\n");

    /* create fixed point nonlinear solver object */
    *NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
    if (check_retval((void*)*NLS, "SUNNonlinSol_FixedPoint", 0)) { return (1); }

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if (check_retval(&retval, "CVodeSetNonlinearSolver", 1)) { return (1); }
  }
  else
  {
    printf("NEWTON\n");

    /* create Newton nonlinear solver object */
    *NLS = SUNNonlinSol_Newton(y, sunctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) { return (1); }

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if (check_retval(&retval, "CVodeSetNonlinearSolver", 1)) { return (1); }

    printf("Linear Solver           : ");

    switch (miter)
    {
    case DENSE_DQ:
      printf("Dense, Difference Quotient Jacobian\n");

      /* Create dense SUNMatrix for use in linear solves */
      *A = SUNDenseMatrix(P1_NEQ, P1_NEQ, sunctx);
      if (check_retval((void*)*A, "SUNDenseMatrix", 0)) { return (1); }

      /* Create dense SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Dense(y, *A, sunctx);
      if (check_retval((void*)*LS, "SUNLinSol_Dense", 0)) { return (1); }

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

      /* Use a difference quotient Jacobian */
      retval = CVodeSetJacFn(cvode_mem, NULL);
      if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }
      break;

    case DIAG:
      printf("Diagonal Jacobian\n");

      /* Call CVDiag to create/attach the CVODE-specific diagonal solver */
      retval = CVDiag(cvode_mem);
      if (check_retval(&retval, "CVDiag", 1)) { return (1); }
      break;

    }
  }

  return (retval);
}

static void PrintErrOutput(sunrealtype tol_factor)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
   printf("\n\n Error exceeds %Lg * tolerance \n\n", tol_factor);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
   printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#else
   printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#endif

   return;
}

static void PrintFinalStats(void* cvode_mem, int miter, sunrealtype ero)
{
  long int lenrw, leniw, lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf, nje, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  printf("\n Final statistics for this run:\n\n");
  printf(" CVode real workspace length              = %4ld \n", lenrw);
  printf(" CVode integer workspace length           = %4ld \n", leniw);
  printf(" Number of steps                          = %4ld \n", nst);
  printf(" Number of f-s                            = %4ld \n", nfe);
  printf(" Number of setups                         = %4ld \n", nsetups);
  printf(" Number of nonlinear iterations           = %4ld \n", nni);
  printf(" Number of nonlinear convergence failures = %4ld \n", ncfn);
  printf(" Number of error test failures            = %4ld \n\n", netf);

  if (miter != FUNC)
  {
    if (miter != DIAG)
    {
      retval = CVodeGetNumJacEvals(cvode_mem, &nje);
      check_retval(&retval, "CVodeGetNumJacEvals", 1);
      retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);
      retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVodeGetLinWorkSpace", 1);
    }
    else
    {
      nje    = nsetups;
      retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVDiagGetNumRhsEvals", 1);
      retval = CVDiagGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVDiagGetWorkSpace", 1);
    }
    printf(" Linear solver real workspace length      = %4ld \n", lenrwLS);
    printf(" Linear solver integer workspace length   = %4ld \n", leniwLS);
    printf(" Number of Jacobian evaluations           = %4ld \n", nje);
    printf(" Number of f evals. in linear solver      = %4ld \n\n", nfeLS);
  }

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" Error overrun = %.3Lf \n", ero);
#else
  printf(" Error overrun = %.3f \n", ero);
#endif
}

static void PrintErrInfo(int nerr)
{
   printf("\n\n-------------------------------------------------------------");
   printf("\n-------------------------------------------------------------");
   printf("\n\n Number of errors encountered = %d \n", nerr);

   return;
}

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
   int* retval;

   /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
   if (opt == 0 && returnvalue == NULL)
   {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return (1);
   }

   /* Check if retval < 0 */
   else if (opt == 1)
   {
      retval = (int*)returnvalue;
      if (*retval < 0)
      {
         fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                 funcname, *retval);
         return (1);
      }
   }

   /* Check if function returned NULL pointer - no memory allocated */
   else if (opt == 2 && returnvalue == NULL)
   {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return (1);
   }

   return (0);
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
   printf ("This is the start of spray kh.c \n");
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
   CONVERGE_precision_t time_kh, time_cav, time_tur, time_collapse;
   CONVERGE_precision_t radius_cav, coef_area_reduction;
   CONVERGE_precision_t parcel_tke, parcel_eps, parcel_uprime, tke_term1;
   CONVERGE_precision_t time_burst, khact_c_tcav;
   CONVERGE_precision_t scale_ratio, scale_ratio_cav, scale_ratio_tur;
   CONVERGE_index_t kh_model_flag;
   CONVERGE_precision_t fr_prop, pvap, fr_temp, tot_pvap;
   CONVERGE_index_t it_prop, kh_no_enlarge_flag;
   CONVERGE_precision_t oldparent_radius;
   CONVERGE_precision_t liquid_mass_in_cell;
   CONVERGE_vec3_t original_y0_array;
   CONVERGE_precision_t radius_bubble, rMean; // bubble radius and mean radius after thermal breakup
   CONVERGE_precision_t therm_iter;
   CONVERGE_precision_t reduced_dt; // reduced time step for breakup model
   CONVERGE_precision_t Fp, Fs; // pressure force
   CONVERGE_precision_t new_radius_thermal; // surface tension force
   bool breakup; // breakup flag

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

      if(kh_model_flag == 1) {
         therm_iter = 10; // number of iterations for thermal breakup
         reduced_dt = dt / therm_iter; // reduced time step for thermal breakup

         printf ("spray_kh.c sigma = %g vap_den = %g fuel_visc = %g\n", sigma, fuel_den, fuel_visc);
         printf ("Non-zero values means that the inputs for the bubble ODE are changing \n");

         // using Converge variables for the thermal breakup implementation
         P1_T0 = 0.0; // initial start time for ODE
         P1_T1 = 1.5e-7; // initial output time from ODE
         P1_DTOUT = dt; // Using reduced time step for ODE iteration
         P1_TOL_FACTOR = 1.0e10; // Tolerance factor to check how much the error exceeds the tolerance factor
         P1_RBInit = 0.0215e1; // Initial radius of the bubble
         P1_RBDotInit = 0.1; // Initial radius growth rate of the bubble
         P1_CL = sound_speed; // Speed of sound in the liquid
         P1_Pv = fuel_pvap; // Vapor pressure of the liquid
         P1_pAmbient = ambient_pres; // Ambient pressure in the system
         P1_sigma = sigma; // Surface tension of the liquid
         P1_Pr = 0.71; // Prandtl number
         P1_Pr0 = 0.71; // Initial Prandtl number
         P1_mu = fuel_visc; // Dynamic viscosity of the liquid
         P1_kappa = fuel_cond; // Thermal conductivity of the liquid
         P1_rho = fuel_den;; // Density of the liquid
         sunrealtype y0_array[P1_NOUT];  // Array to store y1 values

            // start of thermal breakup implementation

            Problem1(y0_array);

            // Copy the values from sundials into a Converge vector array
            for (int i = 0; i < P1_NOUT; i++) {
               original_y0_array[i] = y0_array[i];
            }

            // Print the original and new values for comparison while troubleshooting
            for (int i = 0; i < P1_NOUT; i++) {
               printf("original_y0_array[%d] = %12.5e \n", i, original_y0_array[i]);
            }

            // Set the last iteration in the vector array to a new Converge variable
            radius_bubble = original_y0_array[9]; // setting radius to the final iteration of the bubble growth ODE

            Fp = 4 * PI * (P1_Pv * CONVERGE_pow(radius_bubble,2) - P1_pAmbient * CONVERGE_pow(parcel_cloud.radius[passed_parcel_idx],2));
            rMean = (radius_bubble + parcel_cloud.radius[passed_parcel_idx]) / 2; // Mean radius
            Fs = 2 * PI * rMean * P1_sigma; // Surface tension force

            if (Fp > Fs) {
               // Breakup occurs
               new_radius_thermal = cbrt(0.5 * (parcel_cloud.radius[passed_parcel_idx] - radius_bubble)); // New radius
               // new_radius_thermal = fmax(new_radius_thermal, 1.0e-20);
               old_num_drop = parcel_cloud.num_drop[passed_parcel_idx];
               parcel_cloud.num_drop[passed_parcel_idx] =
                  old_num_drop * parcel_cloud.radius[passed_parcel_idx] * parcel_cloud.radius[passed_parcel_idx] *
                  parcel_cloud.radius[passed_parcel_idx] / (new_radius_thermal * new_radius_thermal * new_radius_thermal);
               printf ("Thermal breakup worked and updated parcel properties! \n");
            }else {
               // No breakup occurs, handle accordingly
               new_radius_thermal = parcel_cloud.radius[passed_parcel_idx]; // Keep the original radius
               printf ("No thermal breakup happens \n");
            }

            // Use new_radius_thermal in the original KH loop
            new_radius = (new_radius_thermal + (dt / breakup_time) * radius_equil) / (1.0 + (dt / breakup_time));

            /*
            new_radius =
               (parcel_cloud.radius[passed_parcel_idx] + (dt / breakup_time) * radius_equil) / (1.0 + (dt / breakup_time));
         }
         else
         {
            new_radius = parcel_cloud.radius[passed_parcel_idx] - khact_c_tcav * scale_ratio * dt;
            new_radius = fmax(new_radius, 1.0e-20);
         }*/

      // new_radius = (new_radius_thermal + (dt / breakup_time) * radius_equil) / (1.0 + (dt / breakup_time));
      } else {
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


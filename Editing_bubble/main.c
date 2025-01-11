#include <cvode/cvode.h>      /* prototypes for CVODE fcts., consts.          */
#include <cvode/cvode_diag.h> /* access to CVDIAG linear solver               */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector                    */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h> /* definition of sunrealtype                       */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver               */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix                     */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                    */

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to the fixed point SUNNonlinearSolver */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to the newton SUNNonlinearSolver      */

/* helpful macros */

#ifndef SQR
#define SQR(A) ((A) * (A))
#endif

#ifndef CUB
#define CUB(A) ((A) * (A) * (A))
#endif

/* Shared Problem Constants */

#define ATOL SUN_RCONST(1.0e-12)
#define RTOL SUN_RCONST(1.0e-4)

#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define THREE  SUN_RCONST(3.0)
#define FOUR   SUN_RCONST(4.0)
#define EIGHT  SUN_RCONST(8.0)

/* Problem #1 Constants */

#define P1_NEQ        2
#define P1_ETA        SUN_RCONST(3.0)
#define P1_NOUT       50
#define P1_T0         SUN_RCONST(0.0)
#define P1_T1         SUN_RCONST(0.1)
#define P1_DTOUT      SUN_RCONST(0.2)
#define P1_TOL_FACTOR SUN_RCONST(1.0e7)
#define P1_RBInit     SUN_RCONST(0.0215e-6)
#define P1_RBDotInit  SUN_RCONST(0.1)
#define P1_CL         SUN_RCONST(1300)
#define P1_Pv         SUN_RCONST(1160000)
#define P1_pAmbient   SUN_RCONST(116000)
#define P1_sigma      SUN_RCONST(0.25)
#define P1_Pr         SUN_RCONST(0.71)
#define P1_Pr0        SUN_RCONST(0.71)
#define P1_mu         SUN_RCONST(0.24e-3)
#define P1_kappa      SUN_RCONST(0.1)
#define P1_rho        SUN_RCONST(595.59)

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

static int Problem1(void);
static void PrintIntro1(void);
static void PrintHeader1(void);
static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu,
sunrealtype hu);
static int PrepareNextRun(SUNContext sunctx, void* cvode_mem, int lmm,
int miter, N_Vector y, SUNMatrix* A, sunindextype mu,
sunindextype ml, SUNLinearSolver* LS,
SUNNonlinearSolver* NLS);
static void PrintErrOutput(sunrealtype tol_factor);
static void PrintFinalStats(void* cvode_mem, int miter, sunrealtype ero);
static void PrintErrInfo(int nerr);

/* Functions Called by the Solver */

static int f1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */

static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Implementation */

int main(void)
{
int nerr;

nerr = Problem1();
PrintErrInfo(nerr);

return (0);
}

static int Problem1(void)
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
    if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

    y = N_VNew_Serial(P1_NEQ, sunctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }
    /*PrintIntro1();*/

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    NLS = NULL;
    LS  = NULL;
    A   = NULL;

    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

    for (miter = DENSE_DQ; miter <= DENSE_DQ; miter++)
    {
        ero            = ZERO;
        NV_Ith_S(y, 0) = P1_RBInit;
        NV_Ith_S(y, 1) = P1_RBDotInit;

        firstrun = (miter == DENSE_DQ);

        retval = CVodeInit(cvode_mem, f1, P1_T0, y);
        if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        if (check_retval(&retval, "CVodeSStolerances", 1)) { return (1); }

        retval = CVodeReInit(cvode_mem, P1_T0, y);
        if (check_retval(&retval, "CVodeReInit", 1)) { return (1); }

        retval = PrepareNextRun(sunctx, cvode_mem, CV_BDF, miter, y, &A, 0, 0, &LS, &NLS);
        if (check_retval(&retval, "PrepareNextRun", 1)) { return (1); }

        /*PrintHeader1();*/

        for (iout = 1, tout = P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT)
        {
            retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
            check_retval(&retval, "CVode", 1);
            temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
            if (check_retval(&temp_retval, "CVodeGetLastOrder", 1)) { ++nerr; }
            temp_retval = CVodeGetLastStep(cvode_mem, &hu);
            if (check_retval(&temp_retval, "CVodeGetLastStep", 1)) { ++nerr; }
            PrintOutput1(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), qu, hu);
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

        PrintFinalStats(cvode_mem, miter, ero);
    }

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    N_VDestroy(y);
    SUNContext_Free(&sunctx);

    return (nerr);
}

static void PrintIntro1(void)
{
printf("Demonstration program for CVODE package - direct linear solvers\n");
printf("\n\n");
printf("Problem 1: Van der Pol oscillator\n");
printf(" xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
printf(" neq = %d,  reltol = %.2Lg,  abstol = %.2Lg", P1_NEQ, RTOL, ATOL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", P1_NEQ, RTOL, ATOL);
#else
printf(" neq = %d,  reltol = %.2g,  abstol = %.2g", P1_NEQ, RTOL, ATOL);
#endif
}

static void PrintHeader1(void)
{
printf("\n     t           x              xdot         qu     hu \n");

return;
}

static void PrintOutput1(sunrealtype t, sunrealtype y0, sunrealtype y1, int qu,
sunrealtype hu)
{
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

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

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

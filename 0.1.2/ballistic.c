/*
Ballistic: a software to benchmark ballistic models.

AUTHORS: Javier Burguete Tolosa.

Copyright 2018, AUTHORS.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY AUTHORS ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/

/**
 * \file ballistic.c
 * \brief Source file with the main function.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "equation.h"
#include "method.h"
#include "runge-kutta.h"
#include "multi-steps.h"

#define DEBUG 0                 ///< macro to debug the code.

long double convergence_factor;
///< convergence factor.
unsigned int steps;
///< number of numerical steps (1 on Runge-Kutta, >1 on multi-steps).
unsigned int ntrajectories;
///< number of projectil trajectories to calculate.
unsigned int convergence;
///< number of convergence steps.

/**
 * Function to read the basic input data.
 */
int
read_data (FILE * file,         ///< file.
           Equation * eq)       ///< Equation struct.
{
  if (!equation_read (eq, file))
    goto fail;
  if (fscanf (file, "%*s%*s%u", &ntrajectories) != 1)
    goto fail;
  if (fscanf (file, "%*s%*s%u", &convergence) != 1)
    goto fail;
  if (fscanf (file, "%*s%*s%Lf", &convergence_factor) != 1)
    goto fail;
  if (fscanf (file, "%*s%*s%u", &steps) != 1)
    goto fail;
  return 1;

fail:
  printf ("Error reading data\n");
  return 0;
}

/**
 * Function to print the solution values.
 */
void
print_solution (char *label,    ///< label.
                long double *r0,        ///< position vector.
                long double *r1)        ///< velocity vector.
{
  printf ("%s\n", label);
  printf ("x = %.19Le\n", r0[0]);
  printf ("y = %.19Le\n", r0[1]);
  printf ("z = %.19Le\n", r0[2]);
  printf ("vx = %.19Le\n", r1[0]);
  printf ("vy = %.19Le\n", r1[1]);
  printf ("vz = %.19Le\n", r1[2]);
}

/**
 * Function to print the numerical errors.
 */
void
print_error (char *label,       ///< label.
             long double *r1,   ///< numerical solution vector.
             long double *r2)   ///< analytical solution vector.
{
  printf ("%s = %.19Le\n", label, distance (r1, r2));
}

/**
 * Main function
 *
 * \return 0 on success, error code otherwise.
 */
int
main (int argn,                 ///< number of arguments.
      char **argc)              ///< array of argument chars.
{
  MultiSteps ms[1];
  RungeKutta rk[1];
  Equation eq[1];
  Method *m;
  gsl_rng *rng;
  FILE *file;
  long double sr0[3], sr1[3];
  long double t, l0r0, l2r0, l0r1, l2r1, e;
  unsigned int i, j;
#if DEBUG
  fprintf (stderr, "main: start\n");
#endif
  if (argn != 3)
    {
      printf ("The syntax is:\n./runge-kutta input_file output_file\n");
      return 1;
    }
#if DEBUG
  fprintf (stderr, "main: reading data\n");
#endif
  file = fopen (argc[1], "r");
  if (!file)
    {
      printf ("Unable to open the input file\n");
      return 2;
    }
  if (!read_data (file, eq))
    return 3;
#if DEBUG
  fprintf (stderr, "main: initing method\n");
#endif
  if (steps == 1)
    {
      if (!runge_kutta_read (rk, file))
        {
          fclose (file);
          return 4;
        }
      runge_kutta_init_variables (rk);
      m = RUNGE_KUTTA_METHOD (rk);
    }
  else
    {
      if (!multi_steps_init (ms, steps) || !multi_steps_read (ms, file))
        {
          fclose (file);
          return 4;
        }
      multi_steps_init_variables (ms);
      m = MULTI_STEPS_METHOD (ms);
    }
  fclose (file);
  rng = gsl_rng_alloc (gsl_rng_taus);
  file = fopen (argc[2], "w");
  for (j = 0; j < convergence; ++j)
    {
      gsl_rng_set (rng, 0l);
      nevaluations = 0l;
      l0r0 = l2r0 = l0r1 = l2r1 = 0.L;
      for (i = 0; i < ntrajectories; ++i)
        {
#if DEBUG
          fprintf (stderr, "main: initing equation data\n");
#endif
          equation_init (eq, rng);
#if DEBUG
          fprintf (stderr, "main: initing variables\n");
#endif
          equation_solution (eq, r0, r1, 0.);
          equation_acceleration (eq, r0, r1, r2, 0.L);
#if DEBUG
          fprintf (stderr, "main: running\n");
#endif
          if (steps == 1)
            t = runge_kutta_run (rk, eq);
          else
            t = multi_steps_run (ms, eq);
#if DEBUG
          fprintf (stderr, "main: solutions\n");
          print_solution ("Numerical solution", r0, r1);
          printf ("Time = %.19Le\n", t);
#endif
          switch (eq->land_type)
            {
            case 0:
              equation_solution (eq, sr0, sr1, eq->tf);
              break;
            default:
              t = equation_solve (eq, sr0, sr1);
            }
#if DEBUG
          print_solution ("Analytical solution", sr0, sr1);
          printf ("Time = %.19Le\n", t);
          print_error ("Position error", r0, sr0);
          print_error ("Velocity error", r1, sr1);
#endif
          e = distance (r0, sr0);
          l0r0 = fmaxl (l0r0, e);
          l2r0 += e * e;
          e = distance (r1, sr1);
          l0r1 = fmaxl (l0r1, e);
          l2r1 += e * e;

        }
#if DEBUG
      fprintf (stderr, "main: saving results\n");
#endif
      l2r0 = sqrtl (l2r0 / ntrajectories);
      l2r1 = sqrtl (l2r1 / ntrajectories);
      fprintf (file, "%lu %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le\n",
               nevaluations, l0r0, l2r0, l0r1, l2r1, kt, m->emt);
      switch (eq->size_type)
        {
        case 0:
          dt *= convergence_factor;
          break;
        default:
          kt *= convergence_factor;
        }
      m->emt *= convergence_factor;
      if (steps > 1)
        RUNGE_KUTTA_METHOD (MULTI_STEPS_RUNGE_KUTTA (ms))->emt
          *= convergence_factor;
    }
  fclose (file);
#if DEBUG
  fprintf (stderr, "main: deleting method\n");
#endif
  if (steps == 1)
    runge_kutta_delete (rk);
  else
    multi_steps_delete (ms);
  gsl_rng_free (rng);
#if DEBUG
  fprintf (stderr, "main: end\n");
#endif
  return 0;
}

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
 * \file runge-kutta.c
 * \brief Source file to define the Runge-Kutta method data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "equation.h"
#include "method.h"
#include "runge-kutta.h"

#define DEBUG_RUNGE_KUTTA 0     ///< macro to debug the Runge-Kutta functions.

///> 1st array of 1st order Runge-Kutta a coefficients.
static const long double rk_a1_1[1] = { 1.L };

///> 1st array of 1st order Runge-Kutta c coefficients.
static const long double rk_c1_1[1] = { 1.L };

///> matrix of 1st order Runge-Kutta a coefficients.
static const long double *rk_a1[1] = { rk_a1_1 };

///> matrix of 1st order Runge-Kutta c coefficients.
static const long double *rk_c1[1] = { rk_c1_1 };

///> array of 1st order Runge-Kutta t coefficients.
static const long double rk_t1[1] = { 1.L };

///> array of 1st order Runge-Kutta error coefficients.
static const long double rk_e1[1] = { -1.L };


///> 1st array of 1st order Runge-Kutta a coefficients.
static const long double rk_a2_1[1] = { 1.L };

///> 1st array of 1st order Runge-Kutta c coefficients.
static const long double rk_c2_1[1] = { 1.L };

///> 2nd array of 2nd order Runge-Kutta a coefficients.
static const long double rk_a2_2[2] = { 0.5L, 0.5L };

///> 2nd array of 2nd order Runge-Kutta c coefficients.
static const long double rk_c2_2[2] = { 0.L, 1.L };

///> matrix of 2nd order Runge-Kutta a coefficients.
static const long double *rk_a2[2] = { rk_a2_1, rk_a2_2 };

///> matrix of 2nd order Runge-Kutta c coefficients.
static const long double *rk_c2[2] = { rk_c2_1, rk_c2_2 };

///> array of 2nd order Runge-Kutta t coefficients.
static const long double rk_t2[2] = { 1.L, 1.L };

///> array of 2nd order Runge-Kutta error coefficients.
static const long double rk_e2[2] = { 0.5L, -0.5L };


///> 1st array of 3rd order Runge-Kutta a coefficients.
static const long double rk_a3_1[1] = { 1.L };

///> 1st array of 3rd order Runge-Kutta c coefficients.
static const long double rk_c3_1[1] = { 1.L };

///> 2nd array of 3rd order Runge-Kutta a coefficients.
static const long double rk_a3_2[2] = { 0.75L, 0.25L };

///> 2nd array of 3rd order Runge-Kutta c coefficients.
static const long double rk_c3_2[2] = { 0.L, 1.L };

///> 3rd array of 3rd order Runge-Kutta a coefficients.
static const long double rk_a3_3[3] = { 1.L / 3.L, 0.L, 2.L / 3.L };

///> 3rd array of 3rd order Runge-Kutta c coefficients.
static const long double rk_c3_3[3] = { 0.L, 0.L, 1.L };

///> matrix of 3rd order Runge-Kutta a coefficients.
static const long double *rk_a3[3] = { rk_a3_1, rk_a3_2, rk_a3_3 };

///> matrix of 3rd order Runge-Kutta c coefficients.
static const long double *rk_c3[3] = { rk_c3_1, rk_c3_2, rk_c3_3 };

///> array of 3rd order Runge-Kutta t coefficients.
static const long double rk_t3[3] = { 1.L, 0.5L, 1.L };

///> array of 3rd order Runge-Kutta error coefficients.
static const long double rk_e3[3] = { 1.L / 12.L, 1.L / 12.L, -1.L / 6.L };


///> 1st array of 4th order Runge-Kutta a coefficients.
static const long double rk_a4_1[1] = { 1.L };

///> 1st array of 4th order Runge-Kutta c coefficients.
static const long double rk_c4_1[1] = { 0.5L };

///> 2nd array of 4th order Runge-Kutta a coefficients.
static const long double rk_a4_2[2] = { 0.5L, 0.5L };

///> 2nd array of 4th order Runge-Kutta c coefficients.
static const long double rk_c4_2[2] = { -0.5L, 1.L };

///> 3rd array of 4th order Runge-Kutta a coefficients.
static const long double rk_a4_3[3] = { 0.5L, -0.5L, 1.L };

///> 3rd array of 4th order Runge-Kutta c coefficients.
static const long double rk_c4_3[3] = { 0.5L, 1.L, 1.L };

///> 4th array of 4th order Runge-Kutta a coefficients.
static const long double rk_a4_4[4] =
  { 5.L / 12.L, 0.25L, 1.L / 6.L, 1.L / 6.L };
///> 4th array of 4th order Runge-Kutta c coefficients.
static const long double rk_c4_4[4] = { 0.1L, 1.L, 1.L, 1.L };

///> matrix of 4th order Runge-Kutta a coefficients.
static const long double *rk_a4[4] = { rk_a4_1, rk_a4_2, rk_a4_3, rk_a4_4 };

///> matrix of 4th order Runge-Kutta c coefficients.
static const long double *rk_c4[4] = { rk_c4_1, rk_c4_2, rk_c4_3, rk_c4_4 };

///> array of 4th order Runge-Kutta t coefficients.
static const long double rk_t4[4] = { 0.5L, 0.5L, 1.L, 1.L };


/**
 * Function to init the coefficients of the 1st order Runge-Kutta method.
 */
void
runge_kutta_init_1 (RungeKutta * rk)    ///< RungeKutta struct.
{
  method_init (RUNGE_KUTTA_METHOD (rk), 1, 1);
  rk->a = rk_a1;
  rk->c = rk_c1;
  rk->t = rk_t1;
  rk->e = rk_e1;
}

/**
 * Function to init the coefficients of the 2nd order Runge-Kutta method.
 */
void
runge_kutta_init_2 (RungeKutta * rk)    ///< RungeKutta struct.
{
  method_init (RUNGE_KUTTA_METHOD (rk), 2, 2);
  rk->a = rk_a2;
  rk->c = rk_c2;
  rk->t = rk_t2;
  rk->e = rk_e2;
}

/**
 * Function to init the coefficients of the 3rd order Runge-Kutta method.
 */
void
runge_kutta_init_3 (RungeKutta * rk)    ///< RungeKutta struct.
{
  method_init (RUNGE_KUTTA_METHOD (rk), 3, 3);
  rk->a = rk_a3;
  rk->c = rk_c3;
  rk->t = rk_t3;
  rk->e = rk_e3;
}

/**
 * Function to init the coefficients of the 4th order Runge-Kutta method.
 */
void
runge_kutta_init_4 (RungeKutta * rk)    ///< RungeKutta struct.
{
  method_init (RUNGE_KUTTA_METHOD (rk), 4, 4);
  rk->a = rk_a4;
  rk->c = rk_c4;
  rk->t = rk_t4;
}

/**
 * Function to init the variables used by the Runge-Kutta methods.
 */
void
runge_kutta_init_variables (RungeKutta * rk)
{
  method_init_variables (RUNGE_KUTTA_METHOD (rk));
}

/**
 * Function to perform a step of the Runge-Kutta method.
 */
void
runge_kutta_step (RungeKutta * rk,      ///< RungeKutta struct.
                  Equation * eq,        ///< Equation struct.
                  long double t,        ///< actual time.
                  long double dt)       ///< time step size.
{
  Method *m;
  const long double *a;
  const long double *c;
  unsigned int i, j, n;
#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_step: start\n");
#endif
  m = RUNGE_KUTTA_METHOD (rk);
  memcpy (m->r0[0], r0, 3 * sizeof (long double));
  memcpy (m->r1[0], r1, 3 * sizeof (long double));
  memcpy (m->r2[0], r2, 3 * sizeof (long double));
#if DEBUG_RUNGE_KUTTA
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "runge_kutta_step: r0[0][%u]=%Lg\n", i, m->r0[0][i]);
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "runge_kutta_step: r1[0][%u]=%Lg\n", i, m->r1[0][i]);
#endif
  n = m->nsteps;
  for (i = 1; i <= n; ++i)
    {
      a = rk->a[i - 1];
      c = rk->c[i - 1];
#if DEBUG_RUNGE_KUTTA
      for (j = 0; j < i; ++j)
        fprintf (stderr, "runge_kutta_step: a%u-%u=%Lg\n", i, j, a[j]);
      for (j = 0; j < i; ++j)
        fprintf (stderr, "runge_kutta_step: c%u-%u=%Lg\n", i, j, c[j]);
#endif
      m->r0[i][0] = a[0] * (m->r0[0][0] + dt * c[0] * m->r1[0][0]);
      m->r0[i][1] = a[0] * (m->r0[0][1] + dt * c[0] * m->r1[0][1]);
      m->r0[i][2] = a[0] * (m->r0[0][2] + dt * c[0] * m->r1[0][2]);
      m->r1[i][0] = a[0] * (m->r1[0][0] + dt * c[0] * m->r2[0][0]);
      m->r1[i][1] = a[0] * (m->r1[0][1] + dt * c[0] * m->r2[0][1]);
      m->r1[i][2] = a[0] * (m->r1[0][2] + dt * c[0] * m->r2[0][2]);
      for (j = 1; j < i; ++j)
        {
          m->r0[i][0] += a[j] * (m->r0[j][0] + dt * c[j] * m->r1[j][0]);
          m->r0[i][1] += a[j] * (m->r0[j][1] + dt * c[j] * m->r1[j][1]);
          m->r0[i][2] += a[j] * (m->r0[j][2] + dt * c[j] * m->r1[j][2]);
          m->r1[i][0] += a[j] * (m->r1[j][0] + dt * c[j] * m->r2[j][0]);
          m->r1[i][1] += a[j] * (m->r1[j][1] + dt * c[j] * m->r2[j][1]);
          m->r1[i][2] += a[j] * (m->r1[j][2] + dt * c[j] * m->r2[j][2]);
        }
      equation_acceleration (eq, m->r0[i], m->r1[i], m->r2[i],
                             t + rk->t[i - 1] * dt);
#if DEBUG_RUNGE_KUTTA
      fprintf (stderr, "runge_kutta_step: t%u=%Lg\n", i, rk->t[i - 1]);
#endif
    }
  --i;
  memcpy (r0, m->r0[i], 3 * sizeof (long double));
  memcpy (r1, m->r1[i], 3 * sizeof (long double));
  memcpy (r2, m->r2[i], 3 * sizeof (long double));
#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_step: end\n");
#endif
}

/**
 * Function to estimate the error on a Runge-Kutta step.
 */
void
runge_kutta_error (RungeKutta * rk,     ///< Runge-Kutta struct.
                   long double dt)      ///< time step size.
{
  long double e0[3], e1[3];
  Method *m;
  unsigned int i;
#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_error: start\n");
#endif
  m = RUNGE_KUTTA_METHOD (rk);
  e0[0] = e0[1] = e0[2] = e1[0] = e1[1] = e1[2] = 0.L;
  for (i = 0; i < m->nsteps; ++i)
    {
      e0[0] += dt * rk->e[i] * m->r1[i][0];
      e0[1] += dt * rk->e[i] * m->r1[i][1];
      e0[2] += dt * rk->e[i] * m->r1[i][2];
      e1[0] += dt * rk->e[i] * m->r2[i][0];
      e1[1] += dt * rk->e[i] * m->r2[i][1];
      e1[2] += dt * rk->e[i] * m->r2[i][2];
    }
  m->e0 = sqrtl (e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2]);
  m->e1 = sqrtl (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  m->et0 += m->e0;
  m->et1 += m->e1;
#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_error: e0=%Lg et0=%Lg\n", m->e0, m->et0);
  fprintf (stderr, "runge_kutta_error: e1=%Lg et1=%Lg\n", m->e1, m->et1);
  fprintf (stderr, "runge_kutta_error: end\n");
#endif
}

/**
 * Function to run the Runge-Kutta method bucle.
 *
 * \return final time. 
 */
long double
runge_kutta_run (RungeKutta * rk,       ///< RungeKutta struct.
                 Equation * eq) ///< Equation struct.
{
  Method *m;
  long double t, to, dt, dto, et0o, et1o;

#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_run: start\n");
#endif

  // variables backup 
  m = RUNGE_KUTTA_METHOD (rk);
  memcpy (ro0, r0, 3 * sizeof (long double));
  memcpy (ro1, r1, 3 * sizeof (long double));
  memcpy (ro2, r2, 3 * sizeof (long double));

  // temporal bucle
  for (t = 0.L; 1;)
    {

      // time step size
      if (t > 0.L && m->error_dt)
        {
          dto = dt;
          dt = method_dt (m, dt);

          // revert the step if big error
          if (dt < m->beta * dto)
            {
              t = to;
              m->et0 = et0o;
              m->et1 = et1o;
              memcpy (r0, ro0, 3 * sizeof (long double));
              memcpy (r1, ro1, 3 * sizeof (long double));
              memcpy (r2, ro2, 3 * sizeof (long double));
            }
        }
      else
        dt = equation_step_size (eq);

      // checking trajectory end
      to = t;
      if (equation_land (eq, &t, &dt))
        break;
#if DEBUG_RUNGE_KUTTA
      fprintf (stderr, "runge_kutta_run: t=%Lg dt=%Lg\n", t, dt);
#endif

      // backup of variables
      memcpy (ro0, r0, 3 * sizeof (long double));
      memcpy (ro1, r1, 3 * sizeof (long double));
      memcpy (ro2, r2, 3 * sizeof (long double));

      // Runge-Kutta step
      runge_kutta_step (rk, eq, to, dt);

      // error estimate
      if (m->error_dt)
        {
          et0o = m->et0;
          et1o = m->et1;
          runge_kutta_error (rk, dt);
        }
    }
#if DEBUG_RUNGE_KUTTA
  fprintf (stderr, "runge_kutta_run: end\n");
#endif
  return t;
}

/**
 * Function to free the memory used by a RungeKutta struct.
 */
void
runge_kutta_delete (RungeKutta * rk)    ///< RungeKutta struct.
{
  method_delete (RUNGE_KUTTA_METHOD (rk));
}

/**
 * Function to read the Runge-Kutta method data on a file.
 *
 * \return 1 on success, 0 on error.
 */
int
runge_kutta_read (RungeKutta * rk,      ///< RungeKutta struct.
                  FILE * file)  ///< file.
{
  unsigned int type;
  if (fscanf (file, "%*s%*s%u", &type) != 1)
    goto fail;
  if (!method_read (RUNGE_KUTTA_METHOD (rk), file))
    goto fail;
  switch (type)
    {
    case 1:
      runge_kutta_init_1 (rk);
      break;
    case 2:
      runge_kutta_init_2 (rk);
      break;
    case 3:
      runge_kutta_init_3 (rk);
      break;
    case 4:
      runge_kutta_init_4 (rk);
      break;
    default:
      goto fail;
    }
  return 1;

fail:
  printf ("Error reading Runge-Kutta data\n");
  return 0;
}

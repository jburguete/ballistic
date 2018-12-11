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
 * \file multi-steps.c
 * \brief Source file to define the multi-steps method data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <libxml/parser.h>
#include "config.h"
#include "utils.h"
#include "equation.h"
#include "method.h"
#include "runge-kutta.h"
#include "multi-steps.h"

#define DEBUG_MULTI_STEPS 0     ///< macro to debug the multi-steps functions.

///> array of a coefficients of the 2nd order multi-steps method. 
const long double ms_a2[3] = { 0.75L, 0.L, 0.25L };

///> array of c coefficients of the 2nd order multi-steps method. 
const long double ms_c2[3] = { 2.L, 0.L, 0.L };

///> array of error a coefficients of the 2nd order multi-steps method. 
const long double ms_ea2[3] = { 0.25L, 0.L, -0.25L };

///> array of error b coefficients of the 2nd order multi-steps method. 
const long double ms_eb2[3] = { 0.5L, 0.L, 0.L };


///> array of a coefficients of the 3rd order multi-steps method. 
const long double ms_a3[4] = { 16.L / 27.L, 0.L, 0.L, 11.L / 27.L };

///> array of c coefficients of the 3rd order multi-steps method. 
const long double ms_c3[4] = { 3.L, 0.L, 0.L, 12.L / 11.L };

///> array of error a coefficients of the 3rd order multi-steps method. 
const long double ms_ea3[4] = { 17.L / 108.L, 0.L, 0.25L, -11.L / 27.L };

///> array of error b coefficients of the 3rd order multi-steps method. 
const long double ms_eb3[4] = { -5.L / 18.L, 0.L, 0.L, -4.L / 9.L };

/**
 * Function to init the multi-steps method coefficients.
 *
 * \return 1 on success, 0 on error.
 */
int
multi_steps_init (MultiSteps * ms,      ///< MultiSteps struct.
                  unsigned int nsteps)  ///< number of steps.
{
  Method *m;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_init: start\n");
#endif
  m = MULTI_STEPS_METHOD (ms);
  switch (nsteps)
    {
    case 2:
      method_init (m, 3, 2);
      ms->a = ms_a2;
      ms->c = ms_c2;
      ms->ea = ms_ea2;
      ms->eb = ms_eb2;
      break;
    case 3:
      method_init (m, 4, 3);
      ms->a = ms_a3;
      ms->c = ms_c3;
      ms->ea = ms_ea3;
      ms->eb = ms_eb3;
      break;
    default:
      printf ("Error reading multi-steps data\n");
#if DEBUG_MULTI_STEPS
      fprintf (stderr, "multi_steps_init: error\n");
      fprintf (stderr, "multi_steps_init: end\n");
#endif
      return 0;
    }
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_init: success\n");
  fprintf (stderr, "multi_steps_init: end\n");
#endif
  return 1;
}

/**
 * Function to init the variables used by the multi-steps methods.
 */
void
multi_steps_init_variables (MultiSteps * ms)
{
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_init_variables: start\n");
#endif
  runge_kutta_init_variables (MULTI_STEPS_RUNGE_KUTTA (ms));
  method_init_variables (MULTI_STEPS_METHOD (ms));
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_init_variables: end\n");
#endif
}

/**
 * Function to perform a step of the multi-steps method.
 */
static void
multi_steps_step (MultiSteps * ms,      ///< MultiSteps struct.
                  Equation * eq,        ///< Equation struct.
                  long double t,        ///< actual time.
                  long double dt)       ///< time step size.
{
  long double msr0[3], msr1[3];
  Method *m;
  const long double *a;
  const long double *c;
  unsigned int i, n;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_step: start\n");
#endif
  m = MULTI_STEPS_METHOD (ms);
#if DEBUG_MULTI_STEPS
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r0[0][%u]=%Lg\n", i, m->r0[0][i]);
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r1[0][%u]=%Lg\n", i, m->r1[0][i]);
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r2[0][%u]=%Lg\n", i, m->r2[0][i]);
#endif
  n = m->nsteps;
  a = ms->a;
  c = ms->c;
#if DEBUG_MULTI_STEPS
  for (i = 0; i < n; ++i)
    fprintf (stderr, "multi_steps_step: a%u=%Lg\n", i, a[i]);
  for (i = 0; i < n; ++i)
    fprintf (stderr, "multi_steps_step: c%u=%Lg\n", i, c[i]);
#endif
  msr0[0] = a[0] * (r0[0] + dt * c[0] * r1[0]);
  msr0[1] = a[0] * (r0[1] + dt * c[0] * r1[1]);
  msr0[2] = a[0] * (r0[2] + dt * c[0] * r1[2]);
  msr1[0] = a[0] * (r1[0] + dt * c[0] * r2[0]);
  msr1[1] = a[0] * (r1[1] + dt * c[0] * r2[1]);
  msr1[2] = a[0] * (r1[2] + dt * c[0] * r2[2]);
  for (i = 1; i < n; ++i)
    {
#if DEBUG_MULTI_STEPS
      fprintf (stderr, "multi_steps_step: r0[%u][0]=%Lg\n", i, m->r0[i][0]);
      fprintf (stderr, "multi_steps_step: r0[%u][1]=%Lg\n", i, m->r0[i][1]);
      fprintf (stderr, "multi_steps_step: r0[%u][2]=%Lg\n", i, m->r0[i][2]);
      fprintf (stderr, "multi_steps_step: r1[%u][0]=%Lg\n", i, m->r1[i][0]);
      fprintf (stderr, "multi_steps_step: r1[%u][1]=%Lg\n", i, m->r1[i][1]);
      fprintf (stderr, "multi_steps_step: r1[%u][2]=%Lg\n", i, m->r1[i][2]);
      fprintf (stderr, "multi_steps_step: r2[%u][0]=%Lg\n", i, m->r2[i][0]);
      fprintf (stderr, "multi_steps_step: r2[%u][1]=%Lg\n", i, m->r2[i][1]);
      fprintf (stderr, "multi_steps_step: r2[%u][2]=%Lg\n", i, m->r2[i][2]);
#endif
      msr0[0] += a[i] * (m->r0[i][0] + dt * c[i] * m->r1[i][0]);
      msr0[1] += a[i] * (m->r0[i][1] + dt * c[i] * m->r1[i][1]);
      msr0[2] += a[i] * (m->r0[i][2] + dt * c[i] * m->r1[i][2]);
      msr1[0] += a[i] * (m->r1[i][0] + dt * c[i] * m->r2[i][0]);
      msr1[1] += a[i] * (m->r1[i][1] + dt * c[i] * m->r2[i][1]);
      msr1[2] += a[i] * (m->r1[i][2] + dt * c[i] * m->r2[i][2]);
    }
#if DEBUG_MULTI_STEPS
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r0[0][%u]=%Lg\n", i, msr0[i]);
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r1[0][%u]=%Lg\n", i, msr1[i]);
#endif
  memcpy (m->r0[0], r0, 3 * sizeof (long double));
  memcpy (m->r1[0], r1, 3 * sizeof (long double));
  memcpy (m->r2[0], r2, 3 * sizeof (long double));
  for (i = n; --i > 0;)
    {
      memcpy (m->r0[i], m->r0[i - 1], 3 * sizeof (long double));
      memcpy (m->r1[i], m->r1[i - 1], 3 * sizeof (long double));
      memcpy (m->r2[i], m->r2[i - 1], 3 * sizeof (long double));
    }
  memcpy (r0, msr0, 3 * sizeof (long double));
  memcpy (r1, msr1, 3 * sizeof (long double));
  equation_acceleration (eq, r0, r1, r2, t + dt);
#if DEBUG_MULTI_STEPS
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r0[0][%u]=%Lg\n", i, r0[i]);
  for (i = 0; i < 3; ++i)
    fprintf (stderr, "multi_steps_step: r1[0][%u]=%Lg\n", i, r1[i]);
  fprintf (stderr, "multi_steps_step: end\n");
#endif
}

/**
 * Function to estimate the error on a multi-steps method step.
 */
static void
multi_steps_error (MultiSteps * ms,     ///< MultiSteps struct.
                   long double dt)      ///< time step size.
{
  long double e0[3], e1[3];
  Method *m;
  unsigned int i;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_error: start\n");
#endif
  m = MULTI_STEPS_METHOD (ms);
  e0[0] = e0[1] = e0[2] = e1[0] = e1[1] = e1[2] = 0.L;
  for (i = 0; i < m->nsteps; ++i)
    {
      e0[0] += ms->ea[i] * m->r0[i][0] + dt * ms->eb[i] * m->r1[i][0];
      e0[1] += ms->ea[i] * m->r0[i][1] + dt * ms->eb[i] * m->r1[i][1];
      e0[2] += ms->ea[i] * m->r0[i][2] + dt * ms->eb[i] * m->r1[i][2];
      e1[0] += ms->ea[i] * m->r1[i][0] + dt * ms->eb[i] * m->r2[i][0];
      e1[1] += ms->ea[i] * m->r1[i][1] + dt * ms->eb[i] * m->r2[i][1];
      e1[2] += ms->ea[i] * m->r1[i][2] + dt * ms->eb[i] * m->r2[i][2];
    }
  m->e0 = sqrtl (e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2]);
  m->e1 = sqrtl (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  m->et0 += m->e0;
  m->et1 += m->e1;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_error: e0=%Lg et0=%Lg\n", m->e0, m->et0);
  fprintf (stderr, "multi_steps_error: e1=%Lg et1=%Lg\n", m->e1, m->et1);
  fprintf (stderr, "multi_steps_error: end\n");
#endif
}

/**
 * Function to run the multi-steps method bucle.
 *
 * \return final time. 
 */
long double
multi_steps_run (MultiSteps * ms,       ///< MultiSteps struct.
                 Equation * eq) ///< Equation struct.
{
  RungeKutta *rk;
  Method *m, *mrk;
  long double t, dt, to, dto, et0o, et1o;
  unsigned int i, n;

#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_run: start\n");
#endif

  // variables backup 
  rk = MULTI_STEPS_RUNGE_KUTTA (ms);
  mrk = RUNGE_KUTTA_METHOD (rk);
  m = MULTI_STEPS_METHOD (ms);
  memcpy (ro0, r0, 3 * sizeof (long double));
  memcpy (ro1, r1, 3 * sizeof (long double));
  memcpy (ro2, r2, 3 * sizeof (long double));

  // Runge-Kutta first steps
  n = m->nsteps;
  for (t = 0.L, i = n; --i > 0;)
    {

      // time step size
      if (t > 0.L && mrk->error_dt)
        {
          dto = dt;
          dt = method_dt (mrk, dt);

          // revert the step if big error
          if (dt < mrk->beta * dto)
            {
              ++i;
              t = to;
              mrk->et0 = et0o;
              mrk->et1 = et1o;
              memcpy (r0, ro0, 3 * sizeof (long double));
              memcpy (r1, ro1, 3 * sizeof (long double));
              memcpy (r2, ro2, 3 * sizeof (long double));
            }
        }
      else
        dt = equation_step_size (eq);

      // checking trajectory end
      to = t;
      if (equation_land (eq, to, &t, &dt))
        goto end;
#if DEBUG_MULTI_STEPS
      fprintf (stderr, "multi_steps_run: t=%Lg dt=%Lg\n", t, dt);
#endif

      // saving step 
      memcpy (m->r0[i], r0, 3 * sizeof (long double));
      memcpy (m->r1[i], r1, 3 * sizeof (long double));
      memcpy (m->r2[i], r2, 3 * sizeof (long double));
      memcpy (ro0, r0, 3 * sizeof (long double));
      memcpy (ro1, r1, 3 * sizeof (long double));
      memcpy (ro2, r2, 3 * sizeof (long double));

      // Runge-Kutta step
      runge_kutta_step (rk, eq, to, dt);

      // error estimate
      if (mrk->error_dt)
        {
          et0o = mrk->et0;
          et1o = mrk->et1;
          runge_kutta_error (rk, dt);
        }
    }

  // saving last step 
  memcpy (m->r0[0], r0, 3 * sizeof (long double));
  memcpy (m->r1[0], r1, 3 * sizeof (long double));
  memcpy (m->r2[0], r2, 3 * sizeof (long double));

  // initing errors
  m->et0 = mrk->et0;
  m->et1 = mrk->et1;

  // temporal bucle
  while (1)
    {

      // time step size
      if (m->error_dt)
        {
          dto = dt;
          dt = method_dt (m, dt);

          // revert the step if big error
          if (dt < mrk->beta * dto)
            {
              ++i;
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
      dto = dt;
      if (equation_land (eq, to, &t, &dt))
        break;
#if DEBUG_MULTI_STEPS
      fprintf (stderr, "multi_steps_run: t=%Lg dt=%Lg\n", t, dt);
#endif

      // variables backup
      memcpy (ro0, r0, 3 * sizeof (long double));
      memcpy (ro1, r1, 3 * sizeof (long double));
      memcpy (ro2, r2, 3 * sizeof (long double));

      // multi-steps step
      if (dto == dt)
        multi_steps_step (ms, eq, to, dt);
      else
        runge_kutta_step (rk, eq, to, dt);

      // error estimate
      if (m->error_dt)
        {
          et0o = m->et0;
          et1o = m->et1;
          multi_steps_error (ms, dt);
        }
    }
end:
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_run: end\n");
#endif
  return t;
}

/**
 * Function to free the memory used by a MultiSteps struct.
 */
void
multi_steps_delete (MultiSteps * ms)    ///< RungeKutta struct.
{
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_delete: start\n");
#endif
  runge_kutta_delete (MULTI_STEPS_RUNGE_KUTTA (ms));
  method_delete (MULTI_STEPS_METHOD (ms));
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_delete: end\n");
#endif
}

/**
 * Function to read on a file the multi-steps method input data.
 *
 * \return 1 on success, 0 on error.
 */
int
multi_steps_read (MultiSteps * ms,      ///< MultiSteps struct.
                  FILE * file)  ///< input file.
{
  int e;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_read: start\n");
#endif
  if (!method_read (MULTI_STEPS_METHOD (ms), file))
    e = 0;
  else
    e = runge_kutta_read (MULTI_STEPS_RUNGE_KUTTA (ms), file);
#if DEBUG_MULTI_STEPS
  if (e)
    fprintf (stderr, "multi_steps_read: success\n");
  else
    fprintf (stderr, "multi_steps_read: error\n");
  fprintf (stderr, "multi_steps_read: end\n");
#endif
  return e;
}

/**
 * Function to read on a XML node the multi-steps method input data.
 *
 * \return 1 on success, 0 on error.
 */
int
multi_steps_read_xml (MultiSteps * ms,  ///< MultiSteps struct.
                      xmlNode * node)   ///< XML node.
{
  int e, error_code;
#if DEBUG_MULTI_STEPS
  fprintf (stderr, "multi_steps_read_xml: start\n");
#endif
  if (xmlStrcmp (node->name, XML_MULTI_STEPS) || !node->children)
    e = 0;
  else
    {
      ms->steps = xml_node_get_uint (node, XML_STEPS, &error_code);
      if (error_code || !method_read_xml (MULTI_STEPS_METHOD (ms), node))
        e = 0;
      else
        e = runge_kutta_read_xml (MULTI_STEPS_RUNGE_KUTTA (ms), node->children);
    }
#if DEBUG_MULTI_STEPS
  if (e)
    fprintf (stderr, "multi_steps_read: success\n");
  else
    fprintf (stderr, "multi_steps_read: error\n");
  fprintf (stderr, "multi_steps_read: end\n");
#endif
  return e;
}

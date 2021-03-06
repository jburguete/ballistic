/*
Ballistic: a software to benchmam ballistic models.

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
 * \brief Source file to define the numerical method data and functions.
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
#include <glib.h>
#include "config.h"
#include "utils.h"
#include "equation.h"
#include "method.h"

#define DEBUG_METHOD 0
///< macro to debug the numerical method functions.

/**
 * Function to init the numerical method.
 */
void
method_init (Method * m,        ///< Method struct.
             unsigned int nsteps,       ///< number of steps.
             unsigned int order)        ///< order of accuracy.
{
#if DEBUG_METHOD
  fprintf (stderr, "method_init: start\n");
#endif
  m->nsteps = nsteps;
  m->order = order;
#if DEBUG_METHOD
  fprintf (stderr, "method_init: end\n");
#endif
}

/**
 * Function to init the variables of the numerical method.
 */
void
method_init_variables (Method * m)      ///< Method struct.
{
  unsigned int i, n;
#if DEBUG_METHOD
  fprintf (stderr, "method_init_variables: start\n");
  fprintf (stderr, "method_init_variables: nsteps=%u\n", m->nsteps);
#endif
  n = m->nsteps + 1;
  m->r0 = (long double **) g_slice_alloc (n * sizeof (long double *));
  m->r1 = (long double **) g_slice_alloc (n * sizeof (long double *));
  m->r2 = (long double **) g_slice_alloc (n * sizeof (long double *));
  for (i = 0; i < n; ++i)
    {
      m->r0[i] = (long double *) g_slice_alloc (3 * sizeof (long double));
      m->r1[i] = (long double *) g_slice_alloc (3 * sizeof (long double));
      m->r2[i] = (long double *) g_slice_alloc (3 * sizeof (long double));
    }
  m->et0 = m->et1 = 0.L;
#if DEBUG_METHOD
  fprintf (stderr, "method_init_variables: end\n");
#endif
}

/**
 * Function to calculate the following numerical step size based on error
 * control.
 *
 * \return next time step size.
 */
long double
method_dt (Method * m,          ///< Method struct.
           long double dt)      ///< actual time step size.
{
  long double dt2;
#if DEBUG_METHOD
  fprintf (stderr, "method_dt: start\n");
#endif
  dt2 = dt * fminl (m->alpha,
                    powl (m->emt * dt / m->e0, 1.L / (m->order - 1.L)));
#if DEBUG_METHOD
  fprintf (stderr, "method_dt: dt=%Lg\n", dt2);
  fprintf (stderr, "method_dt: end\n");
#endif
  return dt2;
}

/**
 * Function to free the memory used by the Method struct.
 */
void
method_delete (Method * m)      ///< Method struct.
{
  unsigned int i, n;
#if DEBUG_METHOD
  fprintf (stderr, "method_delete: start\n");
  fprintf (stderr, "method_delete: nsteps=%u\n", m->nsteps);
#endif
  n = i = m->nsteps + 1;
  do
    {
      --i;
      g_slice_free1 (3 * sizeof (long double), m->r2[i]);
      g_slice_free1 (3 * sizeof (long double), m->r1[i]);
      g_slice_free1 (3 * sizeof (long double), m->r0[i]);
    }
  while (i);
  g_slice_free1 (n * sizeof (long double *), m->r2);
  g_slice_free1 (n * sizeof (long double *), m->r1);
  g_slice_free1 (n * sizeof (long double *), m->r0);
#if DEBUG_METHOD
  fprintf (stderr, "method_delete: end\n");
#endif
}

/**
 * Function to read the numerical method data on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
method_read_xml (Method * m,    ///< Method struct.
                 xmlNode * node)        ///< XML node.
{
  const char *message[] = {
    "Bad dt",
    "Bad alpha",
    "Bad beta",
    "Bad error per time",
    "Unknown error control type"
  };
  int e, error_code;
#if DEBUG_METHOD
  fprintf (stderr, "method_read_xml: start\n");
#endif
  m->error_dt = xml_node_get_uint (node, XML_TIME_STEP, &error_code);
  if (error_code)
    {
      e = 0;
      goto fail;
    }
  switch (m->error_dt)
    {
    case 0:
      m->emt = 0.L;
      break;
    case 1:
      m->alpha = xml_node_get_float (node, XML_ALPHA, &error_code);
      if (error_code)
        {
          e = 1;
          goto fail;
        }
      m->beta = xml_node_get_float (node, XML_BETA, &error_code);
      if (error_code)
        {
          e = 2;
          goto fail;
        }
      m->emt = xml_node_get_float (node, XML_ERROR_TIME, &error_code);
      if (error_code)
        {
          e = 3;
          goto fail;
        }
      break;
    default:
      e = 4;
      goto fail;
    }
#if DEBUG_METHOD
  fprintf (stderr, "method_read_xml: success\n");
  fprintf (stderr, "method_read_xml: end\n");
#endif
  return 1;

fail:
  error_add (message[e]);
#if DEBUG_METHOD
  fprintf (stderr, "method_read_xml: error\n");
  fprintf (stderr, "method_read_xml: end\n");
#endif
  return 0;
}

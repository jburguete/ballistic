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
#include "equation.h"
#include "method.h"

#define DEBUG_METHOD 0          ///< macro to debug the numerical method functions.

/**
 * Function to init the numerical method.
 */
void
method_init (Method * m,        ///< Method struct.
             unsigned int nsteps,       ///< number of steps.
             unsigned int order)        ///< order of accuracy.
{
  m->nsteps = nsteps;
  m->order = order;
}

/**
 * Function to init the variables of the numerical method.
 */
void
method_init_variables (Method * m)      ///< Method struct.
{
  unsigned int i, n;
  n = m->nsteps + 1;
  m->r0 = (long double **) malloc (n * sizeof (long double *));
  m->r1 = (long double **) malloc (n * sizeof (long double *));
  m->r2 = (long double **) malloc (n * sizeof (long double *));
  for (i = 0; i < n; ++i)
    {
      m->r0[i] = (long double *) malloc (3 * sizeof (long double));
      m->r1[i] = (long double *) malloc (3 * sizeof (long double));
      m->r2[i] = (long double *) malloc (3 * sizeof (long double));
    }
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
  dt2 = dt * fmaxl (m->beta,
                    fminl (m->alpha,
                           powl (m->emt * dt / m->e0, 1.L / m->order)));
  return dt2;
}

/**
 * Function to free the memory used by the Method struct.
 */
void
method_delete (Method * m)      ///< Method struct.
{
  unsigned int i;
  i = m->nsteps + 1;
  do
    {
      --i;
      free (m->r2[i]);
      free (m->r1[i]);
      free (m->r0[i]);
    }
  while (i);
  free (m->r2);
  free (m->r1);
  free (m->r0);
#if DEBUG_METHOD
  fprintf (stderr, "method_delete: end\n");
#endif
}

/**
 * Function to read the numerical method data on a file.
 *
 * \return 1 on success, 0 on error.
 */
int
method_read (Method * m,        ///< Method struct.
             FILE * file)       ///< file.
{
  if (fscanf (file, "%*s%*s%u", &m->error_dt) != 1)
    goto fail;
  switch (m->error_dt)
    {
    case 0:
      break;
    case 1:
      if (fscanf (file, "%*s%*s%Lf", &m->alpha) != 1)
        goto fail;
      if (fscanf (file, "%*s%*s%Lf", &m->beta) != 1)
        goto fail;
      if (fscanf (file, "%*s%*s%Lf", &m->emt) != 1)
        goto fail;
      break;
    default:
      goto fail;
    }
  return 1;

fail:
  printf ("Error reading method data\n");
  return 0;
}

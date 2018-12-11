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
 * \file multi-steps.h
 * \brief Header file to define the multi-steps method data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef MULTI_STEPS__H
#define MULTI_STEPS__H 1

/**
 * \struct MultiSteps
 * \brief struct to define a mult-steps method.
 */
typedef struct
{
  Method method[1];             ///< Method struct.
  RungeKutta runge_kutta[1];    ///< Runge-Kutta struct.
  const long double *a;         ///< array of a-coefficients.
  const long double *c;         ///< array of c-coefficients.
  const long double *ea;        ///< array of a-error coefficients.
  const long double *eb;        ///< array of b-error coefficients.
  unsigned int steps;           ///< steps number.
} MultiSteps;

#define MULTI_STEPS_METHOD(ms) ((Method *)ms->method)
///< macro to access to Method struct data on a MultiSteps struct.
#define MULTI_STEPS_RUNGE_KUTTA(ms) ((RungeKutta *)ms->runge_kutta)
///< macro to access to RungeKutta struct data on a MultiSteps struct.

int multi_steps_init (MultiSteps * ms, unsigned int nsteps);
void multi_steps_init_variables (MultiSteps * ms);
long double multi_steps_run (MultiSteps * ms, Equation * eq);
void multi_steps_delete (MultiSteps * ms);
int multi_steps_read (MultiSteps * ms, FILE * file);
int multi_steps_read_xml (MultiSteps * ms, xmlNode * node);

#endif

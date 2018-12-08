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
 * \file runge-kutta.h
 * \brief Header file to define the Runge-Kutta method data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef RUNGE_KUTTA__H
#define RUNGE_KUTTA__H 1

/**
 * \struct RungeKutta
 * \brief struct to define a Runge-Kutta method.
 */
typedef struct
{
  Method method[1];             ///< method struct.
  const long double **b;        ///< matrix of b-coefficients.
  const long double *t;         ///< array of t-coefficients.
  const long double *e;         ///< array of error coefficients.
} RungeKutta;

#define RUNGE_KUTTA_METHOD(rk) ((Method *)rk->method)
///< macro to access to Method struct data on a RungeKutta struct.

void runge_kutta_init_variables (RungeKutta * rk);
void runge_kutta_step (RungeKutta * rk, Equation * eq, long double t,
                       long double dt);
void runge_kutta_error (RungeKutta * rk, long double dt);
long double runge_kutta_run (RungeKutta * rk, Equation * eq);
void runge_kutta_delete (RungeKutta * rk);
int runge_kutta_read (RungeKutta * rk, FILE * file);

#endif

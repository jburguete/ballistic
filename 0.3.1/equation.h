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
 * \file equation.h
 * \brief Header file with the equation data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef EQUATION__H
#define EQUATION__H 1

/**
 * \struct Equation
 * \brief struct to define the movement equation.
 */
typedef struct
{
  long double r[3];             ///< initial position vector.
  long double v[3];             ///< initial velocity vector.
  long double w[2];             ///< wind velocity vector.
  long double g;                ///< vertical acceleration constant.
  long double tf;               ///< final time.
  long double lambda;           ///< friction coefficient.
  long double vertical_angle;   ///< initial vertical angle.
  long double max_lambda;       ///< maximum friction coefficient.
  long double min_lambda;       ///< minimum friction coefficient.
  long double max_velocity;     ///< maximum projectil velocity.
  long double min_velocity;     ///< minimum projectil velocity.
  long double max_wind;         ///< maximum wind velocity.
  unsigned int type;            ///< equation type.
  unsigned int land_type;       ///< landing type.
  unsigned int size_type;       ///< time step size type.
} Equation;

extern long double r0[3];
extern long double r1[3];
extern long double r2[3];
extern long double ro0[3];
extern long double ro1[3];
extern long double ro2[3];
extern void (*equation_acceleration) (Equation * eq, long double *r0,
                                      long double *r1, long double *r2,
                                      long double t);
extern void (*equation_solution) (Equation * eq, long double *r0,
                                  long double *r1, long double t);
extern long double (*equation_step_size) (Equation * eq);
extern int (*equation_land) (Equation * eq, long double, long double *t,
                             long double *dt);
extern long double kt;
extern long double dt;
extern unsigned long int nevaluations;

long double equation_solve (Equation * eq, long double *r0, long double *r1);
void equation_init (Equation * eq, gsl_rng * rng);
int equation_read (Equation * eq, FILE * file);
int equation_read_xml (Equation * eq, xmlNode * node);

#endif

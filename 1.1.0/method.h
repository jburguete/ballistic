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
 * \file method.h
 * \brief Header file to define the base of numerical method data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef METHOD__H
#define METHOD__H 1

/**
 * \struct Method
 * \brief struct to define the base of a numerical method.
 */
typedef struct
{
  long double **r0;             ///< array of position vectors.
  long double **r1;             ///< array of velocity vectors.
  long double **r2;             ///< array of acceleration vectors.
  long double e0;               ///< step position error.
  long double e1;               ///< step velocity error.
  long double et0;              ///< total position error.
  long double et1;              ///< total velocity error.
  long double emt;              ///< maximum error per time.
  long double alpha;            ///< error time step size alpha parameter.
  long double beta;             ///< error time step size beta parameter.
  unsigned int nsteps;          ///< number of steps.
  unsigned int order;           ///< order.
  unsigned int error_dt;        ///< type of error time step size control. 
} Method;

void method_init (Method * m, unsigned int nsteps, unsigned int order);
void method_init_variables (Method * m);
long double method_dt (Method * m, long double dt);
void method_delete (Method * m);
int method_read_xml (Method * m, xmlNode * node);

#endif

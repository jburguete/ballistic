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
 * \file utils.c
 * \brief Source file with the useful data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <libxml/parser.h>
#include "config.h"
#include "equation.h"

#define DEBUG_UTILS 1           ///< macro to debug the useful functions.

char *error_message = NULL;     ///< error message string.

/**
 * Function to calculate the distance between two vectors.
 *
 * \return vectors distance.
 */
long double
distance (long double *r1,      ///< 1st vector.
          long double *r2)      ///< 2nd vector.
{
  long double d, dr[3];
#if DEBUG_EQUATION
  fprintf (stderr, "distance: start\n");
  fprintf (stderr, "distance: x1=%Lg y1=%Lg z1=%Lg\n", r1[0], r1[1], r1[2]);
  fprintf (stderr, "distance: x2=%Lg y2=%Lg z2=%Lg\n", r2[0], r2[1], r2[2]);
#endif
  dr[0] = r1[0] - r2[0];
  dr[1] = r1[1] - r2[1];
  dr[2] = r1[2] - r2[2];
  d = sqrtl (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
#if DEBUG_EQUATION
  fprintf (stderr, "distance: d=%Lg\n", d);
  fprintf (stderr, "distance: end\n");
#endif
  return d;
}

/**
 * Function to calculate the solution of a reduced 2nd order equation.
 *
 * This function calculates the solution of a reduced 2nd order equation in the
 * form:
 * \f[x^2+a\,x+b=0\f]
 * in the interval \f$x\in\left[x_1,\;x_2\right]\f$.
 *
 * \return solution value. If the equation can not be solved or the solution is
 * not in the interval the value is undetermined.
 */
long double
solve_quadratic_reduced (long double a, ///< a equation coefficient.
                         long double b, ///< b equation coefficient.
                         long double x1,        ///< lower solution limit.
                         long double x2)        ///< higher solution limit.
{
  long double a2, k, x;
#if DEBUG_EQUATION
  fprintf (stderr, "solve_quadratic_reduced: start\n");
  fprintf (stderr, "solve_quadratic_reduced: a=%Lg b=%Lg\n", a, b);
  fprintf (stderr, "solve_quadratic_reduced: x1=%Lg x2=%Lg\n", x1, x2);
#endif
  a2 = -0.5L * a;
  k = sqrtl (a2 * a2 - b);
  x = a2 + k;
  if (x < x1 || x > x2)
    x = a2 - k;
#if DEBUG_EQUATION
  fprintf (stderr, "solve_quadratic_reduced: x=%Lg\n", x);
  fprintf (stderr, "solve_quadratic_reduced: end\n");
#endif
  return x;
}

/**
 * Function to calculate the solution of a 2nd order equation.
 *
 * This function calculates the solution of a 2nd order equation in the
 * form:
 * \f[a\,x^2+b\,x+c=0\f]
 * in the interval \f$x\in\left[x_1,\;x_2\right]\f$.
 *
 * \return solution value. If the equation can not be solved or the solution is
 * not in the interval the value is undetermined.
 */
long double
solve_quadratic (long double a, ///< a equation coefficient.
                 long double b, ///< b equation coefficient.
                 long double c, ///< c equation coefficient.
                 long double x1,        ///< lower solution limit.
                 long double x2)        ///< higher solution limit.
{
  long double x;
#if DEBUG_EQUATION
  fprintf (stderr, "solve_quadratic: start\n");
  fprintf (stderr, "solve_quadratic: a=%Lg b=%Lg c=%Lg\n", a, b, c);
  fprintf (stderr, "solve_quadratic: x1=%Lg x2=%Lg\n", x1, x2);
#endif
  if (a == 0.L)
    x = -c / b;
  else
    x = solve_quadratic_reduced (b / a, c / a, x1, x2);
#if DEBUG_EQUATION
  fprintf (stderr, "solve_quadratic: x=%Lg\n", x);
  fprintf (stderr, "solve_quadratic: end\n");
#endif
  return x;
}

/**
 * Function to calculate the solution of a reduced 3rd order equation.
 *
 * This function calculates the solution of a reduced 3rd order equation in the
 * form:
 * \f[x^3+a\,x^2+b\,x+c=0\f]
 * in the interval \f$x\in\left[x_1,\;x_2\right]\f$.
 *
 * \return solution value. If the equation can not be solved or the solution is
 * not in the interval the value is undetermined.
 */
long double
solve_cubic_reduced (long double a,     ///< a equation coefficient.
                     long double b,     ///< b equation coefficient.
                     long double c,     ///< c equation coefficient.
                     long double x1,    ///< lower solution limit.
                     long double x2)    ///< higher solution limit.
{
  long double k0, k1, k2;
#if DEBUG_EQUATION
  fprintf (stderr, "solve_cubic_reduced: start\n");
  fprintf (stderr, "solve_cubic_reduced: a=%Lg b=%Lg c=%Lg\n", a, b, c);
  fprintf (stderr, "solve_cubic_reduced: x1=%Lg x2=%Lg\n", x1, x2);
#endif
  a /= 3.L;
  k0 = a * a;
  k1 = b / 3.L - k0;
  k0 = (b * a - c) / 2.L - a * k0;
  k2 = k1 * k1 * k1 + k0 * k0;
  if (k2 < 0.L)
    {
      k1 = sqrtl (-k1);
      k0 = acosl (k0 / (k1 * k1 * k1)) / 3.L;
      k1 *= 2.L;
      k2 = k1 * cosl (k0) - a;
      if (k2 < x1 || k2 > x2)
        {
          k2 = k1 * cosl (k0 + 2.L * M_PIl / 3.L) - a;
          if (k2 < x1 || k2 > x2)
            k2 = k1 * cosl (k0 - 2.L * M_PIl / 3.L) - a;
        }
    }
  else
    {
      k1 = sqrtl (k2);
      k2 = k0 + k1;
      k2 = cbrtl (k2);
      k0 -= k1;
      k2 += cbrtl (k0);
      k2 -= a;
    }
#if DEBUG_EQUATION
  fprintf (stderr, "solve_cubic_reduced: x=%Lg\n", k2);
  fprintf (stderr, "solve_cubic_reduced: end\n");
#endif
  return k2;
}

/**
 * Function to calculate the solution of a 3rd order equation.
 *
 * This function calculates the solution of a 3rd order equation in the
 * form:
 * \f[a\,x^3+b\,x^2+c\,x+d=0\f]
 * in the interval \f$x\in\left[x_1,\;x_2\right]\f$.
 *
 * \return solution value. If the equation can not be solved or the solution is
 * not in the interval the value is undetermined.
 */
long double
solve_cubic (long double a,     ///< a equation coefficient.
             long double b,     ///< b equation coefficient.
             long double c,     ///< c equation coefficient.
             long double d,     ///< d equation coefficient.
             long double x1,    ///< lower solution limit.
             long double x2)    ///< higher solution limit.
{
  long double x;
#if DEBUG_EQUATION
  fprintf (stderr, "solve_cubic: start\n");
  fprintf (stderr, "solve_cubic: a=%Lg b=%Lg c=%Lg d=%Lg\n", a, b, c, d);
  fprintf (stderr, "solve_cubic: x1=%Lg x2=%Lg\n", x1, x2);
#endif
  if (a == 0.L)
    x = solve_quadratic (b, c, d, x1, x2);
  else
    x = solve_cubic_reduced (b / a, c / a, d / a, x1, x2);
#if DEBUG_EQUATION
  fprintf (stderr, "solve_cubic: x=%Lg\n", x);
  fprintf (stderr, "solve_cubic: end\n");
#endif
  return x;
}

/**
 * Function to get an integer number of a XML node property.
 *
 * \return Integer number value.
 */
int
xml_node_get_int (xmlNode * node,       ///< XML node.
                  const xmlChar * prop, ///< XML property.
                  int *error_code)      ///< Error code.
{
  int i = 0;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%d", &i) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return i;
}

/**
 * Function to get an unsigned integer number of a XML node property.
 *
 * \return Unsigned integer number value.
 */
unsigned int
xml_node_get_uint (xmlNode * node,      ///< XML node.
                   const xmlChar * prop,        ///< XML property.
                   int *error_code)     ///< Error code.
{
  unsigned int i = 0;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%u", &i) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return i;
}

/**
 * Function to get an unsigned integer number of a XML node property with a
 *   default value.
 *
 * \return Unsigned integer number value.
 */
unsigned int
xml_node_get_uint_with_default (xmlNode * node, ///< XML node.
                                const xmlChar * prop,   ///< XML property.
                                unsigned int default_value,
                                ///< default value.
                                int *error_code)        ///< Error code.
{
  unsigned int i;
  if (xmlHasProp (node, prop))
    i = xml_node_get_uint (node, prop, error_code);
  else
    {
      i = default_value;
      *error_code = 0;
    }
  return i;
}

/**
 * Function to get a floating point number of a XML node property.
 *
 * \return Floating point number value.
 */
long double
xml_node_get_float (xmlNode * node,     ///< XML node.
                    const xmlChar * prop,       ///< XML property.
                    int *error_code)    ///< Error code.
{
  long double x = 0.L;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%Lf", &x) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return x;
}

/**
 * Function to get a floating point number of a XML node property with a 
 *   default value.
 *
 * \return Floating point number value.
 */
long double
xml_node_get_float_with_default (xmlNode * node,        ///< XML node.
                                 const xmlChar * prop,  ///< XML property.
                                 long double default_value,
                                 ///< default value.
                                 int *error_code)       ///< Error code.
{
  long double x;
  if (xmlHasProp (node, prop))
    x = xml_node_get_float (node, prop, error_code);
  else
    {
      x = default_value;
      *error_code = 0;
    }
  return x;
}

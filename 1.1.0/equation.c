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
 * \file equation.c
 * \brief Source file with the equation data and functions.
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
#include "utils.h"
#include "equation.h"

#define DEBUG_EQUATION 0        ///< macro to debug the equation functions.

long double r0[3];              ///< position vector.
long double r1[3];              ///< velocity vector.
long double r2[3];              ///< acceleration vector.
long double ro0[3];             ///< backup of the position vector.
long double ro1[3];             ///< backup of the velocity vector.
long double ro2[3];             ///< backup of the acceleration vector.
void (*equation_acceleration) (Equation * eq, long double *r0,
                               long double *r1, long double *r2, long double t);
///< pointer to the function to calculate the acceleration.
void (*equation_solution) (Equation * eq, long double *r0, long double *r1,
                           long double t);
///< pointer to the function to calculate the analytical solution.
long double (*equation_step_size) (Equation * eq);
///< pointer to the function to calculate the time step size.
int (*equation_land) (Equation * eq, long double, long double *t,
                      long double *dt);
///< pointer to the function to finalize the trajectory.
long double kt;
///< stability time step size coefficient.
long double dt;
///< time step size.
unsigned long int nevaluations;
///< number of evaluations of the acceleration function.

/**
 * Function to calculate the acceleration on non-resitance model.
 *
 * This function calculates the acceleration vector on a non-resistance
 * model. The movement equation is:
 * \f{equation}\ddot{\vec{r}}=\vec{g}\f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector.
 */
static void
equation_acceleration_0 (Equation * eq __attribute__ ((unused)),
                         ///< Equation struct.
                         long double *r0 __attribute__ ((unused)),
                         ///< position vector.
                         long double *r1 __attribute__ ((unused)),
                         ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t __attribute__ ((unused)))
  ///< actual time.
{
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_0: start\n");
#endif
  r2[0] = r2[1] = 0.L;
  r2[2] = -G;
  ++nevaluations;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_0: ax=%Lg ay=%Lg az=%Lg\n",
           r2[0], r2[1], r2[2]);
  fprintf (stderr, "equation_acceleration_0: end\n");
#endif
}

/**
 * Function to solve the non-resistance model.
 *
 * This function solves the movement on a resistance model characterized by the
 * movement equation:
 * \f{equation}{\ddot{\vec{r}}=\vec{g},\f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector.
 * The analytical solution is:
 * \f{equation}\dot{\vec{r}}=\dot{\vec{r}}_0+\vec{g}\,t,\f}
 * \f{equation}\vec{r}=\vec{r}_0+\dot{\vec{r}}_0\,t+\frac12\,\vec{g}\,t^2.\f}
 */
static void
equation_solution_0 (Equation * eq,     ///< Equation struct.
                     long double *r0,   ///< position vector.
                     long double *r1,   ///< velocity vector.
                     long double t)     ///< time.
{
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_0: start\n");
#endif
  r1[0] = eq->v[0];
  r1[1] = eq->v[1];
  r1[2] = eq->v[2] - eq->g * t;
  r0[0] = eq->r[0] + eq->v[0] * t;
  r0[1] = eq->r[1] + eq->v[1] * t;
  r0[2] = eq->r[2] + t * (eq->v[2] - t * 0.5L * G);
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_0: vx=%Lg vy=%Lg vz=%Lg\n",
           r1[0], r1[1], r1[2]);
  fprintf (stderr, "equation_solution_0: x=%Lg y=%Lg z=%Lg\n",
           r0[0], r0[1], r0[2]);
  fprintf (stderr, "equation_solution_0: end\n");
#endif
}

/**
 * Function to calculate the acceleration on the 1st resistance model.
 *
 * This function calculates the acceleration vector on a resistance model
 * model characterized by the movement equation:
 * \f[\ddot{\vec{r}}=\vec{g}-\lambda\,\left(\dot{\vec{r}}-\vec{w}\right)\f]
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector,
 * \f$\vec{w}=\left(w_x,\;w_y\;0\right)\f$ the wind velocity vector and
 * \f$\lambda\f$ a resistance coefficient.
 */
static void
equation_acceleration_1 (Equation * eq, ///< Equation struct.
                         long double *r0 __attribute__ ((unused)),
                         ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t __attribute__ ((unused)))
  ///< actual time.
{
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_1: start\n");
#endif
  r2[0] = -eq->lambda * (r1[0] - eq->w[0]);
  r2[1] = -eq->lambda * (r1[1] - eq->w[1]);
  r2[2] = -eq->g - eq->lambda * r1[2];
  ++nevaluations;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_1: ax=%Lg ay=%Lg az=%Lg\n",
           r2[0], r2[1], r2[2]);
  fprintf (stderr, "equation_acceleration_1: end\n");
#endif
}

/**
 * Function to solve the 1st resistance model.
 *
 * This function solves the movement on a resistance model characterized by the
 * movement equation:
 * \f{equation}
 * \ddot{\vec{r}}=\vec{g}-\lambda\,\left(\dot{\vec{r}}-\vec{w}\right),
 * \f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector,
 * \f$\vec{w}=(w_x,\;w_y,\;0)\f$ the wind velocity vector and
 * \f$\lambda\f$ a resistance coefficient.
 * The analytical solution is:
 * \f{equation}
 * \dot{\vec{r}}=\dot{\vec{r}}_0\,\exp\left(-\lambda\,t\right)
 * +\left(\vec{w}+\frac{\vec{g}}{\lambda}\right)
 * \,\left[1-\exp\left(-\lambda\,t\right)\right],
 * \f}
 * \f{equation}
 * \vec{r}=\vec{r}_0+\left(\vec{w}+\frac{\vec{g}}{\lambda}\right)\,t
 * +\frac{\dot{\vec{r}}_0-\vec{w}-\vec{g}/\lambda}{\lambda}
 * \,\left[1-\exp\left(-\lambda\,t\right)\right].
 * \f}
 */
static void
equation_solution_1 (Equation * eq,     ///< Equation struct.
                     long double *r0,   ///< position vector.
                     long double *r1,   ///< velocity vector.
                     long double t)     ///< time.
{
  long double v[2];
  long double li, gl, elt, k;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_1: start\n");
#endif
  v[0] = eq->v[0] - eq->w[0];
  v[1] = eq->v[1] - eq->w[1];
  elt = expl (-eq->lambda * t);
  r1[0] = eq->w[0] + v[0] * elt;
  r1[1] = eq->w[1] + v[1] * elt;
  li = 1.L / eq->lambda;
  gl = eq->g * li;
  r1[2] = (eq->v[2] + gl) * elt - gl;
  k = li * (1.L - elt);
  r0[0] = eq->r[0] + eq->w[0] * t + v[0] * k;
  r0[1] = eq->r[1] + eq->w[1] * t + v[1] * k;
  r0[2] = eq->r[2] - gl * t + (eq->v[2] + gl) * k;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_1: vx=%Lg vy=%Lg vz=%Lg\n",
           r1[0], r1[1], r1[2]);
  fprintf (stderr, "equation_solution_1: x=%Lg y=%Lg z=%Lg\n",
           r0[0], r0[1], r0[2]);
  fprintf (stderr, "equation_solution_1: end\n");
#endif
}

/**
 * Function to calculate the acceleration on 2nd resistance model.
 *
 * This function calculates the acceleration vector on a resistance model
 * model characterized by the movement equations:
 * \f[\left.\begin{array}{r}
 * \ddot{x}=-\lambda\,\left|\dot{x}-w_x\right|\,\left(\dot{x}-w_x\right),\\
 * \ddot{y}=-\lambda\,\left|\dot{y}-w_y\right|\,\left(\dot{y}-w_y\right),\\
 * \ddot{z}=-g-\lambda\,\left|\dot{z}\right|\,\dot{z},
 * \end{array}\right\}\f]
 * with g the gravitational constant, \f$w_x\f$ and \f$w_y\f$ the wind velocity
 * vector components and \f$\lambda\f$ a resistance coefficient.
 */
static void
equation_acceleration_2 (Equation * eq, ///< Equation struct.
                         long double *r0 __attribute__ ((unused)),
                         ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t __attribute__ ((unused)))
  ///< actual time.
{
  long double v[2];
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_2: start\n");
#endif
  v[0] = r1[0] - eq->w[0];
  v[1] = r1[1] - eq->w[1];
  r2[0] = -eq->lambda * fabsl (v[0]) * v[0];
  r2[1] = -eq->lambda * fabsl (v[1]) * v[1];
  r2[2] = -eq->g - eq->lambda * fabsl (r1[2]) * r1[2];
  ++nevaluations;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_2: ax=%Lg ay=%Lg az=%Lg\n",
           r2[0], r2[1], r2[2]);
  fprintf (stderr, "equation_acceleration_2: end\n");
#endif
}

/**
 * Function to solve the 2nd resistance mode.
 *
 * This function solves the movement on a resistance model characterized by the
 * movement equations:
 * \f{equation}\left.\begin{array}{r}
 * \ddot{x}=-\lambda\,\left|\dot{x}-w_x\right|\,\left(\dot{x}-w_x\right)\\
 * \ddot{y}=-\lambda\,\left|\dot{y}-w_y\right|\,\left(\dot{y}-w_y\right)\\
 * \ddot{z}=-g-\lambda\,\left|\dot{z}\right|\,\dot{z}
 * \end{array}\right\}\f}
 * with g the gravitational constant, \f$w_x\f$ and \f$w_y\f$ the wind velocity
 * vector components and \f$\lambda\f$ a resistance coefficient.
 * The analytical solution is:
 * \f{equation}
 * \dot{x}=w_x+\frac{\dot{x}_0-w_x}{1+\lambda\,\left|\dot{x}_0-w_x\right|\,t},
 * \f}
 * \f{equation}
 * \dot{y}=w_y+\frac{\dot{y}_0-w_y}{1+\lambda\,\left|\dot{y}_0-w_y\right|\,t},
 * \f}
 * \f{equation}
 * \dot{z}=\left\{\begin{array}{cl}
 * \dot{z}_0\leq 0\Rightarrow & \sqrt{\frac{g}{\lambda}}\,
 * \frac{\dot{z}_0\,\cosh\left(\sqrt{g\,\lambda}\,t\right)
 * -\sqrt{g/\lambda}\,\sinh\left(\sqrt{g\,\lambda}\,t\right)}
 * {\sqrt{g/\lambda}\,\cosh\left(\sqrt{g\,\lambda}\,t\right)
 * -\dot{z}_0\,\sinh\left(\sqrt{g\,\lambda}\,t\right)},\\
 * \dot{z}_0>0,\;
 * t\leq\frac{\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * \sqrt{\frac{g}{\lambda}}
 * \,\tan\left[\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)
 * -\sqrt{g\,\lambda}\,t\right],\\
 * \dot{z}_0>0,\;
 * t>\frac{\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * -\sqrt{\frac{g}{\lambda}}\,\tanh\left[\sqrt{g\,\lambda}\,t
 *  -\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)\right],
 * \end{array}\right.
 * \f}
 * \f{equation}
 * x=x_0+w_x\,t+\frac{\dot{x}_0-w_x}{\lambda\,\left|\dot{x}_0-w_x\right|}
 * \,\ln\left(1+\lambda\,\left|\dot{x}_0-w_x\right|\,t\right),
 * \f}
 * \f{equation}
 * y=y_0+w_y\,t+\frac{\dot{y}_0-w_y}{\lambda\,\left|\dot{y}_0-w_y\right|}
 * \,\ln\left(1+\lambda\,\left|\dot{y}_0-w_y\right|\,t\right),
 * \f}
 * \f{equation}
 * z=\left\{\begin{array}{cl}
 * \dot{z}_0\leq 0\Rightarrow & z_0
 * -\frac{\ln\left[\cosh\left(\sqrt{g\,\lambda}\,t\right)
 * -\dot{z}_0\,\sqrt{\lambda/g}\,
 * \sinh\left(\sqrt{g\,\lambda}\,t\right)\right]}{\lambda},\\
 * \dot{z}_0>0,\;
 * t\leq\frac{\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * z_0+\frac{1}{\lambda}\,\ln\left\{\frac{\cos\left[
 * \arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)-\sqrt{g\,\lambda}\,t\right]}
 * {\cos\left[\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)\right]}\right\},\\
 * \dot{z}_0>0,\;
 * t>\frac{\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * z_0-\frac{\ln\left\{
 * \cos\left[\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)\right]
 * \,\cosh\left[\sqrt{g\,\lambda}\,t
 * -\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)\right]\right\}}{\lambda}
 * \end{array}\right.
 * \f}
 */
static void
equation_solution_2 (Equation * eq,     ///< Equation struct.
                     long double *r0,   ///< position vector.
                     long double *r1,   ///< velocity vector.
                     long double t)     ///< actual time.
{
  long double v[2], k[2];
  long double tc, g_l, gl, glt, lt, li, alpha;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_2: start\n");
#endif
  v[0] = eq->v[0] - eq->w[0];
  v[1] = eq->v[1] - eq->w[1];
  lt = eq->lambda * t;
  k[0] = 1.L + lt * fabsl (v[0]);
  k[1] = 1.L + lt * fabsl (v[1]);
  r1[0] = eq->w[0] + v[0] / k[0];
  r1[1] = eq->w[1] + v[1] / k[1];
  li = 1.L / eq->lambda;
  k[0] = li * logl (k[0]);
  k[1] = li * logl (k[1]);
  r0[0] = eq->r[0] + eq->w[0] * t;
  r0[1] = eq->r[1] + eq->w[1] * t;
  if (v[0] >= 0.L)
    r0[0] += k[0];
  else
    r0[0] -= k[0];
  if (v[1] >= 0.L)
    r0[1] += k[1];
  else
    r0[1] -= k[1];
  gl = sqrtl (eq->g * eq->lambda);
  glt = gl * t;
  g_l = sqrtl (eq->g / eq->lambda);
  r0[2] = eq->r[2];
  if (eq->v[2] <= 0.L)
    {
      k[0] = coshl (glt);
      k[1] = sinhl (glt);
      r1[2] = g_l * (eq->v[2] * k[0] - g_l * k[1])
        / (g_l * k[0] - eq->v[2] * k[1]);
      r0[2] -= li * logl (k[0] - eq->v[2] * k[1] / g_l);
    }
  else
    {
      alpha = atanl (eq->v[2] / g_l);
      tc = alpha / gl;
      if (t <= tc)
        {
          r1[2] = g_l * tanl (alpha - glt);
          r0[2] += li * logl (cosl (alpha - glt) / cosl (alpha));
        }
      else
        {
          t -= tc;
          glt = gl * t;
          r1[2] = -g_l * tanhl (glt);
          r0[2] -= li * logl (cosl (alpha) * coshl (glt));
        }
    }
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_2: vx=%Lg vy=%Lg vz=%Lg\n",
           r1[0], r1[1], r1[2]);
  fprintf (stderr, "equation_solution_2: x=%Lg y=%Lg z=%Lg\n",
           r0[0], r0[1], r0[2]);
  fprintf (stderr, "equation_solution_2: end\n");
#endif
}

/**
 * Function to calculate the acceleration on a forced model.
 *
 * This function calculates the acceleration vector on a forced model. The
 * movement equation is:
 * \f{equation}
 * \ddot{\vec{r}}=\vec{g}+\vec{w}\,\exp\left(-\lambda\,t\right)
 * \f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector, \f$\lambda\f$ the
 * force decay factor and \f$\vec{w}=\left(w_x,\;w_y\;0\right)\f$ the force
 * vector.
 */
static void
equation_acceleration_3 (Equation * eq, ///< Equation struct.
                         long double *r0 __attribute__ ((unused)),
                         ///< position vector.
                         long double *r1 __attribute__ ((unused)),
                         ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t) ///< actual time.
{
  long double elt;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_3: start\n");
#endif
  elt = expl (-eq->lambda * t);
  r2[0] = eq->w[0] * elt;
  r2[1] = eq->w[1] * elt;
  r2[2] = -G;
  ++nevaluations;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_acceleration_3: ax=%Lg ay=%Lg az=%Lg\n",
           r2[0], r2[1], r2[2]);
  fprintf (stderr, "equation_acceleration_3: end\n");
#endif
}

/**
 * Function to solve the forced model.
 *
 * This function solves the movement on a forced model characterized by the
 * movement equation:
 * \f{equation}
 * \ddot{\vec{r}}=\vec{g}+\vec{w}\,\exp\left(-\lambda\,t\right)
 * \f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector, \f$\lambda\f$ the
 * force decay factor and \f$\vec{w}=\left(w_x,\;w_y\;0\right)\f$ the force
 * vector.
 * The analytical solution is:
 * \f{equation}
 * \dot{\vec{r}}=\dot{\vec{r}}_0+\vec{g}\,t
 * +\frac{\vec{w}}{\lambda}\left[1-\exp\left(-\lambda\,t\right)\right],\f}
 * \f{equation}\vec{r}=\vec{r}_0
 * +\left(\dot{\vec{r}}_0+\frac{\vec{w}}{\lambda}\right)\,t
 * +\frac12\,\vec{g}\,t^2
 * +\frac{\vec{w}}{\lambda^2}\,\left[1-\exp\left(-\lambda\,t\right)\right].
 * \f}
 */
static void
equation_solution_3 (Equation * eq,     ///< Equation struct.
                     long double *r0,   ///< position vector.
                     long double *r1,   ///< velocity vector.
                     long double t)     ///< time.
{
  long double li, k;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_3: start\n");
#endif
  li = 1.L / eq->lambda;
  k = li * (1.L - expl (-eq->lambda * t));
  r1[0] = eq->v[0] + eq->w[0] * k;
  r1[1] = eq->v[1] + eq->w[1] * k;
  r1[2] = eq->v[2] - eq->g * t;
  k *= li;
  r0[0] = eq->r[0] + (eq->v[0] + eq->w[0] * li) * t - eq->w[0] * k;
  r0[1] = eq->r[1] + (eq->v[1] + eq->w[1] * li) * t - eq->w[1] * k;
  r0[2] = eq->r[2] + t * (eq->v[2] - t * 0.5L * G);
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solution_3: vx=%Lg vy=%Lg vz=%Lg\n",
           r1[0], r1[1], r1[2]);
  fprintf (stderr, "equation_solution_3: x=%Lg y=%Lg z=%Lg\n",
           r0[0], r0[1], r0[2]);
  fprintf (stderr, "equation_solution_3: end\n");
#endif
}

/**
 * Function to solve numerically the equation by the mean point method.
 *
 * \return solution time.
 */
long double
equation_solve (Equation * eq,  ///< Equation struct.
                long double *r0,        ///< position vector solution.
                long double *r1)        ///< velocity vector solution.
{
  long double r02[3], r12[3];
  long double t1, t2, t3;
  unsigned int i;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solve: start\n");
#endif
  t1 = 0.L;
  t2 = 1.L;
  equation_solution (eq, r02, r12, t2);
  while (r02[2] > 0.L)
    {
      t2 *= 2.L;
      equation_solution (eq, r02, r12, t2);
    }
  for (i = 0; i < 64; ++i)
    {
      t3 = 0.5L * (t1 + t2);
      equation_solution (eq, r02, r12, t3);
      if (r02[2] > 0.L)
        t1 = t3;
      else
        t2 = t3;
    }
  memcpy (r0, r02, 3 * sizeof (long double));
  memcpy (r1, r12, 3 * sizeof (long double));
#if DEBUG_EQUATION
  fprintf (stderr, "equation_solve: vx=%Lg vy=%Lg vz=%Lg\n",
           r1[0], r1[1], r1[2]);
  fprintf (stderr, "equation_solve: x=%Lg y=%Lg z=%Lg\n", r0[0], r0[1], r0[2]);
  fprintf (stderr, "equation_solve: t=%Lg\n", t3);
  fprintf (stderr, "equation_solve: end\n");
#endif
  return t3;
}

/**
 * Function to set a constant time step size.
 *
 * \return time step size.
 */
static long double
equation_step_size_0 (Equation * eq __attribute__ ((unused)))
///< Equation struct.
{
#if DEBUG_EQUATION
  fprintf (stderr, "equation_step_size_0: start\n");
  fprintf (stderr, "equation_step_size_0: dt=%Lg\n", dt);
  fprintf (stderr, "equation_step_size_0: end\n");
#endif
  return dt;
}

/**
 * Function to set the time step size based on stability condition for the 1st
 * resistance model.
 *
 * \return time step size.
 */
static long double
equation_step_size_1 (Equation * eq)    ///< Equation struct.
{
  long double dt;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_step_size_1: start\n");
#endif
  dt = kt / fabsl (eq->lambda);
#if DEBUG_EQUATION
  fprintf (stderr, "equation_step_size_1: dt=%Lg\n", dt);
  fprintf (stderr, "equation_step_size_1: end\n");
#endif
  return dt;
}

/**
 * Function to set the time step size based on stability condition for the 2nd
 * resistance model.
 *
 * \return time step size.
 */
static long double
equation_step_size_2 (Equation * eq)    ///< Equation struct.
{
  long double dt;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_step_size_2: start\n");
#endif
  dt = kt / (fabsl (eq->lambda) *
             fmaxl (fabsl (r1[0] - eq->w[0]),
                    fmaxl (fabsl (r1[1] - eq->w[1]), fabsl (r1[2]))));
#if DEBUG_EQUATION
  fprintf (stderr, "equation_step_size_2: dt=%Lg\n", dt);
  fprintf (stderr, "equation_step_size_2: end\n");
#endif
  return dt;
}

/**
 * Function to finish the trajectory based on final time.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_0 (Equation * eq, ///< Equation struct.
                 long double to,        ///< old time.
                 long double *t,        ///< next time.
                 long double *dt)       ///< time step size.
{
  long double tf;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_0: start\n");
  fprintf (stderr, "equation_land_0: to=%Lg\n", to);
#endif
  tf = eq->tf;
  if (to >= tf)
    {
#if DEBUG_EQUATION
      fprintf (stderr, "equation_land_0: landing\n");
      fprintf (stderr, "equation_land_0: end\n");
#endif
      return 1;
    }
  *t = to + *dt;
  if (*t >= tf)
    {
      *dt = tf - to;
      *t = tf;
    }
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_0: t=%Lg dt=%Lg\n", *t, *dt);
  fprintf (stderr, "equation_land_0: no landing\n");
  fprintf (stderr, "equation_land_0: end\n");
#endif
  return 0;
}

/**
 * Function to finish the trajectory based on 1st order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_1 (Equation * eq __attribute__ ((unused)),
                 ///< Equation struct.
                 long double to,        ///< old time.
                 long double *t,        ///< next time.
                 long double *dt)       ///< time step size.
{
  long double h;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_1: start\n");
  fprintf (stderr, "equation_land_1: to=%Lg\n", to);
#endif
  if (r0[2] > 0.)
    {
      *t = to + *dt;
#if DEBUG_EQUATION
      fprintf (stderr, "equation_land_1: t=%Lg dt=%Lg\n", *t, *dt);
      fprintf (stderr, "equation_land_1: no landing\n");
      fprintf (stderr, "equation_land_1: end\n");
#endif
      return 0;
    }
  h = r0[2] / r1[2];
  r0[0] -= h * r1[0];
  r0[1] -= h * r1[1];
  r0[2] -= h * r1[2];
  r1[0] -= h * r2[0];
  r1[1] -= h * r2[1];
  r1[2] -= h * r2[2];
  *t = to - h;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_1: t=%Lg dt=%Lg\n", *t, *dt);
  fprintf (stderr, "equation_land_1: landing\n");
  fprintf (stderr, "equation_land_1: end\n");
#endif
  return 1;
}

/**
 * Function to finish the trajectory based on 2nd order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_2 (Equation * eq __attribute__ ((unused)),
                 ///< Equation struct.
                 long double to,        ///< old time.
                 long double *t,        ///< next time.
                 long double *dt)       ///< time step size.
{
  long double h;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_2: start\n");
  fprintf (stderr, "equation_land_2: to=%Lg\n", to);
#endif
  if (r0[2] > 0.)
    {
      *t = to + *dt;
#if DEBUG_EQUATION
      fprintf (stderr, "equation_land_2: t=%Lg dt=%Lg\n", *t, *dt);
      fprintf (stderr, "equation_land_2: no landing\n");
      fprintf (stderr, "equation_land_2: end\n");
#endif
      return 0;
    }
  h = solve_quadratic (0.5L * r2[2], -r1[2], r0[2], 0.L, *dt);
  r0[0] -= h * (r1[0] - h * 0.5L * r2[0]);
  r0[1] -= h * (r1[1] - h * 0.5L * r2[1]);
  r0[2] -= h * (r1[2] - h * 0.5L * r2[2]);
  r1[0] -= h * r2[0];
  r1[1] -= h * r2[1];
  r1[2] -= h * r2[2];
  *t = to - h;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_2: t=%Lg dt=%Lg\n", *t, *dt);
  fprintf (stderr, "equation_land_2: landing\n");
  fprintf (stderr, "equation_land_2: end\n");
#endif
  return 1;
}

/**
 * Function to finish the trajectory based on 3rd order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_3 (Equation * eq __attribute__ ((unused)),
                 ///< Equation struct.
                 long double to,        ///< old time.
                 long double *t,        ///< next time.
                 long double *dt)       ///< time step size.
{
  long double h, r3[3];
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_3: start\n");
  fprintf (stderr, "equation_land_3: to=%Lg\n", to);
#endif
  if (r0[2] > 0.)
    {
      *t = to + *dt;
#if DEBUG_EQUATION
      fprintf (stderr, "equation_land_3: t=%Lg dt=%Lg\n", *t, *dt);
      fprintf (stderr, "equation_land_3: no landing\n");
      fprintf (stderr, "equation_land_3: end\n");
#endif
      return 0;
    }
  r3[0] = (r2[0] - ro2[0]) / *dt;
  r3[1] = (r2[1] - ro2[1]) / *dt;
  r3[2] = (r2[2] - ro2[2]) / *dt;
  h = solve_cubic (-1.L / 6.L * r3[2], 0.5L * r2[2], -r1[2], r0[2], 0.L, *dt);
  r0[0] -= h * (r1[0] - h * (0.5L * r2[0] - h * 1.L / 6.L * r3[0]));
  r0[1] -= h * (r1[1] - h * (0.5L * r2[1] - h * 1.L / 6.L * r3[1]));
  r0[2] -= h * (r1[2] - h * (0.5L * r2[2] - h * 1.L / 6.L * r3[2]));
  r1[0] -= h * (r2[0] - h * 0.5L * r3[0]);
  r1[1] -= h * (r2[1] - h * 0.5L * r3[1]);
  r1[2] -= h * (r2[2] - h * 0.5L * r3[2]);
  *t = to - h;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_land_3: t=%Lg dt=%Lg\n", *t, *dt);
  fprintf (stderr, "equation_land_3: landing\n");
  fprintf (stderr, "equation_land_3: end\n");
#endif
  return 1;
}

/**
 * Function to init the equation variables.
 */
void
equation_init (Equation * eq,   ///< Equation struct.
               gsl_rng * rng)   ///< gsl_rng struct.
{
  long double v, ha, va;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_init: start\n");
#endif
  switch (eq->type)
    {
    case 1:
    case 2:
    case 3:
      eq->lambda = eq->min_lambda
        + (eq->max_lambda - eq->min_lambda) * gsl_rng_uniform (rng);
    }
  v = eq->min_velocity
    + (eq->max_velocity - eq->min_velocity) * gsl_rng_uniform (rng);
  ha = 2.L * M_PIl * gsl_rng_uniform (rng);
  eq->r[0] = eq->r[1] = 0.L;
  va = eq->vertical_angle * M_PIl / 180.L;
  eq->v[2] = v * sinl (va);
  eq->v[1] = v * cosl (va);
  eq->v[0] = eq->v[1] * cosl (ha);
  eq->v[1] *= sinl (ha);
  v = eq->max_wind * gsl_rng_uniform (rng);
  ha = 2.L * M_PIl * gsl_rng_uniform (rng);
  eq->w[0] = v * cosl (ha);
  eq->w[1] = v * sinl (ha);
#if DEBUG_EQUATION
  fprintf (stderr, "equation_init: end\n");
#endif
}

/**
 * Function to read the equation data on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
equation_read_xml (Equation * eq,       ///< Equation struct.
                   xmlNode * node,      ///< XML node.,
									 unsigned int initial)   ///< type of initial conditions.
{
  const char *message[] = {
    "Bad XML node",
    "Bad type",
    "Unknown type",
    "Bad x",
    "Bad y",
    "Bad z",
    "Bad vx",
    "Bad vy",
    "Bad vz",
    "Bad wx",
    "Bad wy",
    "Bad lambda",
		"Bad minimum velocity",
		"Bad maximum velocity",
		"Bad vertical angle",
		"Bad maximum wind",
		"Bad minimum lambda",
		"Bad maximum lambda",
    "Bad g",
    "Bad time step type",
    "Bad dt",
    "Bad kt",
    "Unknown time step type",
    "Bad land type",
    "Bad t",
    "Unknown land type"
  };
  int e, error_code;
#if DEBUG_EQUATION
  fprintf (stderr, "equation_read_xml: start\n");
  fprintf (stderr, "equation_read_xml: name=%s\n", node->name);
#endif
  if (xmlStrcmp (node->name, XML_EQUATION))
    {
      e = 0;
      goto exit_on_error;
    }
  eq->type = xml_node_get_uint (node, XML_TYPE, &error_code);
  if (error_code)
    {
      e = 1;
      goto exit_on_error;
    }
  switch (eq->type)
    {
    case 0:
      equation_acceleration = equation_acceleration_0;
      equation_solution = equation_solution_0;
      break;
    case 1:
      equation_acceleration = equation_acceleration_1;
      equation_solution = equation_solution_1;
      break;
    case 2:
      equation_acceleration = equation_acceleration_2;
      equation_solution = equation_solution_2;
      break;
    case 3:
      equation_acceleration = equation_acceleration_3;
      equation_solution = equation_solution_3;
      break;
    default:
      e = 2;
      goto exit_on_error;
    }
  eq->r[0] = xml_node_get_float_with_default (node, XML_X, 0.L, &error_code);
  if (error_code)
    {
      e = 3;
      goto exit_on_error;
    }
  eq->r[1] = xml_node_get_float_with_default (node, XML_Y, 0.L, &error_code);
  if (error_code)
    {
      e = 4;
      goto exit_on_error;
    }
  eq->r[2] = xml_node_get_float (node, XML_Z, &error_code);
  if (error_code || eq->r[2] < 0.L)
    {
      e = 5;
      goto exit_on_error;
    }
	if (initial)
	  {
      eq->v[0] 
				= xml_node_get_float_with_default (node, XML_VX, 0.L, &error_code);
      if (error_code)
        {
          e = 6;
          goto exit_on_error;
        }
      eq->v[1] 
				= xml_node_get_float_with_default (node, XML_VY, 0.L, &error_code);
      if (error_code)
        {
          e = 7;
          goto exit_on_error;
        }
      eq->v[2] 
				= xml_node_get_float_with_default (node, XML_VZ, 0.L, &error_code);
      if (error_code)
        {
          e = 8;
          goto exit_on_error;
        }
      eq->w[0] 
				= xml_node_get_float_with_default (node, XML_WX, 0.L, &error_code);
      if (error_code)
        {
          e = 9;
          goto exit_on_error;
        }
      eq->w[1]
			 	= xml_node_get_float_with_default (node, XML_WY, 0.L, &error_code);
      if (error_code)
        {
          e = 10;
          goto exit_on_error;
        }
      switch (eq->type)
        {
        case 1:
        case 2:
        case 3:
          eq->lambda = xml_node_get_float_with_default (node, XML_LAMBDA, 0.L, 
							                                          &error_code);
          if (error_code)
            {
              e = 11;
              goto exit_on_error;
            }
				}
		}
	else
	  {
      eq->min_velocity
				= xml_node_get_float_with_default (node, XML_VMIN, 0.L, &error_code);
			if (error_code || eq->min_velocity < 0.L)
			  {
					e = 12;
					goto exit_on_error;
				}
      eq->max_velocity = xml_node_get_float (node, XML_VMAX, &error_code);
			if (error_code || eq->max_velocity <= 0.L)
			  {
					e = 13;
					goto exit_on_error;
				}
      eq->vertical_angle
				= xml_node_get_float (node, XML_VERTICAL_ANGLE, &error_code);
			if (error_code)
			  {
					e = 14;
					goto exit_on_error;
				}
      eq->max_wind 
				= xml_node_get_float_with_default (node, XML_WMAX, 0.L, &error_code);
			if (error_code || eq->max_wind < 0.L)
			  {
					e = 15;
					goto exit_on_error;
				}
      switch (eq->type)
        {
        case 1:
        case 2:
        case 3:
          eq->min_lambda
						= xml_node_get_float_with_default (node, XML_LAMBDA_MIN, 0.L, 
								                               &error_code);
				  if (error_code)
					  {
						  e = 16;
						  goto exit_on_error;
						}
          eq->max_lambda
						= xml_node_get_float_with_default (node, XML_LAMBDA_MAX, 0.L, 
								                               &error_code);
				  if (error_code || eq->max_lambda < eq->min_lambda)
					  {
						  e = 17;
						  goto exit_on_error;
						}
			   }	
		}
  eq->g = xml_node_get_float_with_default (node, XML_G, G, &error_code);
  if (error_code)
    {
      e = 18;
      goto exit_on_error;
    }
  eq->size_type = xml_node_get_uint (node, XML_TIME_STEP, &error_code);
  if (error_code)
    {
      e = 19;
      goto exit_on_error;
    }
  switch (eq->size_type)
    {
    case 0:
      equation_step_size = equation_step_size_0;
      dt = xml_node_get_float (node, XML_DT, &error_code);
      if (error_code)
        {
          e = 20;
          goto exit_on_error;
        }
      break;
    case 1:
      switch (eq->type)
        {
        case 1:
          equation_step_size = equation_step_size_1;
          break;
        case 2:
          equation_step_size = equation_step_size_2;
        }
      kt = xml_node_get_float (node, XML_KT, &error_code);
      if (error_code)
        {
          e = 21;
          goto exit_on_error;
        }
      break;
    default:
      e = 22;
      goto exit_on_error;
    }
  eq->land_type = xml_node_get_uint (node, XML_LAND, &error_code);
  if (error_code)
    {
      e = 23;
      goto exit_on_error;
    }
  switch (eq->land_type)
    {
    case 0:
      equation_land = equation_land_0;
      eq->tf = xml_node_get_float (node, XML_T, &error_code);
      if (error_code || eq->tf < 0.)
        {
          e = 24;
          goto exit_on_error;
        }
      break;
    case 1:
      equation_land = equation_land_1;
      break;
    case 2:
      equation_land = equation_land_2;
      break;
    case 3:
      equation_land = equation_land_3;
      break;
    default:
      e = 25;
      goto exit_on_error;
    }
#if DEBUG_EQUATION
  fprintf (stderr, "equation_read_xml: success\n");
  fprintf (stderr, "equation_read_xml: end\n");
#endif
  return 1;
exit_on_error:
#if DEBUG_EQUATION
  fprintf (stderr, "equation_read_xml: error\n");
#endif
  error_add (message[e]);
#if DEBUG_EQUATION
  fprintf (stderr, "equation_read_xml: end\n");
#endif
  return 0;
}

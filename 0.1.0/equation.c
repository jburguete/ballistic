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
#include "equation.h"

long double r0[3];              ///< position vector.
long double r1[3];              ///< velocity vector.
long double r2[3];              ///< acceleration vector.
long double ro0[3];             ///< backup of the position vector.
long double ro1[3];             ///< backup of the velocity vector.
long double ro2[3];             ///< backup of the acceleration vector.
void (*equation_acceleration) (Equation * eq, long double *r0, long double *r1,
                               long double *r2, long double t);
///< pointer to the function to calculate the acceleration.
void (*equation_solution) (Equation * eq, long double *r0, long double *r1,
                           long double t);
///< pointer to the function to calculate the analytical solution.
long double (*equation_step_size) (Equation * eq);
///< pointer to the function to calculate the time step size.
int (*equation_land) (Equation * eq, long double *t, long double *dt);
///< pointer to the function to finalize the trajectory.
long double kt;
///< stability time step size coefficient.
long double dt;
///< time step size.
unsigned long int nevaluations;
///< number of evaluations of the acceleration function.

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
  dr[0] = r1[0] - r2[0];
  dr[1] = r1[1] - r2[1];
  dr[2] = r1[2] - r2[2];
  d = sqrtl (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
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
  a2 = -0.5L * a;
  k = sqrtl (a2 * a2 - b);
  x = a2 + k;
  if (x >= x1 && x <= x2)
    return x;
  return a2 - k;
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
  if (a == 0.L)
    return -c / b;
  return solve_quadratic_reduced (b / a, c / a, x1, x2);
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
  if (a == 0.L)
    return solve_quadratic (b, c, d, x1, x2);
  return solve_cubic_reduced (b / a, c / a, d / a, x1, x2);
}

/**
 * Function to calculate the acceleration on non-resitance model.
 *
 * This function calculates the acceleration vector on a non-resistance
 * model. The movement equation is:
 * \f{equation}\ddot{\vec{r}}=\vec{g}\f}
 * with \f$\vec{g}=(0,\;0,\;-g)\f$ the gravity field vector.
 */
static void
equation_acceleration_0 (Equation * eq, ///< Equation struct.
                         long double *r0,       ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t) ///< actual time.
{
  r2[0] = r2[1] = 0.L;
  r2[2] = -G;
  ++nevaluations;
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
  r1[0] = eq->v[0];
  r1[1] = eq->v[1];
  r1[2] = eq->v[2] - G * t;
  r0[0] = eq->r[0] + eq->v[0] * t;
  r0[1] = eq->r[1] + eq->v[1] * t;
  r0[2] = eq->r[2] + t * (eq->v[2] - t * 0.5L * G);
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
                         long double *r0,       ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t) ///< actual time.
{
  r2[0] = -eq->lambda * (r1[0] - eq->w[0]);
  r2[1] = -eq->lambda * (r1[1] - eq->w[1]);
  r2[2] = -G - eq->lambda * r1[2];
  ++nevaluations;
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
 * \f$\vec{w}=(w_x,\;w_y\;0)\f$ the wind velocity vector and
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
  v[0] = eq->v[0] - eq->w[0];
  v[1] = eq->v[1] - eq->w[1];
  r1[0] = eq->w[0] + v[0] * expl (-eq->lambda * t);
  r1[1] = eq->w[1] + v[1] * expl (-eq->lambda * t);
  r1[2] = (eq->v[2] + G / eq->lambda) * expl (-eq->lambda * t) - G / eq->lambda;
  r0[0] = eq->r[0] + eq->w[0] * t
    + v[0] / eq->lambda * (1.L - expl (-eq->lambda * t));
  r0[1] = eq->r[1] + eq->w[1] * t
    + v[1] / eq->lambda * (1.L - expl (-eq->lambda * t));
  r0[2] = eq->r[2] - G / eq->lambda * t
    + (eq->v[2] + G / eq->lambda) / eq->lambda * (1.L - expl (-eq->lambda * t));
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
                         long double *r0,       ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t) ///< actual time.
{
  long double v[2];
  v[0] = r1[0] - eq->w[0];
  v[1] = r1[1] - eq->w[1];
  r2[0] = -eq->lambda * fabsl (v[0]) * v[0];
  r2[1] = -eq->lambda * fabsl (v[1]) * v[1];
  r2[2] = -G - eq->lambda * fabsl (r1[2]) * r1[2];
  ++nevaluations;
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
 * t\leq\frac{\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * \sqrt{\frac{g}{\lambda}}
 * \,\tan\left[\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)
 * -\sqrt{g\,\lambda}\,t\right],\\
 * \dot{z}_0>0,\;
 * t>\frac{\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)}
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
 * -\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\,
 *  \sinh\left(\sqrt{g\,\lambda}\,t\right)\right]}{\lambda},\\
 * \dot{z}_0>0,\;
 * t\leq\frac{\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)}
 * {\sqrt{g\,\lambda}}\Rightarrow & 
 * z_0+\frac{\ln\left\{1-\sqrt{g\,\lambda}\,t
 * \,\sec\left[\arctan\left(\dot{z}_0\,\sqrt{\lambda/g}\right)\right]\right\}}
 * {\lambda},\\
 * \dot{z}_0>0,\;
 * t>\frac{\arctan\left(\dot{z}_0\,\sqrt{\frac{\lambda}{g}}\right)}
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
  long double r[2], v[2];
  long double tc, g_l, gl;
  v[0] = eq->v[0] - eq->w[0];
  v[1] = eq->v[1] - eq->w[1];
  r1[0] = eq->w[0] + v[0] / (1.L + eq->lambda * fabsl (eq->v[0]) * t);
  r1[1] = eq->w[1] + v[1] / (1.L + eq->lambda * fabsl (eq->v[1]) * t);
  r[0] = eq->r[0] + eq->w[0] * t;
  r[1] = eq->r[1] + eq->w[1] * t;
  if (v[0] >= 0.L)
    r0[0] = r[0]
      + 1.L / eq->lambda * logl (1.L + eq->lambda * fabsl (v[0]) * t);
  else
    r0[0] = r[0]
      - 1.L / eq->lambda * logl (1.L + eq->lambda * fabsl (v[0]) * t);
  if (v[1] >= 0.L)
    r0[1] = r[1]
      + 1.L / eq->lambda * logl (1.L + eq->lambda * fabsl (v[1]) * t);
  else
    r0[1] = r[1]
      - 1.L / eq->lambda * logl (1.L + eq->lambda * fabsl (v[1]) * t);
  gl = sqrtl (G * eq->lambda);
  g_l = sqrtl (G / eq->lambda);
  r0[2] = eq->r[2];
  if (eq->v[2] <= 0.L)
    {
      r1[2] = g_l * (eq->v[2] * coshl (gl * t) - g_l * sinhl (gl * t))
        / (g_l * coshl (gl * t) - eq->v[2] * sinhl (gl * t));
      r0[2] -=
        logl (coshl (gl * t) - eq->v[2] * sinh (gl * t) / g_l) / eq->lambda;
      return;
    }
  tc = atanl (eq->v[2] / g_l) / gl;
  if (t <= tc)
    {
      r1[2] = tanl (atanl (eq->v[2] / g_l) - gl * t) * g_l;
      r0[2] += logl (1.L - gl * t / cosl (atanl (eq->v[2] / g_l))) / eq->lambda;
      return;
    }
  t -= tc;
  r1[2] = -g_l * tanhl (gl * t);
  r0[2] -= logl (cosl (atanl (eq->v[2] / g_l)) * coshl (gl * t)) / eq->lambda;
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
                         long double *r0,       ///< position vector.
                         long double *r1,       ///< velocity vector.
                         long double *r2,       ///< acceleration vector.
                         long double t) ///< actual time.
{
  r2[0] = eq->w[0] * expl (-eq->lambda * t);
  r2[1] = eq->w[1] * expl (-eq->lambda * t);
  r2[2] = -G;
  ++nevaluations;
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
  r1[0] = eq->v[0] + eq->w[0] / eq->lambda * (1.L - expl (-eq->lambda * t));
  r1[1] = eq->v[1] + eq->w[1] / eq->lambda * (1.L - expl (-eq->lambda * t));
  r1[2] = eq->v[2] - G * t;
  r0[0] = eq->r[0] + (eq->v[0] + eq->w[0] / eq->lambda) * t
    - eq->w[0] / (eq->lambda * eq->lambda) * (1.L - expl (-eq->lambda * t));
  r0[1] = eq->r[1] + (eq->v[1] + eq->w[1] / eq->lambda) * t
    - eq->w[1] / (eq->lambda * eq->lambda) * (1.L - expl (-eq->lambda * t));
  r0[2] = eq->r[2] + t * (eq->v[2] - t * 0.5L * G);
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
  t1 = 0.L;
  t2 = 1.L;
  equation_solution (eq, r02, r12, t2);
  while (r02[2] > 0.L)
    {
      t2 *= 2.L;
      equation_solution (eq, r02, r12, t2);
    }
  for (i = 0; i < 128; ++i)
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
  return t3;
}

/**
 * Function to set a constant time step size.
 *
 * \return time step size.
 */
static long double
equation_step_size_0 (Equation * eq)    ///< Equation struct.
{
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
  return kt / fabsl (eq->lambda);
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
  return kt / (fabsl (eq->lambda) *
               fmaxl (fabsl (r1[0] - eq->w[0]),
                      fmaxl (fabsl (r1[1] - eq->w[1]), fabsl (r1[2]))));
}

/**
 * Function to finish the trajectory based on final time.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_0 (Equation * eq, ///< Equation struct.
                 long double *t,        ///< time.
                 long double *dt)       ///< time step size.
{
  long double t2, tf;
  tf = eq->tf;
  if (*t >= tf)
    return 1;
  t2 = *t + *dt;
  if (t2 >= tf)
    {
      *dt = tf - *t;
      *t = tf;
    }
  else
    *t = t2;
  return 0;
}

/**
 * Function to finish the trajectory based on 1st order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_1 (Equation * eq, ///< Equation struct.
                 long double *t,        ///< time.
                 long double *dt)       ///< time step size.
{
  long double h;
  if (r0[2] > 0.)
    {
      *t += *dt;
      return 0;
    }
  h = r0[2] / r1[2];
  *dt -= h;
  r0[0] -= h * r1[0];
  r0[1] -= h * r1[1];
  r0[2] -= h * r1[2];
  r1[0] -= h * r2[0];
  r1[1] -= h * r2[1];
  r1[2] -= h * r2[2];
  *t += *dt;
  return 1;
}

/**
 * Function to finish the trajectory based on 2nd order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_2 (Equation * eq, ///< Equation struct.
                 long double *t,        ///< time.
                 long double *dt)       ///< time step size.
{
  long double h;
  if (r0[2] > 0.)
    {
      *t += *dt;
      return 0;
    }
  h = solve_quadratic (0.5L * r2[2], -r1[2], r0[2], 0.L, *dt);
  *dt -= h;
  r0[0] -= h * (r1[0] - h * 0.5L * r2[0]);
  r0[1] -= h * (r1[1] - h * 0.5L * r2[1]);
  r0[2] -= h * (r1[2] - h * 0.5L * r2[2]);
  r1[0] -= h * r2[0];
  r1[1] -= h * r2[1];
  r1[2] -= h * r2[2];
  *t += *dt;
  return 1;
}

/**
 * Function to finish the trajectory based on 3rd order landing.
 *
 * \return 1 on finish, 0 on continuing.
 */
static int
equation_land_3 (Equation * eq, ///< Equation struct.
                 long double *t,        ///< time.
                 long double *dt)       ///< time step size.
{
  long double h, r3[3];
  if (r0[2] > 0.)
    {
      *t += *dt;
      return 0;
    }
  r3[0] = (r2[0] - ro2[0]) / *dt;
  r3[1] = (r2[1] - ro2[1]) / *dt;
  r3[2] = (r2[2] - ro2[2]) / *dt;
  h = solve_cubic (-1.L / 6.L * r3[2], 0.5L * r2[2], -r1[2], r0[2], 0.L, *dt);
  *dt -= h;
  r0[0] -= h * (r1[0] - h * (0.5L * r2[0] - h * 1.L / 6.L * r3[0]));
  r0[1] -= h * (r1[1] - h * (0.5L * r2[1] - h * 1.L / 6.L * r3[1]));
  r0[2] -= h * (r1[2] - h * (0.5L * r2[2] - h * 1.L / 6.L * r3[2]));
  r1[0] -= h * (r2[0] - h * 0.5L * r3[0]);
  r1[1] -= h * (r2[1] - h * 0.5L * r3[1]);
  r1[2] -= h * (r2[2] - h * 0.5L * r3[2]);
  *t += *dt;
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
}

/**
 * Function to read the equation data on a file.
 *
 * \return 1 on success, 0 on error.
 */
int
equation_read (Equation * eq,   ///< Equation struct.
               FILE * file)     ///< file.
{
  if (fscanf (file, "%*s%*s%u", &eq->type) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%Lf", &eq->r[2]) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%Lf", &eq->min_velocity) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%Lf", &eq->max_velocity) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%Lf", &eq->vertical_angle) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%Lf", &eq->max_wind) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%u", &eq->size_type) != 1)
    return 0;
  if (fscanf (file, "%*s%*s%u", &eq->land_type) != 1)
    return 0;
  switch (eq->type)
    {
    case 0:
      equation_acceleration = equation_acceleration_0;
      equation_solution = equation_solution_0;
      break;
    case 1:
      if (fscanf (file, "%*s%*s%Lf", &eq->min_lambda) != 1)
        return 0;
      if (fscanf (file, "%*s%*s%Lf", &eq->max_lambda) != 1)
        return 0;
      equation_acceleration = equation_acceleration_1;
      equation_solution = equation_solution_1;
      break;
    case 2:
      if (fscanf (file, "%*s%*s%Lf", &eq->min_lambda) != 1)
        return 0;
      if (fscanf (file, "%*s%*s%Lf", &eq->max_lambda) != 1)
        return 0;
      equation_acceleration = equation_acceleration_2;
      equation_solution = equation_solution_2;
      break;
    case 3:
      if (fscanf (file, "%*s%*s%Lf", &eq->min_lambda) != 1)
        return 0;
      if (fscanf (file, "%*s%*s%Lf", &eq->max_lambda) != 1)
        return 0;
      equation_acceleration = equation_acceleration_3;
      equation_solution = equation_solution_3;
      break;
    default:
      return 0;
    }
  switch (eq->size_type)
    {
    case 0:
      equation_step_size = equation_step_size_0;
      if (fscanf (file, "%*s%*s%Lf", &dt) != 1)
        return 0;
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
      if (fscanf (file, "%*s%*s%Lf", &kt) != 1)
        return 0;
      break;
    default:
      return 0;
    }
  switch (eq->land_type)
    {
    case 0:
      equation_land = equation_land_0;
      if (fscanf (file, "%*s%*s%Lf", &eq->tf) != 1)
        return 0;
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
      return 0;
    }
  return 1;
}

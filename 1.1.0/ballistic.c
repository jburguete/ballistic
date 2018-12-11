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
 * \file ballistic.c
 * \brief Source file with the main function.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <libxml/parser.h>
#include <glib.h>
#include "config.h"
#include "utils.h"
#include "equation.h"
#include "method.h"
#include "runge-kutta.h"
#include "multi-steps.h"

#define DEBUG_BALLISTIC 0       ///< macro to debug the ballistic functions.

long double convergence_factor;
///< convergence factor.
unsigned int ntrajectories;
///< number of projectil trajectories to calculate.
unsigned int convergence;
///< number of convergence steps.

/**
 * Function to read the basic input data.
 *
 * \return 1 on success, 0 on error.
 */
static int
convergence_read_xml (xmlNode * node)      ///< XML node.
{
	const char *message[] = {
		"Bad trajectories number",
		"Bad convergence steps",
		"Bad convergence factor"
	};
	int e, error_code;

#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_read_xml: start\n");
#endif
  ntrajectories = xml_node_get_uint (node, XML_TRAJECTORIES, &error_code);
  if (error_code || !ntrajectories)
	  {
			e = 0;
      goto fail;
		}
  convergence = xml_node_get_uint (node, XML_CONVERGENCE, &error_code);
  if (error_code || !convergence)
	  {
			e = 1;
      goto fail;
		}
  convergence_factor = xml_node_get_float (node, XML_FACTOR, &error_code);
  if (error_code || convergence_factor <= 0.)
	  {
			e = 2;
      goto fail;
		}
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_read_xml: success\n");
  fprintf (stderr, "convergence_read_xml: end\n");
#endif
  return 1;

fail:
	error_add (message[e]);
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_read_xml: error\n");
  fprintf (stderr, "convergence_read_xml: end\n");
#endif
  return 0;
}

/**
 * Function to open a numerical method on a XML node.
 *
 * \return 0 on error, 1 on Runge-Kutta method, 2 on multi-steps method.
 */
int
method_open_xml (MultiSteps * ms,
		             RungeKutta * rk,
								 xmlNode * node)
{
	char *message[] = {
    "No numerical method XML node",
    "Bad Runge-Kutta data",
    "Bad multi-steps data",
    "Unknown numerical method"
	};
	int e, m;

#if DEBUG_BALLISTIC
  fprintf (stderr, "method_open_xml: start\n");
#endif
  if (!node)
    {
      e = 0;
      goto fail;
    }
  if (!xmlStrcmp (node->name, XML_RUNGE_KUTTA))
    {
      if (!runge_kutta_read_xml (rk, node))
        {
          e = 1;
          goto fail;
        }
      runge_kutta_init_variables (rk);
			m = 1;
    }
  else if (!xmlStrcmp (node->name, XML_MULTI_STEPS))
    {
      if (!multi_steps_read_xml (ms, node))
        {
          e = 2;
          goto fail;
        }
      multi_steps_init_variables (ms);
			m = 2;
    }
  else
    {
      e = 3;
      goto fail;
    }
#if DEBUG_BALLISTIC
  fprintf (stderr, "method_open_xml: success\n");
  fprintf (stderr, "method_open_xml: end\n");
#endif
	return m;

fail:
	error_add (message[e]);
#if DEBUG_BALLISTIC
  fprintf (stderr, "method_open_xml: error\n");
  fprintf (stderr, "method_open_xml: end\n");
#endif
	return 0;
}

/**
 * Function to show the error message.
 */
void
show_error ()
{
  printf ("ERROR!\n%s", error_message);
  g_free (error_message);
	error_message = NULL;
}

/**
 * Function to print the solution values.
 */
void
print_solution (char *label,    ///< label.
                long double *r0,        ///< position vector.
                long double *r1)        ///< velocity vector.
{
  printf ("%s\n", label);
  printf ("x = %.19Le\n", r0[0]);
  printf ("y = %.19Le\n", r0[1]);
  printf ("z = %.19Le\n", r0[2]);
  printf ("vx = %.19Le\n", r1[0]);
  printf ("vy = %.19Le\n", r1[1]);
  printf ("vz = %.19Le\n", r1[2]);
}

/**
 * Function to print the numerical errors.
 */
void
print_error (char *label,       ///< label.
             long double *r1,   ///< numerical solution vector.
             long double *r2)   ///< analytical solution vector.
{
  printf ("%s = %.19Le\n", label, distance (r1, r2));
}

/**
 * Function to perform a convergence analysis of a method.
 *
 * \return 0 on success, error code on error.
 */
int
convergence_run (xmlNode * node,   ///< XML node.
                 char *output)  ///< results file name.
{
	const char *message[] = {
		NULL,
		"Bad convergence data",
		"No equation XML node",
		"Unknown numerical method",
		"Bad numerical method data"
	};
  MultiSteps ms[1];
  RungeKutta rk[1];
  Equation eq[1];
  Method *m;
  gsl_rng *rng;
  FILE *file;
  long double sr0[3], sr1[3];
  long double t, l0r0, l2r0, l0r1, l2r1, e;
	int er, me;
  unsigned int i, j;
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_run: start\n");
#endif
	er = 0;
  if (!convergence_read_xml (node))
	  {
      er = 1;
			goto fail;
		}
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_run: initing equation\n");
#endif
	node = node->children;
	if (!node)
	  {
			er = 2;
			goto fail;
		}
  if (!equation_read_xml (eq, node, 0))
	  {
			er = 3;
      goto fail;
		}
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_run: initing method\n");
#endif
	node = node->next;
	me = method_open_xml (ms, rk, node);
	switch (me)
	  {
		case 1:
      m = RUNGE_KUTTA_METHOD (rk);
			break;
		case 2:
      m = MULTI_STEPS_METHOD (ms);
			break;
		default:
			er = 4;
			goto fail;
		}
  rng = gsl_rng_alloc (gsl_rng_taus2);
  file = fopen (output, "w");
  for (j = 0; j < convergence; ++j)
    {
      gsl_rng_set (rng, 0l);
      nevaluations = 0l;
      l0r0 = l2r0 = l0r1 = l2r1 = 0.L;
      for (i = 0; i < ntrajectories; ++i)
        {
#if DEBUG_BALLISTIC
          fprintf (stderr, "convergence_run: initing equation data\n");
#endif
          equation_init (eq, rng);
#if DEBUG_BALLISTIC
          fprintf (stderr, "convergence_run: initing variables\n");
#endif
          equation_solution (eq, r0, r1, 0.);
          equation_acceleration (eq, r0, r1, r2, 0.L);
#if DEBUG_BALLISTIC
          fprintf (stderr, "convergence_run: running\n");
#endif
          if (me == 1)
            t = runge_kutta_run (rk, eq);
          else
            t = multi_steps_run (ms, eq);
#if DEBUG_BALLISTIC
          fprintf (stderr, "convergence_run: solutions\n");
          print_solution ("Numerical solution", r0, r1);
          printf ("Time = %.19Le\n", t);
#endif
          switch (eq->land_type)
            {
            case 0:
              equation_solution (eq, sr0, sr1, eq->tf);
              break;
            default:
              t = equation_solve (eq, sr0, sr1);
            }
#if DEBUG_BALLISTIC
          print_solution ("Analytical solution", sr0, sr1);
          printf ("Time = %.19Le\n", t);
          print_error ("Position error", r0, sr0);
          print_error ("Velocity error", r1, sr1);
#endif
          e = distance (r0, sr0);
          l0r0 = fmaxl (l0r0, e);
          l2r0 += e * e;
          e = distance (r1, sr1);
          l0r1 = fmaxl (l0r1, e);
          l2r1 += e * e;

        }
#if DEBUG_BALLISTIC
      fprintf (stderr, "convergence_run: saving results\n");
#endif
      l2r0 = sqrtl (l2r0 / ntrajectories);
      l2r1 = sqrtl (l2r1 / ntrajectories);
      fprintf (file, "%lu %.19Le %.19Le %.19Le %.19Le %.19Le %.19Le\n",
               nevaluations, l0r0, l2r0, l0r1, l2r1, kt, m->emt);
      switch (eq->size_type)
        {
        case 0:
          dt *= convergence_factor;
          break;
        default:
          kt *= convergence_factor;
        }
      m->emt *= convergence_factor;
      if (me == 2)
        RUNGE_KUTTA_METHOD (MULTI_STEPS_RUNGE_KUTTA (ms))->emt
          *= convergence_factor;
    }
  fclose (file);
  printf ("Time = %.19Le\n", t);
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_run: deleting method\n");
#endif
  if (me == 1)
    runge_kutta_delete (rk);
  else
    multi_steps_delete (ms);
  gsl_rng_free (rng);
fail:
	if (er)
	  {
			error_add (message[er]);
		  show_error ();
		}
#if DEBUG_BALLISTIC
  fprintf (stderr, "convergence_run: end\n");
#endif
  return 0;
}

/**
 * Function to calculate a ballistic trajectory.
 *
 * \return 0 on success, error code on error.
 */
int
ballistic_run (xmlNode * node)    ///< XML node.
{
  const char *message[] = {
    NULL,
    "Unable to open the XML root element",
    "Bad XML file",
    "No equation XML node",
    "Bad equation data",
  };
  MultiSteps ms[1];
  RungeKutta rk[1];
  Equation eq[1];
  long double sr0[3], sr1[3];
  long double t;
  int e, m;
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: start\n");
#endif
  e = 0;
  node = node->children;
  if (!node)
    {
      e = 3;
      goto end;
    }
  if (!equation_read_xml (eq, node, 1))
    {
      e = 4;
      goto end;
    }
  node = node->next;
  if (!node)
    {
      e = 5;
      goto end;
    }
  if (!xmlStrcmp (node->name, XML_RUNGE_KUTTA))
    {
      if (!runge_kutta_read_xml (rk, node))
        {
          e = 6;
          goto end;
        }
      runge_kutta_init_variables (rk);
			m = 1;
    }
  else if (!xmlStrcmp (node->name, XML_MULTI_STEPS))
    {
      if (!multi_steps_read_xml (ms, node))
        {
          e = 7;
          goto end;
        }
      multi_steps_init_variables (ms);
			m = 2;
    }
  else
    {
      e = 8;
      goto end;
    }
  nevaluations = 0l;
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: initing variables\n");
#endif
  equation_solution (eq, r0, r1, 0.);
  equation_acceleration (eq, r0, r1, r2, 0.L);
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: running\n");
#endif
  switch (m)
	  {
		case 1:
			t = runge_kutta_run (rk, eq);
			break;
		default:
			t = multi_steps_run (ms, eq);
		}
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: solutions\n");
#endif
  print_solution ("Numerical solution", r0, r1);
  printf ("Time = %.19Le\n", t);
  switch (eq->land_type)
    {
    case 0:
      equation_solution (eq, sr0, sr1, eq->tf);
      break;
    default:
      t = equation_solve (eq, sr0, sr1);
    }
  print_solution ("Analytical solution", sr0, sr1);
  printf ("Time = %.19Le\n", t);
  print_error ("Position error", r0, sr0);
  print_error ("Velocity error", r1, sr1);
  printf ("Time = %.19Le\n", t);
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: deleting method\n");
#endif
  switch (m)
	  {
		case 1:
      runge_kutta_delete (rk);
			break;
		default:
      multi_steps_delete (ms);
		}
end:
#if DEBUG_BALLISTIC
  fprintf (stderr, "ballistic_run: end\n");
#endif
	if (e)
    error_add (message[e]);
  return e;
}

/**
 * Main function
 *
 * \return 0 on success, error code otherwise.
 */
int
main (int argn,                 ///< number of arguments.
      char **argc)              ///< array of argument chars.
{
	const char *message[] = {
		NULL,
    "The syntax is:\n./ballistic input_file output_file\n",
		"Unable to open the input file",
		"Bad XML root element",
	  "Bad ballistic run",
	  "Bad convergence run",
		"Unknown model"
	};
  xmlDoc *doc;
	xmlNode *node;
  int e;
#if DEBUG_BALLISTIC
  fprintf (stderr, "main: start\n");
#endif
	e = 0;
	if (argn != 3)
	  {
			e = 1;
			goto end;
		}
  xmlKeepBlanksDefault (0);
  doc = xmlParseFile (argc[1]);
  if (!doc)
	  {
			e = 2;
			goto end;
		}
  node = xmlDocGetRootElement (doc);
  if (!node)
    {
      e = 3;
      goto end;
    }
  if (!xmlStrcmp (node->name, XML_BALLISTIC))
	  {
      e = ballistic_run (node);
			if (e)
			  {
					e = 4;
					goto end;
				}
		}
	else if (!xmlStrcmp (node->name, XML_CONVERGENCE))
	  {
      e = convergence_run (node, argc[2]);
			if (e)
			  {
					e = 5;
					goto end;
				}
		}
	else
	  e = 6;
  xmlFreeDoc (doc);
end:
  if (e)
    {
			error_add (message[e]);
      show_error ();
    }
  return e;
}

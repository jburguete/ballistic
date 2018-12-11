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
 * \file utils.h
 * \brief Header file with the useful data and functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef UTILS__H
#define UTILS__H 1

extern char *error_message;

long double distance (long double *r1, long double *r2);
long double solve_quadratic_reduced (long double a, long double b,
                                     long double x1, long double x2);
long double solve_quadratic (long double a, long double b, long double c,
                             long double x1, long double x2);
long double solve_cubic_reduced (long double a, long double b, long double c,
                                 long double x1, long double x2);
long double solve_cubic (long double a, long double b, long double c,
                         long double d, long double x1, long double x2);
int xml_node_get_int (xmlNode * node, const xmlChar * prop, int *error_code);
unsigned int xml_node_get_uint (xmlNode * node, const xmlChar * prop,
                                int *error_code);
unsigned int xml_node_get_uint_with_default (xmlNode * node,
                                             const xmlChar * prop,
                                             unsigned int default_value,
                                             int *error_code);
long double xml_node_get_float (xmlNode * node, const xmlChar * prop,
                                int *error_code);
long double xml_node_get_float_with_default (xmlNode * node,
                                             const xmlChar * prop,
                                             long double default_value,
                                             int *error_code);

#endif

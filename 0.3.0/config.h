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
 * \file config.h
 * \brief Header file with the configuration data.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2018.
 */
#ifndef CONFIG__H
#define CONFIG__H 1

#define G 9.81L                 ///< gravitational constant.

#define XML_ALPHA       (const xmlChar*)"alpha"
///< XML alpha label.
#define XML_BETA        (const xmlChar*)"beta"
///< XML beta label.
#define XML_DT          (const xmlChar*)"dt"
///< XML dt label.
#define XML_EQUATION    (const xmlChar*)"equation"
///< XML equation label.
#define XML_ERROR_TIME  (const xmlChar*)"error_time"
///< XML error-time label.
#define XML_G           (const xmlChar*)"g"
///< XML g label.
#define XML_KT          (const xmlChar*)"kt"
///< XML kt label.
#define XML_LAMBDA      (const xmlChar*)"lambda"
///< XML lambda label.
#define XML_LAND        (const xmlChar*)"land"
///< XML land label.
#define XML_MULTI_STEPS (const xmlChar*)"multi-steps"
///< XML multi-steps label.
#define XML_RUNGE_KUTTA (const xmlChar*)"runge_kutta"
///< XML runge-kutta label.
#define XML_STEPS       (const xmlChar*)"steps"
///< XML steps label.
#define XML_T           (const xmlChar*)"t"
///< XML t label.
#define XML_TIME_STEP   (const xmlChar*)"time-step"
///< XML time-step label.
#define XML_TYPE        (const xmlChar*)"type"
///< XML type label.
#define XML_VX          (const xmlChar*)"vx"
///< XML vx label.
#define XML_VY          (const xmlChar*)"vy"
///< XML vy label.
#define XML_VZ          (const xmlChar*)"vz"
///< XML vz label.
#define XML_WX          (const xmlChar*)"wx"
///< XML wx label.
#define XML_WY          (const xmlChar*)"wy"
///< XML wy label.
#define XML_X           (const xmlChar*)"x"
///< XML x label.
#define XML_Y           (const xmlChar*)"y"
///< XML y label.
#define XML_Z           (const xmlChar*)"z"
///< XML z label.

#endif

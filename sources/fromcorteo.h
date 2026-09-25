/*********************************************************************

    Copyright 2019, Christian Borschel

    This file is part of iradina.

    iradina is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iradina is distributed WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with iradina.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************/


/*********************************************************************/
/* This file contains functions and declarations, taken or adapted with
   small changes from Corteo (version Corteo20090527). Corteo was
   written by Francois Schiettekatte.
   Corteo was released under the GNU General Public License as
   published by the Free Software Foundation, version 3.
   
   You may obtain the original corteo source code from:
   http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo
*/
/*********************************************************************/

/*********************************************************************/
/* Note on the Corteo scattering database:
   If you want to compile iradina to support the scattering matrix version
   with 4 mantissa bits (Corteo before 2013), you need to include the
   following line in this file (see below):
   #include "indexvalues.h"
   If you want to compile iradina to support the version with 6 mantissa
   bits (liek corteo20130715), change the line to:
   #include "indexvalues6bit.h"
*/
/*********************************************************************/

#ifndef FROMCORTEO_H
#define FROMCORTEO_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>

#include "iradina.h"

/*************** Adapted from corteo.h ***********************/

/* prime numbers for lists length */
#define MAXLOGLIST 104729
#define MAXAZILIST 72211
#define MAXRANLIST 1000003


float       randomlist[MAXRANLIST];     /* list of evenly distributed but randomly ordered values between 0 and 1 */
float       sqrtrandomlist[MAXRANLIST]; /* sqrt of randomlist */
float       sqrtloglist[MAXLOGLIST];    /* list of evenly distributed but randomly ordered values of sqrt of -log of 1/MAXLOGLIST to 1 */
float       sqrtloglist1[MAXLOGLIST];   /* 1/sqrtloglist */
float       sinAzimAngle[MAXAZILIST];   /* list cos and sin components of angles... */
float       cosAzimAngle[MAXAZILIST];   /*   ...this angle are evenly distributed but randomly ordered between 0 and 2*PI */

void fillCosSinTable(float * cosTable, float * sinTable, float mr); /* corteo function to calc cos and sin of scattering angles */

void rotate(float *l, float *m, float *n, unsigned int * iazimAngle, float costheta, float sintheta); /* corteo function to rotate the flying direction as function of scattering angle */

void computelists();

/*************** Adapted from constants.h ***********************/

#define SCREENCONST 0.46848 /* 0.8853*0.5291772108 // screening length constant [A] */
#define E2 14.3996445f      /* e^2 / 4 pi eps0 = e^2 c^2 in [eV][A] */
#define PI 3.1415926535897932384626433832795
#define AMUbyE 1.036426883E-8 /* amu/e = 1.660538782E-27/1.602176487E-19 */
#define ELEMENTARY_CHARGE 1.602176487E-19 

#define NUMBERELEMENTS 101
float MostAbundantIsotope[NUMBERELEMENTS];
float AtomicMass[NUMBERELEMENTS];
char  AtomicNames[NUMBERELEMENTS][3];

/*************** Adapted from corteoindex.h ***********************/

/* Global variables: number of errors during index fuction calls (in work if INDEX_BOUND_CHECKING is defined) */
unsigned long EminErr, EmaxErr, SminErr, SmaxErr, DminErr, DmaxErr; 

#define INDEX_BOUND_CHECKING  /* Perform bound cheking during index evaluation (whole program 10-15% slower).
                                 No check may result out-of-bound matrix access: crash or wrong value, */

/* #include "indexvalues6bit.h" */
#include "indexvalues.h"

/* functions that compute an index from a value */
unsigned int Eindex(float Eval);
unsigned int Sindex(float Sval);
unsigned int Dindex(float Dval);

/* functions that return the value corresponding to an index */
float Eval(unsigned int index);
float Sval(unsigned int index);
float Dval(unsigned int index);

/*************** Adapted from corteomatrix.h ***********************/

int calcMatrix(int showProgress, char* DirectoryData);
int loadMatrix(char* DirectoryData);
float Matrix(unsigned long i);
void setMatrix(unsigned long i, float val);
double THETA(double epsilon, double s, unsigned int nsum);
double crossSectionScreenPot(double E, unsigned int Z1, unsigned int Z2, double massRatio, double thetaCM, unsigned int screeningType);

/*************** Adapted from corteoutil.h ***********************/

float myInvSqrtTableExp[256];
float mySqrtTableExp[256];
float myInvSqrtTable[1<<16];
float mySqrtTable[1<<16];

void  mySqrtTableFill();
float myInvSqrt(float val);
float mySqrt(float val);

void randomizelist(float *list, unsigned int maxlist);
void ignoreline(FILE *ifp);
float d2f(double val);
float sqrtdf(double val); 
float a2f(char * s);

/*************** Adapted from randomx.h ***********************/

int seed1,seed2;
double randomx();

/************** Adapted from corteo20130715 (for 64 bit compatibility) ***/
unsigned int check_type_representation();

#endif

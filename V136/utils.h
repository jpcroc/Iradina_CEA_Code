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


/*********************************************************/
/* This module contains some auxiliary functions         */
/*********************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <time.h>


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <stdarg.h>

#include "iradina.h"
#include "target.h"
#include "fileio.h"
#include "transport.h"

#include "fromcorteo.h"


#define PI 3.1415926535897932384626433832795
#define MAX_STOPPING_ENTRIES_dpass 100 /* Maximum number of values for stopping table */
/* Some predefined lists: */
/* inverf list loaded from file with a fixed length (using corteo's data): */
#define MAXERFLIST 79999

#define ZMAXdpass 92

typedef union {
    char c[800];
    double d[MAX_STOPPING_ENTRIES_dpass];
} Data;

double Enerdpass[MAX_STOPPING_ENTRIES_dpass];



float inverse_erf_list[MAXERFLIST]; /* list of evenly distributed but randomly ordered values of inverse erf of -1+1/MAXERFLIST to 1-1/MAXERFLIST */
unsigned int erflist_pointer; /* points to next erf element to use */
unsigned int iazimAngle;      /* Points to next azimutal angle to choose */
unsigned int iranlist;        /* Points to next entry in the random list */
unsigned int iranloglist;     /* Points to next entry in the random sqrt logarithmic list */

int stopping_target_index;    /* In case of -s option: target index */

int conv_create_separate_elements; /* if 1 then, separate elements are created for each material, so some elements may appear more than once */

int store_joined_output;      /* if 1, then put output data into one file */
int dpa_output;               /* outputs are in dpa (read only if ion_dose > 0 ; dpa values are meaningful only at a given dose) */
float* dscaled;               /* Introduced by J.P.Croc., required for output scaling in case of output in dpa */
char* message_buffer;         /* String buffer for message to print out to console */

FILE* elec_stopping_fp;    /* pointer to range_3d_file */
/************************** functions ****************/
void message(int level, char* msg, ...);            /* print out message to console. */
void message_error(int err_number, char* msg, ...); /* print out standardized error message to stderr */
#ifdef DEBUG_MODE
void message_debug(char* file, int line, char* msg, ...);
#endif
 
int make_double_array(char* values, int count, double* d_array);
/* Read comma-seprated values from string and put them into
   the double array, which has #count entries */

int make_float_array(char* values, int count, float* f_array);
/* Read comma-seprated values from string and put them into
   the float array, which has #count entries. The float array
   must exist already */

int make_int_array(char* values, int count, int* i_array);
/* Read comma-seprated values from string and put them into
   the int array, which has #count entries */

int store_results(char* BaseName,int ion_number);
/* Store the results of the simulation (arrays with ditribution of
   implanted ions, defects etc.. */

int handle_cmd_line_options(int argc, char* argv[]);
/* ... */

int print_help_text();

int InitConfiguration(char* ConfigFileName); /* Read configuration from file, initialize variables etc. */

int PrepareStoppingTables(); /* Read stopping data from file and fill arrays etc. */

int PrepareStoppingTablesdpass(); /* Read DPASS stopping data from file and fill arrays etc. */

int PrepareStoppingTablesfile(); /* Read DPASS stopping data from file and fill arrays etc. */

int PrepareStragglingTables(int model); /* Create straggling tables etc. */
/* Note: there is a difference between the material and the elemental version of iradina:
   Omega is stored in the tables in the material version.
   In contrast, Omega^2 is stored in the tables in the element version. */

int load_Chu_straggling_values(); /* Load tabulated Chu straggling data */

int sum_up_material_arrays(); /* Interstitials and so on are stored for each element from each material
				 separately, but may also be interesting in sum. So this function does
				 all the summing up. */

int fill_zero(int* array, int count); /* Fills an array with zeros */

void add_int_array(int* dest, int* source, int count); /* adds array source to array dest */

float fullMAGIC(float red_ip, float red_energy, float screening_length, float(*Potential)(float r,int Z1, int Z2));

float UnivPotential(float r,int Z1, int Z2);    /* return universal potential */

int count_existing_elements(int* elementarray); /* returns the number of ones in the provided array */

int calculate_normalization_factor(int num_of_ions); /* for converting units to 1/cm^2 per 1/cm^2 */

int write_status_file(char* status_text, int ion_number); /* creates a file, that describes iradina's running state */

double MAGIC(double B, double epsilon);
/* Returns cos(theta/2) for the given scattering parameters*/

double ZBL_and_deri(double R, double* Vprime);
/* returns the ZBL potential, and via the pointer Vprime its derivative */

void CalculateRelativeTargetAtomPosition(float vx,float vy, float vz,float *px, float *py, float *pz, unsigned int iazimAngle);
/* This calculates the direction in which the target nucleus is found.
   v is the projectile velocity vector, the IP vector is returned in p components */

float signf(float f);   /* return signum(f) */
double signd(double d); /* return signum(d) */

int MaterialToElementConverter(char* OutputFile); /* Convert normal material-based input files to element-based input files */

void get_float_one_bit_smaller(float* fltInput,float* fltOutput); /* returns the largest float that is smaller than the fltInput */

int print_version_info(FILE* fp); /* print some machine-readable info on this version of iradina. */

int print_some_simulation_parameters(FILE* fp,int ion_number); /* print some information on the current simulation to the stream pointed to by fp */

int print_stopping_table(int ionZ, double ionM, int target, double e_min, double e_max, double e_step); /* print stopping table for testing */

int prepare_KP_tables2 ();
/*CROC : some initialization for modified Kinchin-Pease quick calculation of damage*/

int DensityScaleArray( int* unscl, float* scl, float conc);


#endif


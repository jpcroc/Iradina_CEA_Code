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

/**************************************************************************************************/
/* Kinchin Pease calculation of damage was added by Jean-Paul Crocombette (CROC) from CEA Saclay. */
/* Some modifications were introduced by Christian Van Wambeke (CVW) from CEA Saclay in order to  */
/* support the CEA GUI.                                                                           */
/**************************************************************************************************/

/* ***********************************************************************
   This module declares some general function of iradina and some global
   variables.
   Some switsches are defined which make iradina safer but slower. */
/************************************************************************/


#ifndef IRADINA_H
#define IRADINA_H

/************************************************************************/
/* if the following parameter is defined, then iradina may account for
   special geometries deviating from the rectangular grid. 
   The special geometry is accounted for by the functions in geometry.c */

//#define INCLUDE_SPECIAL_GEOMETRY

/* if this parameter is not defined, then iradina is a little faster,
   because some queries will not be compiled into the code */
/************************************************************************/


/* Define some more options that need to be known during compilation.
   Switching these on makes iradina slower but safer */
/* #define SAFE_ROTATION */          /* makes iradina slower but SAFER */
#define SAFE_SQR_ROTATION            /* makes iradina slower but SAFER */
//#define CHECK_NAN_VECTORS            /* makes iradina not so slow as the par above but almost as safe against nan velocities */
#define RENORM_VELOCITIES            /* reduce vel unit vec to 1 */
#define PROJ_HANGUP_SAFETY 100000000 /* if a recoil or ion hangs up, this ensure that it will be stopped after so many steps.
					This should be a large number if you simulate large structures. */

//#define DEBUG_MODE                 /* If defined, then some DEBUG messages are printed. */
//#define DEBUG_MODE2                /* If defined, then some other DEBUG messages are printed. */
//#define DEBUG_MODE3                /* ... */
//#define DEBUG_MODE4
#define MONITOR_ION 1000000    /* Lowest ion number for which debug messages are printed */

/* include some standard c libraries: */
#include <stdio.h>
#include <time.h>

/* include other modules: */
#include "utils.h"
#include "target.h"
#include "fileio.h"
#include "transport.h"

#ifdef INCLUDE_SPECIAL_GEOMETRY
#include "geometry.h"
#endif

#define COMPILEDATE __DATE__
#define COMPILETIME __TIME__
#define VERSIONDATE "2019-AUG-23"
#define VERSION 1
#define SUBVERSION 2
#define SUBSUBVERSION 3
#define RELEASESTRING " "
#define VERSIONCOMMENT "  "
//#define RELEASESTRING ""
//#define VERSIONCOMMENT ""

#define MAX_FILENAME_LENGTH 1024

/* -- Some global variables -- */

int print_level;               /* Determines how much stuff is printed to the console.
				  0 means normal, positive values mean more, negative less */
int mem_usage_only;            /* if 1, then do not simulate anything, just estimate memory usage */
int mem_usage_details;         /* if 1, print more details for memory usage */
unsigned long int mem_usage;   /* if mem_usage active, then sum it up here (in bytes) */

char* ConfigFileName;          /* Name of the general input config file */
char* TargetStructureFileName; /* Filename of the file that define the structure of the target */
char* OutputFileBaseName;      /* All outputfiles begin with this name, so that's put them in a distinct directory if included in this name */
char* MaterialsFileName;       /* Name of the file that defines the materials in the target */
char* ConversionFileName;      /* Name of the converted input file (when converting from materials to elements) */
char* ElementsFileName;        /* Name of the file that defines the elements in the target.
				  This is not needed for standard material based operation,
				  but for conversion of one to another... so it is defined in both cases. */
char* DirectoryData;           /* Name of the data corteo directory */
char* TargetDensityMultFileName; /* */

int wait_before_end;           /* if 1, the program only quits after pressing return key */
int normalize_output;          /* if 1, the output of the program will be normalized to (1/cm^3)/(ions/cm^2) */
double unit_conversion_factor; /* This converts output to unit of (1/cm^3)per(1/cm^2) */
double minE; 
int create_status_file;        /* if this is 1, then iradina output its status regularly to a file that can be read by another program */
char* start_id_string;         /* string that helps other programs to identify iradina output */
int do_not_store_damage;       /* if this is 1, then far fewer data are stored to disc... saves memory */
int no_headers_in_files;       /* if this is 1, then files have no header. Default=0 */
int store_info_file;           /* if 1, then iradina creates an output file with some information on the simulation and the other outputfiles */
FILE* ball_stopping_fp;    /* pointer to file */
/******************************************************************************/
/* Function declarations:                                                     */
/******************************************************************************/

int main();                    /* The program */
int display_startup_message(); /* does what you think it does */



struct tm * timeinfo;

#endif

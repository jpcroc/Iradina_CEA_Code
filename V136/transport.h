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


/*****************************************************************************/
/* This module contains the functions to simulate the ion transport, etc.    */
/*****************************************************************************/


#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <stdio.h>
#include <math.h>
#include <time.h>
/* #include <unistd.h> */

#include "fileio.h"
#include "utils.h"
#include "iradina.h"
#include "target.h"
#include "fromcorteo.h"


#ifdef INCLUDE_SPECIAL_GEOMETRY
#include "geometry.h"
#endif

/* -- Declare some global variables (not nice, but fast)-- */


/* Parameters of the ion beam */

int simulation_type;     /* How to do the simulation.
			    0 = Full Damage Cascade, follow recoils, free flight path statistically distr.
			    3 = Ions only. Recoils are not followed, No damage profiles stored (--- not deeply TESTED!)
			    4 = KP quick calculation of damage, mono-elemental formula similar to SRIM (added by J.-P. Crocombette)
			    5 = KP quick calculation of damage, material averaging, more physical as 4 but different compared to SRIM (added by J.-P. Crocombette) */

int detailed_sputtering; /* Needs to be set to 1, if you want to get better results regarding sputtering.
			    If set to 0, the program doesn't care about surface binding energy at all
			    However, this takes some calculation time, because in the 3d geometry, it
			    it takes some time to find out whether a projectile would be ejected to
			    vacuum or not. */
int store_energy_deposit; /* if 1, array with deposited energy are created and stored.
			     if 0, not */

int flight_length_type;  /* Flight lengths between collisions can be selected
			    by three different approaches:
			    0: poisson distributed flight lengths mit average of interatomic spacing
			    1: always mean interatomic spacing of current material
			    2: constant flight length. Is has be be specified in units of nm.
			    Options 0 and 1 ignore the flight_length_constant parameter.
			    3: SRIM-like flight path */
float flight_length_constant; /* If constant flight length is selected, then this is it. (in nm) */

int ionZ;                  /* Proton number */
float ionM;                /* Mass of the ion */
float ionInitialEnergy;    /* Impinging energy */

float ion_vx;              /* vector of ion velocity, normalized to 1 */
float ion_vy;              /* Note, that the ion is NOT described by its actual velocity vector, */
float ion_vz;              /* but by the flying direction vector of length 1 and its energy */

float min_energy;          /* Minimum energy below which all projectiles are stopped */

int ion_distribution;      /* 0 for random ion entry positions, 1 for centered, 2 for specified position,
			      3 for random square around position, 4 center of simulation volume */
float enter_x,enter_y,enter_z;     /* entry point in nm */
float beam_spread;         /* in nm, only relevant for option 3 */

int max_no_ions;           /* Maximum number of ions */
int display_interval;      /* Display status every so many ions */
int override_max_ions;     /* For command-line argument override of maximum number of ions */
double override_energy;    /* For command-line argument override of ion energy */
int storage_interval;      /* Dump the target arrays into the files every such
			      number of ions */
int status_update_interval; /* After so many ions, the status is written to the status file
			       (if activated by command line arg) */
int store_transmitted_ions;/* 1 if true */
int store_path_limit;      /* Only for so many ions, the exact paths and recoil cascades are stored */
int store_path_limit_recoils;  /* Separate limit for recoils. Default: same as limit of ions, if not set */
int transmission_pointer;  /* point to next free index in array */
int store_exiting_recoils; /* 1 if all exiting recoils should be stored */
int store_exiting_limit;   /* Maximum number of exiting recoils to be stored */
int store_PKA;             /* if 1, store each primary knockon atom to an output file (can be used for later extracting a PKA spectrum externally)  */ 
double ion_dose;           /* Total ion dose. If given: results will be normalized to this dose. Unit: ions/cm^2, measured perpendicular to the ion beam. */ 
double minE,totlE,totlP;
struct transmitted_ion     /* describes an ion that has been transmitted (left the target) */
{
  float x;                /* exit position x */
  float y;                /* y */
  float z;                /* z */
  float vx;               /* exit velocity unit vector */
  float vy;
  float vz;
  float energy;           /* exit energy */
};
struct transmitted_ion* transmit_list; /* List of transmitted ions */

int leaving_ions[6];       /* number of ions, leaving target in each direction */
int store_ion_paths;       /* 1 if the exact paths should be stored (interesting for debugging stuff */
FILE* ion_paths_fp;        /* the file with the ion paths will be open all the time if the paths
			      should be stored. This points to that file */
int store_recoil_cascades; /* 1 if the exact recoils cascades should be stored (interesting for debugging stuff...) */
FILE* recoil_cascades_fp;  /* the file with the recoil cascades will be open all the time if the paths
			      should be stored. This points to that file */
int store_range3d;         /* stores the final positions of implanted ions to range3d.ions (similar to the TRIM file) */
FILE* store_range3d_fp;    /* pointer to range_3d_file */
FILE* store_range3dV_fp;   /* pointer to range_3dV_file */
FILE* store_range3dI_fp;   /* pointer to range_3dI_file */
FILE* store_PKA_fp;        /* pointer to file for storing PKAs */

float chu_values[98][4];   /* Values to calculate stragling according to Chu's model; fit data from Yang et al. NIMB61(1991)149. */
int straggling_model;      /* how to calc straggling */
int stopping_model;      /* how to calc stopping (1=SRIM; 2=DPASS ; 3=READ table) */

int max_annular_coll_volumes;  /* According to W.Eckstein "Computer Simulation if Ion-Solid Interactions",
				  Springer 1991, p.93, multiple collisions in annular volumes should be allowed
				  to occur. This number +1 determines the maximum number of collisions. 0 means
				  just 1. Recommended for sputtering is 2. */

int scattering_calculation;    /* 0: corteo database, 1: MAGIC  */
int transport_type;            /* 0: accurate, 1: Fast (like corteo). For KP: use 1 */
int single_ion_sputter_yields; /* if 1, iradina will store sputter yields for single ions (at the moment those are not
				  seperated by type of sputtered particles */
int* sputter_yield_histogram;  /* array that stores single ion sputter yield histogram */
#define MAX_SPUTTERED 5000     /* maximum number of sputtered particles per single ion possible to store */

time_t sim_start_time;         /* point in time when then transport simulation was started */

/* Some counters for the transport function. These need to be in the header file, because utils.c also has references to ion_c. */
int ion_c;
int disp_c;
int miss_c;
int coll_c;
int repl_c;
int vac_c;
int int_c;
int sputter_c;
int leaving_ions_c;
int leaving_recoils_c;
int escape_solid_c;
long int recoil_counter;   /* consecutive number for identifying each recoil in the cascade-file */

int single_ion_sputter_counter;

int replacer_escaped;

/******************************************************************************/
/* Now come the function declarations: */
/******************************************************************************/

int IrradiateTarget();   /* Let ions impinge on the target */

int FastProjectileTransport(int ProjZ, float ProjM, double ProjE, float Proj_x, float Proj_y, float Proj_z, float Proj_vx, float Proj_vy, float Proj_vz, int is_ion, int OrgMaterial, int OrgElement, int OrgCell, int RD);
/* Calculates the path of a projectile with given porperties (proton number,
   mass, energy, coordiantes, direction) through the target material.
   Can be called recursively to follow recoils */

int FullProjectileTransport(int ProjZ, float ProjM, double ProjE, float Proj_x, float Proj_y, float Proj_z, float Proj_vx, float Proj_vy, float Proj_vz, int is_ion, int OrgMaterial, int OrgElement, int OrgCell, char ProjState, int RD);
/* Like FastProjectileTransport but more accurate. Use this for sputtering
   The ProjState char contains status bits:
   bit 0: if 1, then the projectile is a recoile which has gained less then its displacement energy.
   That means it may fly around a little and kick other atoms, but it probably jumps back to its original locations. 
   bit 1: if 1, the projectile replaced another atom, but is still flying around a bit and creating damage or possibly
   might try to espace from the solid */

float ElectronicStopping(int ionZ, float ionM, float ionE, int material);
/* Calculates the electronic dedx for ion with Z, M, energy Z in given material [eV/A] */
float ElectronicStoppingdpass(int ionZ, float ionM, float ionE, int material);
/* Calculates the electronic dedx for ion with Z, M, energy Z in given material [eV/nm] with DPASS */
float ElectronicStoppingfile(int ionZ, float ionE, int material);
/* Calculates the electronic dedx for ion with Z,  energy Z in given material [eV/nm] from file table */
float ElectronicStraggling(int ionZ, float ionM, float ionE, int material);
  /* Calculates the electronic energy loss straggling for ion with Z, M, energy Z in given material  [eV/A] */

int CheckAndCorrectBoundary(float* dim, float target_size, float* target_size_max, int boundary);
/* Check wheter a position is inside the target for dimension x,y or z.
   If not, the position is corrected according to the boundary conditions.
   0 is returned if new position is inside target, 1 is returned if
   position is still outside the target (ion has left the target) */

int RefractProjectile2(float *vx, float* vy, float* vz, float nx, float ny, float nz, double* energy, float E_surf);
/* performs surface refraction of a projectile.
   The direction and energy (passed as refs) are adjusted.
   n is the surface normal vector and E_surf the surface binding energy.
   The surface vector must point to vacuum!
   returns 0 if projectile went into solid.
   returns 1 if projectile tried to leave solid and succeded.
   returns 2 if projectile tried to leave solid but did not succeed. */

int CalcSurfaceNormal(int old_cell, int new_cell, float* nx, float* ny, float *nz);
/* if the ion moves from one cell to another, it is necessary to know the surface
   normal between the to cells. This function calculates the nomalized surface vector
   for two cell indices.
   For rectangular cells the surface normal should consist of integers (1s or 0s),
   but for general geometries it might differ, so we will allow float values */

#endif

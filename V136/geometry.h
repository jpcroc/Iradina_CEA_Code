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


/*****************************************************************************/
/* This module can be specialized to incorporate special geometries          */
/*                                                                           */
/* This is the nanoparticle version of the geometry module                   */
/*                                                                           */
/*****************************************************************************/


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "target.h"
#include "iradina.h"

#define SPECIAL_GEOMETRY_NAME "nanoparticle"

/* The c-file must implement all the functions defined below */

int InitSpecialGeometry();
/* iradina calls this function if the special-geometry-parameter is activated.
   It can be used to load additional parameters from other files, if desired */

int special_CheckSurfaceAtom(int cellindex, float x, float y, float z, float spacing);
/* will usually call the standard function from target.h*/

int GetMaterialFromPosition(int cell, float x, float y, float z);
/* returns material at the given position. Can point to the usual procedure or might implement something else */

float special_CalcDirectionalFractionSqr(float x, float y, float z, float vx, float vy, float vz);
/* This function must return the square of the directional fraction of the energy which is perpendicular to the surface.
   This function should return 0, if the velocity does not point into vacuum but into the solid. This way we
   know if we must apply the E_surf or E_latt */

float special_CalcSurfaceNormal(int new_cell, int cell,float x, float y, float z, float* nx, float* ny, float* nz, float factor);
  /* if the ion moves from material or outside, we need surface normal. Factor can be used to change its direction */

/*****************************************************************/
/* The functions below are special for the nanoparticle module   */
/*****************************************************************/

float NP_diameter;
float NP_radius;
float NP_radius_sqr;

float target_half_x;
float target_half_y;
float target_half_z;

int NP_DataBlockReader(char* BlockName); /* Needs to be called from the ini file reader
					    while the nanoparticle structure file is read */
int NP_DataReader(char* ParName, char* ParValue);/* Needs to be called from the ini file reader
						    while the NP structure file is read */
 
#endif

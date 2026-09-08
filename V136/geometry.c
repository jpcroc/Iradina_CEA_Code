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
/* Note: since iradina version 1.0.6, the nanoparticle is centered within    */
/* the simulation volume also if the volume is much larger than the particle.*/
/*                                                                           */
/*****************************************************************************/


#include "geometry.h"

int InitSpecialGeometry(){
	/* iradina calls this function if the special-geometry-parameter is activated.
	It can be used to load additional parameters from other files, if desired */

	int result;
	char filename[50];

	strcpy(filename,"np_structure.in");

	/* Load additional nanoparticle structure file */
	result=IniFileReader(NP_DataBlockReader, NP_DataReader, filename);
	if(result!=0){
		message_error(result,"reading NP structure file failed %s.\n",filename);
		return result;
	}

	/* Calculate some more parameters */
	NP_radius     = NP_diameter*0.5;
	NP_radius_sqr = NP_radius * NP_radius;

	target_half_x = target_size_x/2.0;
	target_half_y = target_size_y/2.0;
	target_half_z = target_size_z/2.0;
	message(2,"Simulation volume half sizes are:\t%g\t%g\t%g nm\n",target_half_x,target_half_y,target_half_z);

	return 0;
}

int special_CheckSurfaceAtom(int cellindex, float x, float y, float z, float spacing){
	/* like GetMaterial, but compare to a reduced radius */
	float temp;
	temp = (NP_radius-spacing);
	temp = temp * temp;
	if( ((x-target_half_x)*(x-target_half_x)+(y-target_half_y)*(y-target_half_y)+(z-target_half_z)*(z-target_half_z)) < temp ) { 
		/* if( ((x-NP_radius)*(x-NP_radius)+(y-NP_radius)*(y-NP_radius)+(z-NP_radius)*(z-NP_radius)) < temp ) { */
		return 0;
	} else {
		return 1;
	}
}

int GetMaterialFromPosition(int cell, float x, float y, float z){
	if( ((x-target_half_x)*(x-target_half_x)+(y-target_half_y)*(y-target_half_y)+(z-target_half_z)*(z-target_half_z)) < NP_radius_sqr ) { /* we're inside the NP */
		/* if( ((x-NP_radius)*(x-NP_radius)+(y-NP_radius)*(y-NP_radius)+(z-NP_radius)*(z-NP_radius)) < NP_radius_sqr ) {*/ 
		return 0;
	} else {
		return 1;
	}
}

float special_CalcDirectionalFractionSqr(float x, float y, float z, float vx, float vy, float vz){
	/* This function returns the square directional fraction of the energy which is perpendicular to the surface.
	This function should return 0, if the velocity does not point into vacuum but into the solid. This way we
	know if we must apply the E_surf or E_latt */
	float temp;

	/*  x -= NP_radius;
	y -= NP_radius;
	z -= NP_radius;*/
	x -= target_half_x;
	y -= target_half_y;
	z -= target_half_z;

	temp = x*vx+y*vy+z*vz;
	if(temp<0){ /* v-vector points into the NP */
		return 0;
	} else {
		return temp*temp/(x*x+y*y+z*z);
	}
}

float special_CalcSurfaceNormal(int new_cell, int cell,float x, float y, float z, float* nx, float* ny, float* nz, float factor){
	/* if the ion moves from material or outside, we need surface normal */
	/* Factor can be used to change its direction */
	float temp;
	/*x -= NP_radius;
	y -= NP_radius;
	z -= NP_radius;*/
	x -= target_half_x;
	y -= target_half_y;
	z -= target_half_z;

	temp = factor/sqrtf(x*x + y*y + z*z); /* inverse factorized length */
	*nx=x*temp;
	*ny=y*temp;
	*nz=z*temp;
	return 0;
}

int NP_DataBlockReader(char* BlockName){
	return 0;
}

int NP_DataReader(char* ParName, char* ParValue){
	/* reads the data from the NP structure input file */
	if(strcmp(ParName,"NP_diameter")==0){ /* The NP diameter */
		sscanf(ParValue,"%f",&NP_diameter);
		message(1,"NP diameter:\t %g nm\n",NP_diameter);
		message(1,"NP material:\t %s\n",ListOfMaterials[0].Name);
	}
	return 0;
}

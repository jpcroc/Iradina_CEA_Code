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
/* This module contains the target-related functions, etc.                   */
/*****************************************************************************/


#ifndef TARGET_H
#define TARGET_H

#include <stdio.h>
#include <math.h>

#include "fileio.h"
#include "utils.h"
#include "iradina.h"
#include "transport.h"

/* -- Declare some global variables (not nice, but easier)-- */

/*******************************************************************************/
/*	Different "materials" can be defined. A material can be single-elemental
	as well as a compound.
	The target consists of (possibly a large number of) little cells. Each of
	these is filled with one of the defined materials */
/*******************************************************************************/


#define MAX_EL_PER_MAT 10         /* Maximum number of elements per material */
#define MAX_NO_MATERIALS 20       /* Maximum number of different materials */

#define MAX_STOPPING_ENTRIES 1000 /* Maximum number of values for stopping table */
#define MAX_STOPPING_ENTRIES_dpass 100 /* Maximum number of values for stopping table */
#define MAX_ELEMENT_NO 93         /* Maximum number of elements in simulation +1 */


struct material                             /* all properties of a material */
{
	char  Name[50];                           /* A name can be defined */
	int   Is_Vacuum;                          /* 1 if the material is vacuum, 0 else */
	float Density;                            /* Total atomic density in at/cm^3 */
	float DensityNM;                          /* Total atomic density in at/nm^3 */
	float AtomicDistance;                     /* Average inter-atomic distance [nm] */
	float LayerDistance;                      /* Layer distcance assuming simple cubic structure[nm] */
	float SqrtAtomicDistance;                 /* Square root of average inter-atomic distance [nm] */
	float MeanImpactPar;                      /* Not actually the mean impact paramater but rather 1/sqrt(PI*density*MeanFreePath) */
	float SqrtRecFlDensity;                   /* 1/sqrt(pi*flight_length_constant*density) needed for constant flight lengths */
	int   ElementsZ[MAX_EL_PER_MAT];          /* List of the elements contained in this materials */
	float ElementsM[MAX_EL_PER_MAT];          /* List of the element masses in this materials */
	float ElementsConc[MAX_EL_PER_MAT];       /* Relative concentration of the elements */
	float ElementsDispEnergy[MAX_EL_PER_MAT]; /* Displacement energy for each element in eV */
	float ElementsLattEnergy[MAX_EL_PER_MAT]; /* Lattice energy for each element in eV */
	float ElementsSurfEnergy[MAX_EL_PER_MAT]; /* Surface binding energy for each element in eV */
	float ElementsReplEnergy[MAX_EL_PER_MAT]; /* Replacement threshold energy for each element in eV (default=e_latt) */
	float IonSurfEnergy;                      /* Surface binding energy of the ion in this material. If not provided: use mean */
	int   SputterCounter[MAX_EL_PER_MAT*6];   /* Count sputtered atoms leaving the sample in each of the possible 6 directions
											Note: for all elements all 6 directions in one array, one after another: faster. */
	int   ElementCount;                       /* Number of different elements contained in this material */
	int   CellCount;                          /* Number of cells that consist of this material */
	float MeanZ;                              /* Average Z, weighted by fraction */ 
	float MeanM;                              /* Average M, weighted by fraction */ 
	float MeanF;                              /* Average F (reduced energy conversion factor), weighted by fraction */ 
	float MeanA;                              /* Average A (screening length), weighted by fraction, */ 
	float MeanEd;                             /* Average displacement energy, weighted by fraction */ 
	float MeanMinRedTransfer;                 /* Average minimun energy transfer in reduced units */

	/*CROC  for KP quick calculation*/
	float ed_oE ;              /* ToDo: add comment */
	float k_d ;                /* ToDo: add comment */

	/*	For each element of each material we want to store the recoils, interstitials and so on. The equations are:
		Displacements = Vacancies + Replacement Collisions
		Vacancies     = Interstitials + Atoms that left the target */
	/*	So will need to make an array of pointers (one pointer for each element), pointing to
		the various target arrays. Further down, these four arrays occur once for the target but 
		independently of the elements, and those will hold the sum for all elements */
	int **TargetImplantedRecoilsInt;   /* The implanted recoils stopped as interstitials */
	int **TargetImplantedRecoilsRepl;  /* The implanted recoils stopped as replacement atoms */
	int **TargetElementalVacancies;    /* Vacancies of this element left behind */
	int **TargetElementalDisp;         /* Displacements of this element that took place */
	int **TargetSputteredAtoms;        /* Number of atoms sputtered from this cell */
	float** StoppingZEdpass;                /* points to an array of 92 elements (the Zs) which contain pointers to logarithmically scaled stopping tables for Z in this material.*/
        float** StoppingZEfile;                /* points to an array of 92 elements (the Zs) which contain pointers to logarithmically scaled stopping tables for Z in this material.*/
        float** Efile;                /* points to an array of 92 elements (the Zs) which contain pointers to logarithmically scaled stopping tables for Z in this material.*/
	float** StoppingZE;                /* points to an array of 92 elements (the Zs) which contain pointers to logarithmically scaled stopping tables for Z in this material.*/  
	float** StragglingZE;              /* We need the same for energy loss straggling */
	struct transmitted_ion** ElementalLeavingRecoils; /* Arrays storing recoils that are leaving of each element */
	int *leaving_recoils_pointer;      /* Array of pointers which for each element of this material point to next free leaving position */
        int*  StoppingNLfile
};
struct material ListOfMaterials[MAX_NO_MATERIALS];
int NumberOfMaterials;               /* Points to first free index in the ListOfMaterials */
int number_of_all_target_elements;   /* without ion. */

int existing_elements[MAX_ELEMENT_NO];  /* This array holds a 1 for any element that might exist in the target,
					so all these elements might occur as projectiles. */
int hydrogen_in_target;                 /* Hydrogen is always included in the existing_elements_list, because
					we need its stopping powers. However, we do only need it in the element
					file, if it is really in the target.*/
int ionZ_in_target;                     /* The ion is always included in the existing_elements_list.
					However, we only need it as an extra element in the element file,
					if it is really in the target.*/

struct scattering_matrix     /* This structure holds information for a scattering events of two 
				particles, for each possible energy and impact parameter */
{
	float screening_length;    /* screening length from ZBL85,p45, eq.2-60, but in [nm] */
	float inv_screening_length;/* 1/screening length in 1/nm */
	float mass_ratio;          /* ... */
	float sqrt_mass_ratio;     /* we will need this occasionally */
	float kfactor_m;           /* mass part of the kinematic factor. This is called EC in the TRIM code */
	float red_E_conv;          /* reduced energy conversion factor. This is called FI in TRIM */
	float* CosScat;            /* Cosines of scattering angles, for each energy and p */
	float* SinScat;            /* Sines of scattering angles, for each energy and p */
	/*float max_red_im_par;*/      /* Maximum reduced impact parameter */
	/* opposed to corteo, we need to calc this value during simulation only */
};

/*	The following structures hold scattering matrices for each possible projectile and target combination.
	The ion gets its own matrix, because the ion might be a specific isotope which is not the naturally most abundant.
	For the other elements, we will assume most abundant isotopes only, because memory is not unlimited.
	In case we need to care about specific target isotopes, the MAGIC algorithm should be used instead of the ScatBase approach. */
struct scattering_matrix ion_scattering_matrix[MAX_ELEMENT_NO];
struct scattering_matrix scattering_matrices[MAX_ELEMENT_NO][MAX_ELEMENT_NO]; /* First index: projectile, second: target */


/*******************************************************************************/
/*	The target is a cuboid-shaped box, which consists of (possibly a large number
	of) of small equal-sized cuboid cells. Several arrays describe the target
	composition, the distributions of implanted ions, recoils, vacancies and so
	on. The dimension of the target is defined by the number and size of the
	small cells in each direction.
	Furthermore, boundary conditions need to be defined in each target direction.
	x is the primary incident ion beam direction. (Though the ion beam may hit
	the target under any angle).
	Note, that integer coordinates always refer to a cell x, y and z number,
	while float coordinates always to refer to nm positions in the target.
*/
/*******************************************************************************/

int cell_count_x;     /* Number of cells in x-direction (>=1) */
int cell_count_y;     /* Number of cells in y-direction (>=1) */
int cell_count_z;     /* Number of cells in z-direction (>=1) */
int layer_count_yz;   /* cell_count_y * cell_count_z */
int cell_count;       /* Total number of cells */

float cell_size_x;      /* Size of cells in x-direction in nm */
float cell_size_y;      /* Size of cells in x-direction in nm */
float cell_size_z;      /* Size of cells in x-direction in nm */
double cell_volume;     /* product of the above three */

float target_size_x;    /* Size of target in x-direction in nm */
float target_size_y;    /* They are calculated */
float target_size_z; 
float target_max_x;     /* Maximum allowed position in x-direction in nm */
float target_max_y;     /* These positions are the largest possible values, */
float target_max_z;     /* which are smaller target_size_x etc. */

int boundary_x;         /* 0 if no specific boundary condition */
int boundary_y;         /* 1 means periodic boundary condition */
int boundary_z;

int special_geometry;   /* if 1, then the current material and surface are not checked
			using the rectangular cell grid, but other parameters instead */

/* We need to declare some target arrays (or actually pointers to them) */
int   *TargetComposition;       /* Material of each cell */
float *TargetDensityMult;       /* Multiplicators for material density in each cell */
int   UseDensityMult;           /* 0 if multiplicators not used, 1 if used */
int   *TargetImplantedIons;     /* Implanted ions per cell (interstitials+replacements) */
int   *TargetReplacingIons;     /* Implanted ions that replaced identical target atoms */
int   *TargetTotalVacancies;    /* Vacancies per cell (of all types) */
int   *TargetTotalReplacements; /* Replacements per cell (of all types) */
int   *TargetTotalDisplacements; /* Displacements per cell (of all types) */ /* Disp = Vac + Repl */
int   *TargetTotalInterstitials; /* Sum of interstitials per cell (of all recoil
					types, NOT including the implanted ions) */


int  *TargetTotalSputtered;     /* Sum of sputtered atoms from each cell */
char *TargetVacuumNeighbors;    /* Describes for each cell, which of its neighbors are
				vacuum (bitwise), this is helpful for determining sputtering */

int TotalSputterCounter[6];     /* Number of target atoms leaving the sample in each of the 6 directions */
char* TargetCompositionFileName;
/* char* TargetDensityMultFileName;  */
int TargetCompositionFileType;   /* The file which hold the info of what material is in
					which cell can have two types: either just one column of
					values, which are indexed like the TargetIndexFunction does;
					or: four columns with x,y,z values and the material index */

double *TargetEnergyPhonons;    /* All energy deposited into the phononic system */
double *TargetEnergyElectrons;  /* All energy deposited by electronic stopping.
				Note: these two arrays must be of type double, because small values
				might be added to gigantic number, and should not be lost (happens
				for large number of ions and high energys. Doubles reduce this risk
				as compared to floats) */

/******************************************************************************/
/* Function declarations:                                                     */
/******************************************************************************/

int InitializeMaterials(char* Filename); /* Read materials from input file and do preparatory calculations */
int Materials_DataBlockReader(char* BlockName); /* Needs to be called from the ini file reader
						while the materials input file is read */
int Materials_DataReader(char* ParName, char* ParValue);/* Needs to be called from the ini file reader
							while the materials input file is read */
int InitializeTargetStructure(char* Filename); /* Read target structure (size etc.) from config file
						and reads target concentration array etc.*/

int TargetIndex(int x, int y, int z); /* Calculates index for accessing one-dimensional target
					arrays knowing the integer cell coordinates */

int GetTargetXYZ(int index, int* x, int* y, int* z); /* Calculates the x,y,z coordinates (cell numbers)
							for a given index ... so the inverse function of TargetIndex() */ 

int GetCellIndex(float x, float y, float z); /* Calcualtes index for accessing one-dimensional
						target arrays from a float-value point in space */

int PositionInCell(float x, float y, float z, float* rx, float* ry, float* rz); 
/* Calcualtes the relative position within a cell from an absolute position x,y,z */

int TargetStructure_DataBlockReader(char* BlockName); /* Needs to be called from the ini file reader
						while the traget structure input file is read */
int TargetStructure_DataReader(char* ParName, char* ParValue);/* Needs to be called from the ini file reader
							while the target structure input file is read */

int PrepareScatteringMatrix(struct scattering_matrix * ScatMatrix, float ProjM, int ProjZ, float TarM, int TarZ);
/* Creates the matrix with scattering results and stores it to the structure pointed to by ScatMatrix */
/* This function needs to be here, to exclude some warning */

int DetermineVacuumNeighbors(); /* Creates an array where for each cell a char is stored that
				describes, which neighbors are vacuum */

int CheckSurfaceAtom(int cellindex, float x, float y, float z, float spacing);
/* Returns 1 if the given position corresponds to a surface atom 
	spacing is the maximum distance from atom to vacuum to be
	considered a surface atom */

#endif

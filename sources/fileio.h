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
/* This module contains the some utility functions for file access etc.      */
/*****************************************************************************/


#ifndef FILEIO_H
#define FILEIO_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "target.h"
#include "transport.h"
#include "iradina.h"
#include "utils.h"
#include "fromcorteo.h"


/* gloval variables */
int single_input_file;   /* if 0, then multiple files (normal). If 1, then a single input file is used */

/** functions **/
int IniFileReader(int(*DataBlockReader)(char* BlockName), int(*DataReader)(char* ParName, char* ParValue), char* filename);
/* A general reader for the various input config files. It parses an ini-like
   file. Whenever it finds a new DataBlock it calls the function pointed to
   by *DataBlockReader with the parameter BlockName. When it finds data it
   calls the function pointed to by DataReader with the parametername and its
   value. */ 

int ReadIntFileIntoArray(char* Filename,int* TargetArray, int Count, int FileType);     /* Reads the designated file into the designated array */
int ReadFloatFileIntoArray(char* Filename,float* TargetArray, int Count, int FileType); /* Reads the designated file into the designated array */
int WriteIntArrayToFile(char* Filename,int* SourceArray, int Count, int FileType);      /* Writes designated array into file */
int Write2ArraysToFile(char* Filename,int* SourceArray1,float maxA1,int* SourceArray2, float maxA2,int Count, int FileType);      /* Writes the two designated array into file */
int WriteFloatArrayToFile(char* Filename,float* SourceArray, int Count, int FileType);  /* Writes designated array into file  NO RESCALE BY unit_conveersion_factor*/
int WriteDoubleArrayToFile(char* Filename,double* SourceArray, int Count, int FileType);/* Writes designated array into file */

int ConfigFileDataBlockReader(char* BlockName);          /* reads a data block from the configuration file */
int ConfigFileDataReader(char* ParName, char* ParValue); /* reads a data block from the configuration file */

int StoreTransmissionArray(char* FileName, struct transmitted_ion* trans_array, int tr_pointer); /* Store array of transmitted ions */

FILE* OpenFileContinuous(char* BaseName, char* Extension); /* Opens file basename+extension, keeps it open */

int FloatBlockReader(char* FileName, int Offset, int Count, float* Array);
/* Reads a block of Count float values from file FileName starting at Offset
   and puts these in Array */

int WriteStringToFile(char* Filename, char* string); /* does what you think it does */

int LoadInverseErf();  /* Function adapted from corteo.c, loads a list of inverse error function erfinv(x) values */

int display_a_file(char* Filename); /* prints the contents of a text file to the std out */

int SplitSingleInputFile(char* filename); /* if a single input file is provided, split it */

int CheckSplitInputFile(char* filename);  /* checks if the given config file is a single combined input file (incl. structure, mat, etc...).
					     If so, the file is split into four separate temp files for further "conventional" processing */

int CombineFiles(int count, ...); /* Combine a number of files into one. */

int WriteFileHeader(char* filename, char * title, char * ct4 , char * cn4);
int WriteEnergyFileHeader(char* filename, char * title, char * ct4 , char * cn4);

int WriteResArrayToFile(char* Filename,int* SourceArray, int Count, int FileType, float conc);
/*CROC : deals with dpa output and header of files*/
#endif

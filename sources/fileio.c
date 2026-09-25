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
/* This module contains the some utility functions for file access etc.      */
/*****************************************************************************/

/* error numbers in this module 5000-5999 */

#include "fileio.h"


int IniFileReader(int(*DataBlockReader)(char* BlockName), int(*DataReader)(char* ParName, char* ParValue), char* filename) {
	/* A general reader for the various input config files. It parses
	and ini-like file. Whenever it finds a new DataBlock it calls the
	function pointed to by *DataBlockReader with the parameter BlockName.
	When it finds data it calls the function pointed to by DataReader
	with the parametername and its value. */ 
	/* returns 0 on success, -1 if file cannot be read, minus higher number
	on other errors... */

	FILE* ini_file;  /* Pointer to config file */
	int i=1;         /* To count lines */
	char* temp;      /* To read current line */
	char* temp_cut;  /* temp without its last letter */
	char* value;     /* Value of current parameter */
	int length;      /* To store string_lenth */
	int j;           /* Go through string */
	int equ_sign;    /* Position of =-sign in line */
	int result;      /* for storing returned values */

	equ_sign=0;

	ini_file=fopen(filename,"rt");
	if(ini_file==NULL){
		message_error(-5000,"File %s cannot be opened for reading.\n",filename);
		return -5000;
	}

	/* Reserve some space to read lines */
	temp=(char*)malloc(sizeof(char)*512);
	if(temp==NULL){return -5001;}
	temp_cut=(char*)malloc(sizeof(char)*512);
	if(temp_cut==NULL){return -5002;}
	value=(char*)malloc(sizeof(char)*512);
	if(value==NULL){return -5003;}

	while((i<10000)&&(!feof(ini_file))){ /* read 10000 lines max
					or until file ends */
		if(fgets(temp,511,ini_file)!=NULL){
			/* Now the line has been read */
			/* Check if first character is the comment sign # or lines is empty */
			/* Then do not process the line */
			if((temp[0]!='#')&&(temp[0]!='\n')){ /* Line should be processed */
				length=(int)strlen(temp);

				if(length>1){
					if(temp[0]=='['){/* Check if first sign is [. Then its a new block */
						/* Data Block handling */
						/* Extract the data block name */
						j=1;
						while((j<length)&&(temp[j]!=']')){
							value[j-1]=temp[j];
							j++;
						}
						value[j-1]='\0';
						result=DataBlockReader(value);
						if(result!=0){return result;}
					} else {  /* Assume its a data line */
						/* First, search for the "="-sign */
						j=0;
						while(j<length){
							if(temp[j]=='='){
								equ_sign=j;
								j=length;
							}
							j++;
						}

						/* Seperate the value from the current line assignment */
						for(j=equ_sign+1;j<length;j++){
							value[j-equ_sign-1]=temp[j];
							if(temp[j]=='\015'){value[j-equ_sign-1]='\0';} /* this should prevent some problems
								when reading windows-generated config
								files in unix-like systems */
						}
						value[j-equ_sign-2]='\0';
						/* Cut line before = sign */
						temp[equ_sign]='\0';
						temp_cut[0]='\0';
						/* Copy temp to temp_cut, but without last character */
						for(j=0;j<(equ_sign-1);j++){
							temp_cut[j]=temp[j];
						}
						if(equ_sign>0){temp_cut[equ_sign-1]='\0';}

						/* Let the designated function handle the data read */
						result=DataReader(temp,value);
						if(result!=0){return result;}
					}
				}
			}
		}
	}

	fclose(ini_file);

	free(temp);
	free(temp_cut);
	free(value);

	return 0;
}

int ReadIntFileIntoArray(char* Filename,int* TargetArray, int Count, int FileType){
	/* Reads the designated file into the designated array, which has count
	element. if the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */
	/* If the file does not match the expected format, the program might crash */

	FILE* fp;
	int   i;
	int   iline;
	int   x,y,z,value;
	int   result;

	fp=fopen(Filename,"rt");
	if(fp==NULL){return -5004;}

	i=0;

	switch(FileType){
	case 1: /* just one column */
		iline = 0;
		while((i<Count)&&(!feof(fp))){ /* read count lines max or until file ends */
			iline++;
			result=fscanf(fp,"%i",&value);
			if(result==1){
				(TargetArray[i++])=value;
			} else {
				if(result!=EOF){
					message_error(0,"bad entry in line %i in target composition file '%s'.\n", iline, Filename);
					return 1;
				}
			}
		}
		break;


	default: /* four column file, x,y,z and value. */
		iline = 0;
		while(!feof(fp)){ /* read count lines max or until file ends */
			iline++;
			result=fscanf(fp,"%i %i %i %i",&x,&y,&z,&value);
			if(result==4){
				i=TargetIndex(x,y,z);if((i>=0)&&(i<Count)){TargetArray[i]=value;}
			} else {
				if(result!=EOF){
					message_error(0,"bad entry in line %i in target composition file '%s'.\n", iline, Filename);
					return 1;
				}
			}
		}
		break;
	}
	fclose(fp);
	return 0;
}

int ReadFloatFileIntoArray(char* Filename,float* TargetArray, int Count, int FileType){
	/* Reads the designated file into the designated array, which has count
	element. if the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */
	/* If the file does not match the expected format, the program might crash */

	/* To do: insert check for correct fscanf reading like in the integer version */

	FILE* fp;
	int   i;
	int   x,y,z;
	float value;
	fp=fopen(Filename,"rt");
	if(fp==NULL){return -5005;}
	i=0;
	switch(FileType){
	case 1: /* just one column */
		while((i<Count)&&(!feof(fp))){ /* read count lines max or until file ends */
			fscanf(fp,"%g",&(TargetArray[i++]));
		}
		break;
	default: /* four column file, x,y,z and value. */
		while(!feof(fp)){ /* read count lines max or until file ends */
			fscanf(fp,"%i %i %i %g",&x,&y,&z,&value);
			i=TargetIndex(x,y,z);
			if((i>=0)&&(i<Count)){TargetArray[i]=value;}
		}
		break;
	}
	fclose(fp);
	return 0;
}

int WriteIntArrayToFile(char* Filename,int* SourceArray, int Count, int FileType){
	/* Writes the designated array of Count elements into a file. 
	If the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */

	FILE* fp;
	int   i;
	int   x,y,z;

//	fp=fopen(Filename,"wt");
	if (no_headers_in_files==0) {
		fp=fopen(Filename,"a");
	}else{
		fp=fopen(Filename,"wt");
	}
	if(fp==NULL){return -5006;}
	i=0;

	switch(FileType){
	case 1: /* just one column */
		if(normalize_output>=1){
			for(i=0;i<Count;i++){
				fprintf(fp,"%g\n",(double)SourceArray[i] * unit_conversion_factor);
			}
		} else {
			for(i=0;i<Count;i++){
				fprintf(fp,"%i\n",SourceArray[i]);
			}
		}
		break;
	default: /* four column file, x,y,z and value. */
		if(normalize_output>=1){
			for(i=0;i<Count;i++){
				GetTargetXYZ(i,&x,&y,&z); /* convert linear index to 3 coords */
				fprintf(fp,"%i\t%i\t%i\t%g\n",x,y,z,SourceArray[i] * unit_conversion_factor);
			}
		} else {
			for(i=0;i<Count;i++){
				GetTargetXYZ(i,&x,&y,&z); /* convert linear index to 3 coords */
				fprintf(fp,"%i\t%i\t%i\t%i\n",x,y,z,SourceArray[i]);
			}
		}
		break;
	}
	fclose(fp);
	return 0;
}
int Write2ArraysToFile(char* Filename,int* SourceArray1,float maxA1,int* SourceArray2, float maxA2,int Count, int FileType){      /* Writes designated array into file */
	/* CRC*/
	/* Writes the designated array of Count elements into a file. 
	If the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */


	FILE* fp;
	int   i;
	int   x,y,z;

	if (no_headers_in_files==0) {
		fp=fopen(Filename,"a");
	}else{
		fp=fopen(Filename,"wt");
	}
	if(fp==NULL){return -5006;}
	i=0;

	switch(FileType){
	case 1: 
		for(i=0;i<Count;i++){
			fprintf(fp,"%g\t%g\n",(double)SourceArray1[i]/(double)(maxA1),(double)SourceArray2[i]/(double)(maxA2));
		}
		break;
	default: 
		for(i=0;i<Count;i++){
			GetTargetXYZ(i,&x,&y,&z); 
			fprintf(fp,"%i\t%i\t%i\t%g\t%g\n",x,y,z,(double)SourceArray1[i]/(double)(maxA1),(double)SourceArray2[i]/(double)(maxA2));
		}
		break;
	}
	fclose(fp);
	return 0;
}



int WriteFloatArrayToFile(char* Filename, float* SourceArray, int Count, int FileType){
	/* Writes the designated array of Count elements into a file. 
	If the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */

	FILE* fp;
	int   i;
	int   x,y,z;

	if (no_headers_in_files==0) {
		fp=fopen(Filename,"a");
	}else{
		fp=fopen(Filename,"wt");
	}
	if(fp==NULL){return -5007;}
	i=0;
	switch(FileType){
	case 1: /* just one column */
		for(i=0;i<Count;i++){
			fprintf(fp,"%g\n",SourceArray[i]);
		}
		break;
	default: /* four column file, x,y,z and value. */
		for(i=0;i<Count;i++){
			GetTargetXYZ(i,&x,&y,&z); /* convert linear index to 3 coords */
			fprintf(fp,"%i\t%i\t%i\t%g\n",x,y,z,SourceArray[i]);
		}
		break;
	}
	fclose(fp);
	return 0;
}


int WriteDoubleArrayToFile(char* Filename, double* SourceArray, int Count, int FileType){
	/* Writes the designated array of Count elements into a file. 
	If the file is just one column of values then FileType should be
	set to 1. If the file contains 4 columns (x, y, z, value)
	then set it to 0.
	The caller needs to make sure, that the array is large enough. */

	FILE* fp;
	int   i;
	int   x,y,z;
	if (no_headers_in_files==0) {
		fp=fopen(Filename,"a");
	}else{
		fp=fopen(Filename,"wt");
	}
	if(fp==NULL){return -5008;}
	i=0;
	switch(FileType){
	case 1: /* just one column */
		for(i=0;i<Count;i++){fprintf(fp,"%lg\n",SourceArray[i]);
		}
		break;
	default: /* four column file, x,y,z and value. */
		for(i=0;i<Count;i++){
			GetTargetXYZ(i,&x,&y,&z); /* convert linear index to 3 coords */
			fprintf(fp,"%i\t%i\t%i\t%lg\n",x,y,z,SourceArray[i]*unit_conversion_factor);
		}
		break;
	}
	fclose(fp);
	return 0;
}

int StoreTransmissionArray(char* FileName, struct transmitted_ion* trans_array, int tr_pointer){
	/* Store array of transmitted ions */
	FILE* fp;
	int i;
	fp=fopen(FileName,"wt");
	if(fp==NULL){return -5009;}
	i=0;

	for(i=0;i<tr_pointer;i++){
		fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%i\n",
		trans_array[i].x,
		trans_array[i].y,
		trans_array[i].z,
		trans_array[i].vx,
		trans_array[i].vy,
		trans_array[i].vz,
		trans_array[i].energy,
		i);
	}
	fclose(fp);
	return 0;
}

int ConfigFileDataBlockReader(char* BlockName){
	/* reads a data block from the configuration file */
	/* ... we will ignore data blocks of the configuration file for now */
	return 0;
}

int ConfigFileDataReader(char* ParName, char* ParValue){
	/* reads data from the configuration file */
	/* Needs to be called from the ini file reader
	while the configuration input file is read */
	/* Compare the parameter name to known parameters and then
	read in corresponding value */

	int itemp=0;

	if(strcmp(ParName,"TargetstructureFileName")==0){ /* ... */
		/* TargetStructureFileName=malloc(sizeof(char)*(strlen(ParValue)+1));
	if(TargetStructureFileName==NULL){return -1;} */
		if(single_input_file!=1){
			strcpy(TargetStructureFileName,ParValue);
		}
	}
	if(strcmp(ParName,"MaterialsFileName")==0){ /* ... */
		/*MaterialsFileName=malloc(sizeof(char)*(strlen(ParValue)+1));
	if(MaterialsFileName==NULL){return -1;} */
		if(single_input_file!=1){
			strcpy(MaterialsFileName,ParValue);
		}
	}
	if(strcmp(ParName,"simulation_type")==0){ /* How detailed the simulation should be */
		sscanf(ParValue,"%i",&simulation_type);
		message(1,"Simulation type:\t\t %i\n",simulation_type);
	}

	if(strcmp(ParName,"detailed_sputtering")==0){ /* How detailed the sputtering should be evaluated */
		sscanf(ParValue,"%i",&detailed_sputtering);
		message(1,"Detailed sputtering:\t\t %i\n",detailed_sputtering);
	}

	if(strcmp(ParName,"store_energy_deposit")==0){ /* If energy deposition should be stored, this must be 1 */
		sscanf(ParValue,"%i",&store_energy_deposit);
		message(1,"Detailed energy deposition:\t %i\n",store_energy_deposit);
	}

	if(strcmp(ParName,"min_energy")==0){ /* Minimum energy for moving projectiles */
		sscanf(ParValue,"%f",&min_energy);
		message(1,"Minimum projectile energy:\t %g eV\n",min_energy);
	}

	/* Ion beam parameters */
	if(strcmp(ParName,"ionZ")==0){ /* Read proton number of ion */
		sscanf(ParValue,"%i",&ionZ);
		message(1,"ion beam Z:\t\t\t %i\n",ionZ);
	}
	if(strcmp(ParName,"ionM")==0){ /* Read mass of ion */
		sscanf(ParValue,"%f",&ionM);
		message(1,"ion beam M:\t\t\t %g amu\n",ionM);
	}
	if(strcmp(ParName,"ionE0")==0){ /* Read mass of ion */
		sscanf(ParValue,"%f",&ionInitialEnergy);
		message(1,"ion beam E_0:\t\t\t %g eV\n",ionInitialEnergy);
	}

	if(strcmp(ParName,"ion_vx")==0){ /* velocity vector of unit length */
		sscanf(ParValue,"%f",&ion_vx);
		message(1,"ion beam vx:\t\t\t %g\n",ion_vx);
	}
	if(strcmp(ParName,"ion_vy")==0){ /* velocity vector of unit length */
		sscanf(ParValue,"%f",&ion_vy);
		message(1,"ion beam vy:\t\t\t %g\n",ion_vy);
	}
	if(strcmp(ParName,"ion_vz")==0){ /* velocity vector of unit length */
		sscanf(ParValue,"%f",&ion_vz);
		message(1,"ion beam vz:\t\t\t %g\n",ion_vz);
	}

	if(strcmp(ParName,"enter_x")==0){ /* point of entry */
		sscanf(ParValue,"%f",&enter_x);
		message(1,"Entry point x:\t\t\t %g nm\n",enter_y);
	}
	if(strcmp(ParName,"enter_y")==0){ /* point of entry */
		sscanf(ParValue,"%f",&enter_y);
		message(1,"Entry point y:\t\t\t %g nm\n",enter_y);
	}
	if(strcmp(ParName,"enter_z")==0){ /* point of entry */
		sscanf(ParValue,"%f",&enter_z);
		message(1,"Entry point z:\t\t\t %g nm\n",enter_z);
	}
	if(strcmp(ParName,"beam_spread")==0){
		sscanf(ParValue,"%f",&beam_spread);
		message(1,"Beam spread:\t\t\t %g nm\n",beam_spread);
	}

	/* For compatibility reasons: */
	if(strcmp(ParName,"ion_angle_y")==0){ /* Read orientation of ion beam */
		message_error(-5010,"ion_angle_y option has been disabled!\n");
		return -5010;
	}
	/* For compatibility reasons: */
	if(strcmp(ParName,"ion_angle_z")==0){ /* Read orientation of ion beam */
		message_error(-5011,"ion_angle_z option has been disabled!\n");
		return 5011;
	}


	if(strcmp(ParName,"storage_interval")==0){ /* ... */
		sscanf(ParValue,"%i",&storage_interval);
		message(1,"Store results every \t\t %i ions.\n",storage_interval);
	}
	if(strcmp(ParName,"display_interval")==0){ /* ... */
		sscanf(ParValue,"%i",&display_interval);
		message(1,"Show ion number every \t\t %i ions.\n",display_interval);
	}
	if(strcmp(ParName,"status_update_interval")==0){ /* ... */
		sscanf(ParValue,"%i",&itemp);
		if(itemp>0){status_update_interval=itemp;}
		message(1,"Write status file every \t\t %i ions.\n",status_update_interval);
	}
	if(strcmp(ParName,"store_path_limit")==0){ /* ... */
		sscanf(ParValue,"%i",&store_path_limit);
		message(1,"Store paths or cascades only until ion %i.\n",store_path_limit);
	}
	if(strcmp(ParName,"store_path_limit_recoils")==0){ /* ... */
		sscanf(ParValue,"%i",&store_path_limit_recoils);
		message(1,"Store cascades only until ion %i.\n",store_path_limit_recoils);
	}
	if(strcmp(ParName,"max_no_ions")==0){ /* Read maximum number of ions to simulate */
		sscanf(ParValue,"%i",&max_no_ions);
		message(1,"Max number of ions:\t\t %i\n",max_no_ions);
	}
	if(strcmp(ParName,"store_transmitted_ions")==0){ /* Store transmitted ions? */
		sscanf(ParValue,"%i",&store_transmitted_ions);
		message(1,"Storing transmitted ions?\t %i\n",store_transmitted_ions);
	}
	if(strcmp(ParName,"store_exiting_recoils")==0){ /* Store exiting recoils? */
		sscanf(ParValue,"%i",&store_exiting_recoils);
		message(1,"Storing exiting recoils?\t %i\n",store_exiting_recoils);
	}
	if(strcmp(ParName,"store_exiting_limit")==0){ /* Store transmitted recoils until this number */
		sscanf(ParValue,"%i",&store_exiting_limit);
		message(1,"Max. exiting recoils:\t\t%i\n",store_exiting_limit);
	}

	if(strcmp(ParName,"store_ion_paths")==0){ /* Store ion paths? */
		sscanf(ParValue,"%i",&store_ion_paths);
		message(1,"Storing ions paths?\t\t %i\n",store_ion_paths);
	}
	if(strcmp(ParName,"store_recoil_cascades")==0){ /* Store recoil cascades? */
		sscanf(ParValue,"%i",&store_recoil_cascades);
		message(1,"Storing exact recoils cascades?\t %i\n",store_recoil_cascades);
	}
	if(strcmp(ParName,"store_PKA")==0){ /* Store PKAs? */
		sscanf(ParValue,"%i",&store_PKA);
		message(1,"Logging primary knockons?\t %i\n",store_PKA);
	}
	if(strcmp(ParName,"store_range3d")==0){ /* Store end positions of implanted ions */
		sscanf(ParValue,"%i",&store_range3d);
		message(1,"Storing final ion positions?\t %i\n",store_range3d);
	}
	if(strcmp(ParName,"store_info_file")==0){ /* Store more information */
		sscanf(ParValue,"%i",&store_info_file);
		message(1,"Store information on simulation?\t %i\n",store_info_file);
	}
	if(strcmp(ParName,"store_joined_output")==0){ /* Store output values into one big table */
		sscanf(ParValue,"%i",&store_joined_output);
		message(1,"Store output into joined table?\t %i\n",store_joined_output);
	}
	if(strcmp(ParName,"OutputFileBaseName")==0){ /* ... */
		OutputFileBaseName=malloc(sizeof(char)*(strlen(ParValue)+1));
		if(OutputFileBaseName==NULL){return -5012;}
		strcpy(OutputFileBaseName,ParValue);
	}
	if(strcmp(ParName,"seed1")==0){ /* ... */
		sscanf(ParValue,"%i",&seed1);
		if(seed1==0){message_error(-5031,"invalid random seed1=%i\n",seed1);return -5031;}
		message(1,"Random seed 1 set to:\t\t %i.\n",seed1);
	}
	if(strcmp(ParName,"seed2")==0){ /* ... */
		sscanf(ParValue,"%i",&seed2);
		if(seed2==0){message_error(-5031,"invalid random seed2=%i\n",seed2);return -5031;}
		message(1,"Random seed 2 set to:\t\t %i.\n",seed2);
	}
	if(strcmp(ParName,"ion_distribution")==0){ /* ... */
		sscanf(ParValue,"%i",&ion_distribution);
		message(1,"Ion entry distribution model:\t %i.\n",ion_distribution);
	}

	if(strcmp(ParName,"straggling_model")==0){ /* ... */
		sscanf(ParValue,"%i",&straggling_model);
		message(1,"Straggling model is:\t\t %i\n",straggling_model);
	}

	if(strcmp(ParName,"stopping_model")==0){ /* ... */
		sscanf(ParValue,"%i",&stopping_model);
		message(1,"Stopping model is:\t\t %i\n",stopping_model);
	}

	if(strcmp(ParName,"flight_length_type")==0){ /* ... */
		sscanf(ParValue,"%i",&flight_length_type);
		message(1,"Flight length distribution model:%i\n",flight_length_type);
	}
	if(strcmp(ParName,"flight_length_constant")==0){ /* ... */
		sscanf(ParValue,"%f",&flight_length_constant);
		message(1,"Flight length:\t\t\t %g nm\n",flight_length_constant);
	}

	if(strcmp(ParName,"normalize_output")==0){ /* ... */
		sscanf(ParValue,"%i",&normalize_output);
		message(1,"Normalize Output :\t\t %i \n",normalize_output);
	}

	if(strcmp(ParName,"ion_dose")==0){ /* ... */
		sscanf(ParValue,"%lf",&ion_dose);
		message(1,"Dose:\t\t\t %g \n",ion_dose);
	}
	if(strcmp(ParName,"dpa_output")==0){ /* ... */
		sscanf(ParValue,"%i",&dpa_output);
		message(1,"dpa output :\t\t %i \n",dpa_output);
	}


	if(strcmp(ParName,"multiple_collisions")==0){    /* ... */
		sscanf(ParValue,"%i",&max_annular_coll_volumes);
		message(1,"Multiple collisions :\t\t %i \n",max_annular_coll_volumes);
	}

	if(strcmp(ParName,"scattering_calculation")==0){ /* ... */
		sscanf(ParValue,"%i",&scattering_calculation);
		message(1,"scattering_calculation:\t\t %i \n",scattering_calculation);
	}

	if(strcmp(ParName,"transport_type")==0){ /* ... */
		sscanf(ParValue,"%i",&transport_type);
		message(1,"transport_type:\t\t\t %i \n",transport_type);
	}

	if(strcmp(ParName,"single_ion_sputter_yields")==0){ /* ... */
		sscanf(ParValue,"%i",&single_ion_sputter_yields);
		message(1,"single_ion_sputter_yields:\t %i \n",single_ion_sputter_yields);
	}

	if(strcmp(ParName,"do_not_store_damage")==0){ /* ... */
		sscanf(ParValue,"%i",&do_not_store_damage);
		message(1,"Do not store damage:\t\t %i \n",do_not_store_damage);
	}
	if(strcmp(ParName,"no_headers_in_files")==0){ /* ... */
		sscanf(ParValue,"%i",&no_headers_in_files);
		message(1,"no headers in file:\t\t %i \n",no_headers_in_files);
	}


	return 0;
}

FILE* OpenFileContinuous(char* BaseName, char* Extension){
	/* Opens file basename+extension, keeps it open */
	int BaseNameLength;
	char* strTemp;
	BaseNameLength=strlen(BaseName);
	strTemp=(char*)malloc(sizeof(char)*(BaseNameLength+50));
	strncpy(strTemp,BaseName,BaseNameLength+1);
	strcat(strTemp,Extension);
	return fopen(strTemp,"wt");
}

int FloatBlockReader(char* FileName, int Offset, int Count, float* Array){
	/* Reads a block of Count float values from file FileName starting at Offset
	and puts these in Array */
	/* Might cause a seg fault or crash if array is too short! */

	FILE* fp;

	fp=fopen(FileName,"rb");
	if(fp==NULL){return -5013;  } /* Cannot open file */
	fseek(fp,sizeof(float)*Offset,1);
	fread((void*)Array, sizeof(float), Count,fp); /* read data block */
	fclose(fp);

	return 0;
}


int LoadInverseErf() {
	/* Function adapted from corteo.c: */
	/* load a list of inverse error function erfinv(x) values that contains MAXERFLIST elements
	the list is such that x is uniformly distributed between -1 and 1 excluding these boundaries
	i.e. x = -1+dx, -1+2dx, ...1-2dx, 1-dx, with 1/2dx = MAXERFLIST */
	/* return 0 on success, some other value else */ 

	unsigned int k;
	char erfval[1000];
	FILE * ifp;
	char dataFile[1000];

	sprintf(dataFile, "%s/erfinv.dat", DirectoryData); /* dir + erfinv.dat */
	message(2,"Open %s\n", dataFile);
	
	ifp = fopen(dataFile, "r");
	if(ifp==NULL){return -5014;}

	for(k=0; k<MAXERFLIST; k++) {
		fscanf(ifp,"%s", erfval);
		inverse_erf_list[k] = a2f(erfval);
	}
	fscanf(ifp,"%s", erfval);
	fclose(ifp);
	/* control value at the end should be the number of elements in the list */
	if(atoi(erfval)!=MAXERFLIST){
		return -5015;
	}

	return 0;
}

int WriteStringToFile(char* Filename, char* str){
	/* does what you think it does */
	/* returns 0 on success */
	FILE* fp;
	fp=fopen(Filename,"w");
	if(fp==NULL){
		return -1;
	}
	fprintf(fp,"%s",str);
	fclose(fp);
	return 0;
}

int display_a_file(char* Filename){
	/* prints the contents of a text file to the std out.
	Returns 0 on success. */
	FILE* fp;
	fp=fopen(Filename,"r");
	if(fp==NULL){
		return -1;
	}
	while(!feof(fp)){
		putchar(fgetc(fp));
	}
	return 0;
}

int SplitSingleInputFile(char* filename){
	/* if a single input file is provided, split it to temp files that can be read later */
	/* argument is the name of the combined input file */

	FILE *fp;        /* Pointer to combined input file */
	char* line;      /* To read current line */
	FILE* fout[4];   /* Pointer to the four single output files */
	int filenum;     /* file, into which is currently written */
	int sepline;     /* if the current line separates files, this is 1 */
	int i;
	char* tempFN;    /* to store filename */


	/* Make a copy of the configfilename, so we can change it later */
	tempFN=(char*)malloc(sizeof(char)*1024);
	if(tempFN==NULL){message_error(-5015,"Insufficient memory!\n");return -5015;}
	strncpy(tempFN,ConfigFileName,1023);

	strcpy(ConfigFileName,"temp_configfile.iradina");          /* Name of the general input config file */
	strcpy(TargetStructureFileName,"temp_structfile.iradina"); /* Filename of the file that define the structure of the target */
	strcpy(TargetCompositionFileName,"temp_compfile.iradina"); /* Target composition file */
	strcpy(MaterialsFileName,"temp_matfile.iradina");          /* Name of the file that defines the materials in the target */

	/* remove existing temp files: */
	remove(ConfigFileName);
	remove(TargetStructureFileName);
	remove(TargetCompositionFileName);
	remove(MaterialsFileName);
	single_input_file=1;   /*mark that this is used indeed */

	message(1,"Splitting combined input file...\n");

	/* open input file */
	fp=fopen(tempFN,"r");
	if(fp==NULL){message_error(-5016,"Cannot open combined config file: %s\n",filename);return -5016;}

	/* open output files */
	fout[0]=fopen(ConfigFileName,"w");
	if(fout[0]==NULL){message_error(-5017,"Cannot open file for writing: %s\n",ConfigFileName);return -5017;}
	fout[1]=fopen(TargetStructureFileName,"w");
	if(fout[1]==NULL){message_error(-5018,"Cannot open file for writing: %s\n",TargetStructureFileName);return -5018;}
	fout[2]=fopen(MaterialsFileName,"w");
	if(fout[2]==NULL){message_error(-5019,"Error: Cannot open file for writing: %s\n",MaterialsFileName);return -5019;}
	fout[3]=fopen(TargetCompositionFileName,"w");
	if(fout[3]==NULL){message_error(-5021,"Cannot open file for writing: %s\n",TargetCompositionFileName);return -5021;}

	line=(char*)malloc(sizeof(char)*512);   /* reserve space to read line */
	if(line==NULL){return -5022;}             /* out of mem */

	filenum=0;

	while(!feof(fp)){ /* until file ends */
		sepline=0;
		if(fgets(line,511,fp)!=NULL){ /* Now the line has been read */
			/* Check for file separators */
			if(strncmp(line,"#<<<BEGIN CONFIGFILE",20)==0){ /* it's the config file that follows */
				message(2,"Config file part found.\n");
				filenum=0;sepline=1;
			}
			if(strncmp(line,"#<<<BEGIN STRUCTUREFILE",23)==0){ /* it's the structures file that follows */
				message(2,"Structure file part found.\n");
				filenum=1;sepline=1;
			}
			if(strncmp(line,"#<<<BEGIN MATFILE",17)==0){ /* it's the material file that follows */
				message(2,"Material file part found.\n");
				filenum=2;sepline=1;
			}
			if(strncmp(line,"#<<<BEGIN COMPFILE",18)==0){ /* it's the composition file that follows */
				message(2,"Composition file part found.\n");
				filenum=3;sepline=1;
			}

			if(sepline==0){ /* its a data line, put it into current output file */
				/* fprintf(fout[filenum],line); */
				fputs(line,fout[filenum]);
			}
		}
	}

	/* Close files */
	fclose(fp);
	for(i=0;i<4;i++){fclose(fout[i]);}

	message(1,"Splitting finished successfully.\n");

	return 0;
}

int CheckSplitInputFile(char* filename){
	/* checks if the given config file is a single combined input file (incl. structure, materials, and composition) */
	/* if so, the file is split into four separate temp files for further "conventional" processing */
	FILE* fp;
	char line[512];
	int result;
	result=0;

	/* open file, read first line, check if it is marker for combined file */
	fp=fopen(filename,"r");
	if(fp==NULL){return -5023;};
	if(!feof(fp)){
		if(fgets(line,511,fp)!=NULL){ /* Read line. */
			/*printf("DEBUG: %s\n\n",line);*/
			if(strncmp(line,"#<<<BEGIN CONFIGFILE",20)==0){    /* its a combined file! */
				message(1,"Combined input file detected.\n");
				result=1;
			}
		}
	}
	fclose(fp);
	if(result==1){ /* its combined, split it! */
		SplitSingleInputFile(filename);
	}
	return result;
}

int CombineFiles(int count, ...){
	/* Combine a number of files into one. */
	/* The number of strings must be provided! */
	/* The first string must be the input filename */
	/* Then we have alternating:
	- strings, which are printed into the file
	- filenames, which are appended to the file
*/
	va_list ap;

	int i;
	FILE* fp_out;
	FILE* fp_in;
	char* filename; /* to store filename stuff */
	char* line;

	filename=(char*)malloc(MAX_FILENAME_LENGTH*sizeof(char));
	line=(char*)malloc(2048*sizeof(char));
	if((filename==NULL)||(line==NULL)){
		message_error(-5024,"insufficient memory!\n");
		return -5024;
	}

	va_start(ap, count);
	if(count<1){ /* no output file! exit */
		message_error(-5025,"too few arguments for file combination!\n");
		va_end(ap);
		return -5025;
	}

	strcpy(filename,va_arg(ap,char*)); /* Obtain filename from next argument */
	fp_out=fopen(filename,"w");        /* open output file */
	if(fp_out==NULL){
		message_error(-5026,"Cannot open combined output file %s.\n",filename);
		return -5026;
	}

	/* create combined output file */
	for(i=1;i<count;i++){
		if( (i%2) == 1){                 /* it's a string to include */
			strcpy(line,va_arg(ap,char*));       /* Obtain string from next argument */
			fprintf(fp_out,"%s\n",line);
		} else {                         /* it's a filename to include */
			strcpy(filename,va_arg(ap,char*)); /* Obtain filename from next argument */
			/* open file, read lines, copy them into output file */
			fp_in=fopen(filename,"r");
			if(fp_in==NULL){message_error(-5027,"Error: Cannot open file %s for reading.\n",filename);return -5027;}
			while(!feof(fp_in)){                /* until file ends */
				if(fgets(line,2048,fp_in)!=NULL){ /* Now the line has been read */
					fputs(line,fp_out);             /* put into output file */
				}
			}
			fclose(fp_in);
		}
	}
	fclose(fp_out);
	va_end(ap);
	return 0;
}

int WriteFileHeader(char* filename, char * title, char * ct4 , char * cn4){
	FILE* fp;


	fp=fopen(filename,"wt");
	fprintf(fp,"#NAME: %s \n",filename);
	fprintf(fp,"#TITLE: %s ", title);
	if((ion_dose==-1.0)||(ion_dose==1)){
		fprintf(fp," \n");
	}
	else {
		if(dpa_output==1){
			fprintf(fp,"dpa output") ;
		}
		fprintf(fp," iondose= %g \n",ion_dose) ;
	}
	fprintf(fp,"#DATE: %s", asctime(timeinfo) );
	fprintf(fp,"#COLUMN_NAMES: xc | yc | zc | %s \n",cn4);
	fprintf(fp,"#COLUMN_TITLES: Depth | y width | z width| %s \n",ct4 );
	fprintf(fp,"#COLUMN_UNITS: nm | nm | nm | ");
	if((ion_dose==-1.0)||(ion_dose==1)){
	  if(normalize_output>=1){
	    fprintf(fp,"$cm^{-1}$ per ion \n");
	  }
	  else {
	    fprintf(fp," \n");
	  }
	}
	else {
	  if(dpa_output==1){
	    fprintf(fp,"dpa\n") ;
	  }
	  else {
	    if(normalize_output>=1){
	      fprintf(fp,"$cm^{-1}$ \n");
	    }
	    else {
	      fprintf(fp," \n");
	    }
	  }
	}
	fprintf(fp,"#COLUMN_FACTOR: %g | %g | %g | 1.0 \n", cell_size_x, cell_size_y,cell_size_z);
	fprintf(fp," \n");

	fclose(fp);
	return 0;
}

int WriteEnergyFileHeader(char* filename, char * title, char * ct4 , char * cn4){
	FILE* fp;


	fp=fopen(filename,"wt");
	fprintf(fp,"#NAME: %s \n",filename);
	fprintf(fp,"#TITLE: %s ", title);
	if((ion_dose==-1.0)||(ion_dose==1)){
		fprintf(fp," \n") ;
	}
	else {
		fprintf(fp," iondose= %g \n",ion_dose) ;
	}
	fprintf(fp,"#DATE: %s", asctime(timeinfo) );
	fprintf(fp,"#COLUMN_NAMES: xc | yc | zc | %s \n",cn4);
	fprintf(fp,"#COLUMN_TITLES: Depth | y width | z width| %s \n",ct4 );
	fprintf(fp,"#COLUMN_UNITS: nm | nm | nm | ");
	if((ion_dose==-1.0)||(ion_dose==1)){
		if(normalize_output>=1){
			fprintf(fp,"eV.$cm^{-1}$ per ion \n");
		}
		else {
			fprintf(fp," eV \n");
		}
	}
	else {
		if(normalize_output>=1){
			fprintf(fp,"eV.$cm^{-1}$ \n");
		}
		else {
			fprintf(fp," eV \n");
		}

	}
	fprintf(fp,"#COLUMN_FACTOR: %g | %g | %g | 1.0 \n", cell_size_x, cell_size_y,cell_size_z);
	fprintf(fp," \n");

	fclose(fp);
	return 0;
}
int WriteResArrayToFile(char* Filename,int* SourceArray, int Count, int FileType, float conc){
  /* Writes the designated array of Count elements into a file. 
     If the file is just one column of values then FileType should be
     set to 1. If the file contains 4 columns (x, y, z, value)
     then set it to 0.
     The caller needs to make sure, that the array is large enough. */

  //	FILE* fp;
  //	int   i;
  //	int   x,y,z;
  float* dscaled2 ;
	
  if (dpa_output==0){
    WriteIntArrayToFile(Filename,SourceArray,Count,FileType);
  }
  else{
    dscaled2=(float*)malloc(sizeof(float) * cell_count); 
    DensityScaleArray(SourceArray,dscaled2, conc);
    WriteFloatArrayToFile(Filename,dscaled2,Count,FileType);
    free (dscaled2);
  }


  return 0;
}

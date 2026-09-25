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


#include "iradina.h"

int main(int argc, char* argv[]){
	int result;
	int cmd_result;

	message_buffer=(char*)malloc(4096*sizeof(char)); /* buffer for printout of messages */

	/* Init some default values for some global variables */
	print_level            = 0;
	wait_before_end        = 0;
	mem_usage_only         = 0;
	mem_usage_details      = 0;
	mem_usage              = 0;
	override_max_ions      = 0;
	override_energy        = 0;
	display_interval       = 100;
	store_path_limit       = -1;
	store_path_limit_recoils = -1;
	special_geometry       = 0;
	store_exiting_recoils  = 0;
	store_exiting_limit    = 10;
	store_info_file        = 1;
	store_joined_output    = 0;
	max_annular_coll_volumes=0;
	scattering_calculation = 0;
	transport_type         = 0;
	normalize_output       = -1; /* -1: not defined by user*/
	ion_dose               = -1.0;
	do_not_store_damage    = 0;
	unit_conversion_factor = 1.0;
	create_status_file     = 0;
	single_input_file      = 0;
	hydrogen_in_target     = 0;
	status_update_interval = 100000;
	conv_create_separate_elements = 0;

	/* Obtain memory for filenames */
	ConfigFileName=(char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(ConfigFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	strcpy(ConfigFileName,"Config.in"); /* Default */
	DirectoryData=(char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(DirectoryData==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	strcpy(DirectoryData,"./data"); /* Default as (previously) data dir in current dir*/
	MaterialsFileName=(char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(MaterialsFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	strcpy(MaterialsFileName,"Materials.in"); /* Default */
	ConversionFileName=(char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(ConversionFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	strcpy(ConversionFileName,"Converted.input"); /* Default */
	ElementsFileName=(char*)malloc(sizeof(char)*MAX_FILENAME_LENGTH);
	if(ElementsFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	strcpy(ElementsFileName,"Elements.in"); /* Default */
	
	TargetStructureFileName=(char*)malloc(sizeof(char)*1024);
	if(TargetStructureFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}
	TargetCompositionFileName=(char*)malloc(sizeof(char)*1024);
	if(TargetCompositionFileName==NULL){message_error(-1,"Not enough memory\n"); return -1;}

	/* Handle command line arguments: */
	cmd_result=handle_cmd_line_options(argc,argv);
	if(cmd_result==0){ /* Normal simulation */
	  printf("%i\n",simulation_type);
		display_startup_message();
		/* check for correct type representation on machine: */
		if( (result=check_type_representation()) !=0){
			message_error(-107,"incorrect type on machine (%i).\n",result);return -107;
		}
		if(create_status_file==1){write_status_file("init0", 0);}
		result=InitConfiguration(ConfigFileName);   /* init all things */
		printf("%i\n",simulation_type);
		if(result!=0){
			message_error(result,"ERROR: Initialization failed (%i).\n",result);
			return result;
		}
                if (simulation_type==7){
		  minE=100. ;
		  print_stopping_table(ionZ,ionM,stopping_target_index,minE,ionInitialEnergy,minE);
		  store_ion_paths=0;
		  store_recoil_cascades=0;
		  store_range3d=0;
		  store_PKA=0;
		  ion_distribution=4;
		  //		  max_no_ions=100000;
		  do_not_store_damage=1;
		  printf("#Printing ballistic stopping table \n");
		  printf("#Energy/eV\tStopping/(eV/nm)\n");
		  ball_stopping_fp=OpenFileContinuous(OutputFileBaseName,".ball_loss");
		  if(ball_stopping_fp==NULL){
		    message_error(-1002,"Cannot open range3d-file!\n");
		    return -1002;
		  }
		  if (no_headers_in_files==0) {
		    fprintf(ball_stopping_fp,"#NAME: %s \n", OutputFileBaseName,".ball_loss");
		    fprintf(ball_stopping_fp,"#TITLE: Ballistic Stopping Power \n");
		    fprintf(ball_stopping_fp,"#DATE: %s", asctime(timeinfo) );
		    fprintf(ball_stopping_fp,"#COLUMN_NAMES: E | N1 | N2 | Eloss  \n");
		    fprintf(ball_stopping_fp,"#COLUMN_TITLES: Energy | N3 | N4| Stopping power \n");
		    fprintf(ball_stopping_fp,"#COLUMN_UNITS: eV | N5 | N6 | ev/nm \n");
		    fprintf(ball_stopping_fp,"#COLUMN_FACTOR: %i | %i | %i | 1.0 \n", 1,1,1);
		    fprintf(ball_stopping_fp," \n");
		    //	fclose(elecfp);
		  }
		}
	 
		if(create_status_file==1){write_status_file("init1", 0);}
		if(mem_usage_only==1){ /* Do not simulate, only show memory usage */
		message(-1,"Estimated total memory requirement: %li MByte\n",(mem_usage/0x100000)+1);
		/* Estimate about 1 MB overhead (for arrays of pointers, recursive function calls and so on */
		return 99;
		}
		if(override_max_ions>0){ /* the command line option has been used to override ion number */
		message(1,"Maximum number of ions limited to %i.\n",override_max_ions);
		max_no_ions=override_max_ions;
		}
		if(override_energy>0.0000001){ /* the command line option has been used to override ion energy */
		message(0,"New ion energy: %g.\n",override_energy);
		ionInitialEnergy=override_energy;
		}
		if(create_status_file==1){write_status_file("init2", 0);}
		message(1,"\nEverything is initialized and ready.\n-----------------------------------------------\n\n");
		message(-1,"Starting simulation of irradiation...\n");
		calculate_normalization_factor(max_no_ions);
		message(1,"Normalization factor is: %g\n",unit_conversion_factor);
		fflush(stdout);
		/* Irradiate the target with ions: */
		IrradiateTarget();
		message(1,"\nSimulation is finished.\n-----------------------------------------------\n");
		message(-1,"Storing final results: ...\n ");fflush(stdout);
		if(create_status_file==1){write_status_file("simend", 0);}
		calculate_normalization_factor(max_no_ions);
		store_results(OutputFileBaseName,max_no_ions);
		message(-1,"done.\n\n");
		if(wait_before_end==1){printf("Press enter to exit.\n");getc(stdin);}
		if(create_status_file==1){write_status_file("end", 0);}
		
	} else {  /* Do something else than simulate */
		if(cmd_result<0){ /* An error occured */
		return cmd_result;
		} else { /* Perform some operation that is not simulation */
		switch(cmd_result){
			case 3: /* convert material to element file */
			printf("\nYou are using iradina to convert a material based target definition to an element based one.\n\n");
			result=InitConfiguration(ConfigFileName);   /* init all things, needs to be done to read material file */
			if(result!=0){
				message_error(result,"Initialization error (%i). Aborting.\n",result);
				return result;
			}
			/* Call converter: */
			result=MaterialToElementConverter(ConversionFileName);
			if(result!=0){
				message_error(result,"ERROR: Conversion failed (%i).\n",result);
			return result;
			}
				break;
			case 1: break;	/* help was printed. Do nothing else */
			case 4: /* print stopping table only. Initialization is still required to read target and ion properties! */
				/* check for correct type representation on machine: */
				if( (result=check_type_representation()) !=0){printf("Error 107. Type checker returned %i\n",result);return -107;}
					result=InitConfiguration(ConfigFileName);   /* init all things */
				if(result!=0){
					message_error(result,"Init error: %i. Aborting.\n",result);
				return result;
			}
			/* Calculate and print stopping tables: */
			minE=ionInitialEnergy/100. ;
			print_stopping_table(ionZ,ionM,stopping_target_index,minE,ionInitialEnergy,minE);
			break;
			default:
				message_error(0,"Unknown operation to perform.\n");
			}
		}
	}
	return 0;
}
	
	
int display_startup_message(){
	 /* http://patorjk.com/software/taag/#p=display&h=2&f=Big&t=iradina */
	char *title = "\n  _               _ _             "
		            "\n (_)             | (_)            "
                "\n  _ _ __ __ _  __| |_ _ __   __ _ "
                "\n | | '__/ _` |/ _` | | '_ \\ / _` |"
                "\n | | | | (_| | (_| | | | | | (_| |"
                "\n |_|_|  \\__,_|\\__,_|_|_| |_|\\__,_|";

	time_t rawtime;
//	struct tm * timeinfo;
	FILE* fpout;
	fpout=stdout;
	if(print_level>=-1){
		fprintf(fpout,title);

		fprintf(fpout,"\n\nThis is iradina version %i.%i.%i %s, ",VERSION,SUBVERSION,SUBSUBVERSION,RELEASESTRING);
		fprintf(fpout,"%s\n%s\n",VERSIONDATE,VERSIONCOMMENT);
		fprintf(fpout,"by C. Borschel, 2019\nInstitute for Solid State Physics, University of Jena\n\n");
		fprintf(fpout,"This program comes with ABSOLUTELY NO WARRANTY.\n");
		fprintf(fpout,"This is free software, and you are welcome to redistribute it under certain conditions,\n");
		fprintf(fpout,"see license.txt for details.\n\n");
		fprintf(fpout,"KP additions by JP.Crocombette (cea).\n");
		char *p=getenv("USER");
		fprintf(fpout,"Compiled code on %s at %s by %s.\n\n",__DATE__,__TIME__,p);

	#ifdef INCLUDE_SPECIAL_GEOMETRY
		fprintf(fpout,"This is the %s version of iradina.\n\n",SPECIAL_GEOMETRY_NAME);
	#endif    
	}

	time( &rawtime );
	timeinfo = localtime( &rawtime );
	message(0,"Current time: %s\n", asctime(timeinfo) );
	return 0;
}

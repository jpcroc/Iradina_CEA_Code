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

#include "target.h"


/* Error codes used in this modul: between 550 and 700                       */

int InitializeMaterials(char* Filename){
  /* Read materials from input file and do preparatory calculations, etc.*/
  /* It must be called before calling InitializeTarget. */
  /* returns 0 on success */

  int   i,j=0,k=0;
  float temp;
  int   result;
  int   NumExEl;     /* Number of existing elements */
  unsigned long int lui_temp;
  struct scattering_matrix *psm; /* Pointer to a scattering matrix */

  /* Init some material arrays etc. */
  NumberOfMaterials=0;
  for(i=0;i<MAX_NO_MATERIALS;i++){ /* set some default values for each material: */
    sprintf(ListOfMaterials[i].Name,"new material %i",i); /* some initial name */
    for(j=0;j<MAX_EL_PER_MAT;j++){ /* some default values for each element: */
      ListOfMaterials[i].ElementsReplEnergy[j]=-1.0; /* use -1.0 as dummy value (don't use 0 in case user wants to use 0) */
    }
  }


  /* Read in file with material info */
  result=IniFileReader(Materials_DataBlockReader, Materials_DataReader, Filename);
  if(result!=0){
    return result;
  }

  /* Fill existing elements with 0 */
  for(i=0;i<MAX_ELEMENT_NO;i++){
    existing_elements[i]=0;
  }

  /* Loop through materials, for each:
     - normalize elemental concentrations to 1,
     - calculate mean interatomic distance
     - calculate maximum impact parameter (not reduced one).
     - determine elements existing in target
  */
  number_of_all_target_elements=0;
  for(i=0;i<NumberOfMaterials;i++){
    temp=0; 
    ListOfMaterials[i].CellCount=0;
    for(k=0;k<(MAX_EL_PER_MAT*6);k++){ /* Reset sputtercounters: */
      ListOfMaterials[i].SputterCounter[k]=0;
    }
    for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* Sum up. Also indicate that this element exists in the target. */
      temp+=ListOfMaterials[i].ElementsConc[j];
      existing_elements[ListOfMaterials[i].ElementsZ[j]]=1;
      number_of_all_target_elements++;
    }

    /* If ion SBE is not given, calculate mean value: */
    if( ListOfMaterials[i].IonSurfEnergy<0){
      ListOfMaterials[i].IonSurfEnergy=0;
      for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* Sum up */
	ListOfMaterials[i].IonSurfEnergy+=ListOfMaterials[i].ElementsSurfEnergy[j];
      }
      ListOfMaterials[i].IonSurfEnergy/=ListOfMaterials[i].ElementCount;
      message(0,"No ion SBE for material %i defined. Trying mean value: %f eV\n",i,ListOfMaterials[i].IonSurfEnergy);
    }
    if(isnan(ListOfMaterials[i].IonSurfEnergy)){
      ListOfMaterials[i].IonSurfEnergy=3;
      message(0,"No ion SBE for material %i defined. Using %f eV\n",i,ListOfMaterials[i].IonSurfEnergy);
    }
    ListOfMaterials[i].MeanZ=0;
    ListOfMaterials[i].MeanM=0;
    for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* go through elements, normalize concentration, check for hydorgen existence*/
      ListOfMaterials[i].ElementsConc[j]=ListOfMaterials[i].ElementsConc[j]/temp;
      /* calculate average M,Z */
      ListOfMaterials[i].MeanZ+=ListOfMaterials[i].ElementsConc[j]*ListOfMaterials[i].ElementsZ[j];
      ListOfMaterials[i].MeanM+=ListOfMaterials[i].ElementsConc[j]*ListOfMaterials[i].ElementsM[j];
      if(ListOfMaterials[i].ElementsZ[j]==1){hydrogen_in_target=1;} /* Hydrogen is really within target! */
      if(ListOfMaterials[i].ElementsZ[j]==ionZ){ionZ_in_target=1;} /* An element like the ion is within target! */
      /* in case replacement threshold is not given, set it to e_lattice as default: */
      if(ListOfMaterials[i].ElementsReplEnergy[j]==-1.0){
	if(simulation_type==1){
	  ListOfMaterials[i].ElementsReplEnergy[j]=ListOfMaterials[i].ElementsDispEnergy[j];
	  message(2,"No replacement threshold defined for material %i, element %i. Use DISPLACEMENT energy as default: %f eV\n",i,j,ListOfMaterials[i].ElementsReplEnergy[j]);
	}
	if(simulation_type==0){
	  ListOfMaterials[i].ElementsReplEnergy[j]=ListOfMaterials[i].ElementsLattEnergy[j];
	  message(2,"No replacement threshold defined for material %i, element %i. Use lattice energy as default: %f eV\n",i,j,ListOfMaterials[i].ElementsReplEnergy[j]);
	}				
      }

    }

    /* For testing purposes and comparisons for SRIM, we calculate the mean screening length and energy reduction factor (as in ZBL85) */
    ListOfMaterials[i].MeanA= d2f( 0.1f * SCREENCONST / ( pow(ionZ,0.23) + pow(ListOfMaterials[i].MeanZ,0.23) )  );
    ListOfMaterials[i].MeanF= 10.0f * ListOfMaterials[i].MeanA * ListOfMaterials[i].MeanM / ( ionZ * ListOfMaterials[i].MeanZ * (ionM+ListOfMaterials[i].MeanM) * E2 );

    ListOfMaterials[i].MeanMinRedTransfer= min_energy*ListOfMaterials[i].MeanF * (ionM+ListOfMaterials[i].MeanM)*(ionM+ListOfMaterials[i].MeanM) / (4*ionM*ListOfMaterials[i].MeanM) ;
		
    ListOfMaterials[i].DensityNM          = ListOfMaterials[i].Density*1e-21; /* Calculate density also in at/nm^3 */
    ListOfMaterials[i].AtomicDistance     = 1.0/pow( 4.0*PI*ListOfMaterials[i].Density/3.0, 1.0/3.0)*1e7 ; /* for conversion from cm to nm */
    ListOfMaterials[i].LayerDistance      = 1.0/pow( ListOfMaterials[i].Density, 1.0/3.0)*1e7 ; /* for conversion from cm to nm */
    ListOfMaterials[i].SqrtAtomicDistance = sqrtdf(ListOfMaterials[i].AtomicDistance);
    ListOfMaterials[i].SqrtRecFlDensity   = 1.0 / sqrt(PI * flight_length_constant * ListOfMaterials[i].DensityNM);
    ListOfMaterials[i].MeanImpactPar = 1.0f / sqrtdf(PI*ListOfMaterials[i].DensityNM*ListOfMaterials[i].AtomicDistance);
    /* Later, we need to multiply this by:
       1. sqrt of a random number between 0 and 1     "r2" in corteo manual, select impact par randomly
       2. 1/ln(-random)                               "r1" selects flight length according to poisson statistics
       3. 1/screening_length                           to get reduced quantity
       in order to get the actual reduced impact parameter */
#ifdef DEBUG_MODE
    printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);
    printf("i: %i, MeIP: %g, Product: %g\n,",i,ListOfMaterials[i].MeanImpactPar,PI*ListOfMaterials[i].DensityNM*ListOfMaterials[i].AtomicDistance);
#endif


  }

  //  switch(stopping_model
  int flag=0;
  existing_elements[ionZ]=1; /* the ion occurs as a possible projectile of course */
  /* TO DO: Can we reduce this? For the ion we only need this for stopping... not for the ScatMatrix */
  if (existing_elements[1]==0){
    /* hydrogen stopping is needed to calculate charge state */
    existing_elements[1]=1;    /* hydrogen stopping is needed to calculate charge state */
    flag=1;
  }
  /* TO DO: keep different tables of existing elements for stopping and scattering, to save memory */

  if(mem_usage_only==0){
    /* Load stopping tables for electronic stopping for all possible elements in all materials */
    //    printf("STOPPINGMODEL: %d\n", stopping_model);
    result=PrepareStoppingTables();
    switch (stopping_model) {
        case 0:

	  break;
      case 1:
	result=PrepareStoppingTablesdpass();
	break;
      case 2:
      if (flag==1){
      existing_elements[1]=0;    /* hydrogen stopping is needed to calculate charge state */
      }
      result=PrepareStoppingTablesfile();
      if(result!=0){
	message_error(result,"Cannot create FILE  stopping tables!\n");
	return result;
      }
      existing_elements[1]=1;    /* hydrogen stopping is needed to calculate charge state */
      break;
    

      default:
	printf("wong electronic stopping_model\n");
	exit(1);
  }


	//    result=PrepareStoppingTables();
    if(result!=0){
      switch(result){
      case -1:
	message_error(result,"Not enough memory for stopping tables! (%i)\n",result);
	break;
      case -2:
	message_error(result,"Cannot read data from stopping file (%i)!\n",result);
	break;
      case -3:
	message_error(result,"control data mismatch in stopping file %i\n",j);
	break;
      default:
	message_error(result,"Cannot load stopping tables! (%i)\n",result);
      }
      return result;
    }
    //    result=PrepareStoppingTablesdpass();

    existing_elements[1]=1;    /* hydrogen stopping is needed to calculate charge state */
    
    /* Load straggling tables for electronic energy loss straggling for all possible elements in all materials */
    result=PrepareStragglingTables(straggling_model);
    if(result!=0){
      message_error(result,"Cannot create straggling tables!\n");
      return result;
    }
		
    /* Go through possible target element and create scattering matrices for scattering of ion from each target element */
    message(0,"Prepare scattering matrices for ion on target collisions...");fflush(stdout);
    for(i=0;i<MAX_ELEMENT_NO;i++){
      if(existing_elements[i]==1){ /* Possible target element */
	psm=&(ion_scattering_matrix[i]);
	result=PrepareScatteringMatrix(psm,ionM,ionZ,MostAbundantIsotope[i],i);
	/* Note! This scattering is only valid for the most abundant isotope!
	   If scattering for different isotopes was accounted this would be a bit heavy on memory.
	   If more accuracy is need, we should use MAGIC instead for each collision event */
	if(result!=0){
	  message_error(result,"insufficient memory for ion-target scattering matrices! %i \n",result);
	  return result;
	}
      }
    }
    message(0," finished\n");fflush(stdout);
		
    /* Go through all possible other projectiles (which are the target elements) and create scattering matrices
       for scattering of the projectile with any other target element */
    /* Beware, this twodimensional matrix might get rather large... */
    message(0,"Prepare scattering matrices for recoil on target collisions...");fflush(stdout);
    for(j=0;j<MAX_ELEMENT_NO;j++){
      if(existing_elements[j]==1){ /* Possible projectile element */
	for(i=0;i<MAX_ELEMENT_NO;i++){
	  if(existing_elements[i]==1){ /* Possible target element */
	    result=PrepareScatteringMatrix(&(scattering_matrices[j][i]),MostAbundantIsotope[j],j,MostAbundantIsotope[i],i);
	    /* Not: First index in scattering matrix is projectile, second index is target */
	    if(result!=0){
	      message_error(result,"insufficient memory for target-target scattering matrices! %i \n",result);
	      return result;
	    }
	  }
	}
      }
    }

    /* For each element, possibly create arrays for storing leaving recoils */
    if(store_exiting_recoils==1){
      /* Create arrays first: */
      for(i=0;i<NumberOfMaterials;i++){
	if(ListOfMaterials[i].Is_Vacuum==0){ /* do not store stuff for vacuum */
	  ListOfMaterials[i].ElementalLeavingRecoils = malloc(sizeof(struct transmitted_ion*) * ListOfMaterials[i].ElementCount);
	  ListOfMaterials[i].leaving_recoils_pointer = malloc(sizeof(int) * ListOfMaterials[i].ElementCount);
	  for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* Go through elements */
	    (ListOfMaterials[i].ElementalLeavingRecoils)[j] = malloc(store_exiting_limit * sizeof(struct transmitted_ion));
	    if( (ListOfMaterials[i].ElementalLeavingRecoils)[j] == NULL ){
	      message_error(-626,"Insufficient memory for arrays of leaving recoils!\n");
	      return -626;
	    }
	    (ListOfMaterials[i].leaving_recoils_pointer)[j]=0;
	  }
	}
      }
    }


    message(0," finished\n");fflush(stdout);
  } else { /* only estimate memory usage */

    if(store_exiting_recoils==1){
      lui_temp=0;
      for(i=0;i<NumberOfMaterials;i++){
	if(ListOfMaterials[i].Is_Vacuum==0){ /* do not store stuff for vacuum */
	  lui_temp+=(sizeof(struct transmitted_ion*)+sizeof(int)+store_exiting_limit*sizeof(struct transmitted_ion))*ListOfMaterials[i].ElementCount;
	}
      }
      mem_usage+=lui_temp;
      if(mem_usage_details==1){printf("MEMORY Leaving recoils arrays:      %lu bytes\n",lui_temp);}
    }
    NumExEl=count_existing_elements(existing_elements);
    lui_temp=NumExEl*(NumberOfMaterials-1)*MAX_STOPPING_ENTRIES*sizeof(float);
    mem_usage+=lui_temp;
    NumExEl=count_existing_elements(existing_elements);
    lui_temp=NumExEl*(NumberOfMaterials-1)*MAX_STOPPING_ENTRIES_dpass*sizeof(float);
    mem_usage+=lui_temp;
    if(mem_usage_details==1){printf("MEMORY Stopping tables:             %lu bytes\n",lui_temp);}
    mem_usage+=lui_temp;
    if(mem_usage_details==1){printf("MEMORY Straggling tables:           %lu bytes\n",lui_temp);}
    lui_temp+=NumExEl*2*(DIME*DIMS*sizeof(float));
    mem_usage+=lui_temp;
    if(mem_usage_details==1){printf("MEMORY Ion scattering matrices:     %li bytes\n",lui_temp);}
    lui_temp+=NumExEl*NumExEl*2*(DIME*DIMS*sizeof(float));
    mem_usage+=lui_temp;
    if(mem_usage_details==1){printf("MEMORY Recoil scattering matrices:  %li bytes\n",lui_temp);}
  }
  return 0;
}

int Materials_DataBlockReader(char* MaterialName){
  /* Needs to be called from the ini file reader
     while the materials input file is read and the definition
     of a new material begins */

  /* Check, if there is still space left in the material list */
  if(NumberOfMaterials>=MAX_NO_MATERIALS){/* Material index is full */
    message_error(-627,"Warning! Maximum number of materials (%i) exceeded!\n",MAX_NO_MATERIALS);
    return -627;
  }

  /* Make sure, that materialname is not too long, i.e. cut it if it is */
  if(strlen(MaterialName)>=25){MaterialName[24]='\0';}

  strcpy(ListOfMaterials[NumberOfMaterials].Name,MaterialName); /* Store Name */
  ListOfMaterials[NumberOfMaterials].ElementCount=0;
  ListOfMaterials[NumberOfMaterials].Is_Vacuum=0;
  ListOfMaterials[NumberOfMaterials].Density=0;
  ListOfMaterials[NumberOfMaterials].IonSurfEnergy=-0.001;

  message(1,"\nNew material declared:\t%s\n",MaterialName);


  NumberOfMaterials++; /* increase index */

  return 0;
}

int Materials_DataReader(char* ParName, char* ParValue){
  /* Needs to be called from the ini file reader
     while the materials input file is read */

  /* Check some limits */
  if(NumberOfMaterials>=MAX_NO_MATERIALS){ /* Too large */
    return -628;
  }

  /* Compare the parameter name to known parameters and then read in corresponding value */

  if(strcmp(ParName,"ElementCount")==0){ /* Read number of element into the current material */
    sscanf(ParValue,"%i",&(ListOfMaterials[NumberOfMaterials-1].ElementCount));
    message(1,"Number of elements:\t%i\n",ListOfMaterials[NumberOfMaterials-1].ElementCount);
  }
  if(strcmp(ParName,"IsVacuum")==0){ /*  */
    sscanf(ParValue,"%i",&(ListOfMaterials[NumberOfMaterials-1].Is_Vacuum));
  }
  if(strcmp(ParName,"Density")==0){ /*  */
    sscanf(ParValue,"%f",&(ListOfMaterials[NumberOfMaterials-1].Density));
    if(isnan(ListOfMaterials[NumberOfMaterials-1].Density)){ListOfMaterials[NumberOfMaterials-1].Density=0.0;}
    message(1,"Density:\t\t%g\n",ListOfMaterials[NumberOfMaterials-1].Density);
  }
  if(strcmp(ParName,"ElementsZ")==0){ /* Read list of integer Z values*/
    if(make_int_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsZ)!=0){return -2;}
  }
  if(strcmp(ParName,"ElementsM")==0){ /* Read list of mass values*/
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsM)!=0){return -3;}
  }
  if(strcmp(ParName,"ElementsConc")==0){ /* Read list of concentrations */
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsConc)!=0){return -4;}
  }
  if(strcmp(ParName,"ElementsDispEnergy")==0){ /* Read list of displacement energies */
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsDispEnergy)!=0){return -5;}
  }
  if(strcmp(ParName,"ElementsLattEnergy")==0){ /* Read list of lattice energies */
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsLattEnergy)!=0){return -6;}
  }
  if(strcmp(ParName,"ElementsSurfEnergy")==0){ /* Read list of surface binding energies */
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsSurfEnergy)!=0){return -7;}
  }
  if(strcmp(ParName,"ElementsReplEnergy")==0){ /* Read list of replacement threshold energies */
    if(make_float_array(ParValue,MAX_EL_PER_MAT,ListOfMaterials[NumberOfMaterials-1].ElementsReplEnergy)!=0){return -8;}
  }
  if(strcmp(ParName,"IonSurfEnergy")==0){ /* Read ion surface binding energy */
    sscanf(ParValue,"%f",&(ListOfMaterials[NumberOfMaterials-1].IonSurfEnergy));
    message(1,"Ion SBE:\t\t%g\n",ListOfMaterials[NumberOfMaterials-1].IonSurfEnergy);
  }
  return 0;
}


int InitializeTargetStructure(char* Filename){
  /* Read target structure (size etc.) from config file
     and reads target concentration array etc.*/
  /* first, we need to read the config file about how large the target is,
     how many cells and so on... . Then we can initialze the target
     arrays, and then finally read in the composition matrix */

  int result;
  int i,j;
  unsigned long int lui_temp=0;

  message(0,"Initializing target.\n");

  /* Set some default values */
  UseDensityMult=0;

  /* Read structure definition file */
  result=IniFileReader(TargetStructure_DataBlockReader, TargetStructure_DataReader, Filename);
  if(result!=0){
    return result;
  }
  message(0,"Target structure definition file: %s\n",Filename);

  layer_count_yz=cell_count_y*cell_count_z;
  cell_count=layer_count_yz*cell_count_x;
  target_size_x = cell_count_x*cell_size_x;
  target_size_y = cell_count_y*cell_size_y;
  target_size_z = cell_count_z*cell_size_z;
  get_float_one_bit_smaller(&target_size_x,&target_max_x); /* we sometimes need a float that is one bit smaller than */
  get_float_one_bit_smaller(&target_size_y,&target_max_y); /* the targetsize itself */
  get_float_one_bit_smaller(&target_size_z,&target_max_z);
  /* Debug: print values for testing purposes: 
     printf("tmxyz: %f %f %f\n",target_max_x,target_max_y,target_max_z);
     printf("tmxyz: %x %x %x\n",*((int*)(&(target_max_x))),*((int*)(&(target_max_y))),*((int*)(&(target_max_z))));
     printf("txyz:  %x %x %x\n",*((int*)(&(target_size_x))),*((int*)(&(target_size_y))),*((int*)(&(target_size_z)))); */
  cell_volume   = cell_size_x * cell_size_y * cell_size_z;
  message(0,"\nTarget size is: \n");
  message(0,"x: %i cells, %g nm per cell, %g nm in total.\n",cell_count_x,cell_size_x,target_size_x);
  message(0,"y: %i cells, %g nm per cell, %g nm in total.\n",cell_count_y,cell_size_y,target_size_y);
  message(0,"z: %i cells, %g nm per cell, %g nm in total.\n",cell_count_z,cell_size_z,target_size_z);
  message(0,"Total: %i cells in %g nm^3.\n",cell_count,target_size_x*target_size_y*target_size_z);

  /* Now that we know the size of the target we can initialze some arrays etc. */
  if(mem_usage_only==0){
    TargetComposition=(int*)calloc(cell_count,sizeof(int));     /* Use calloc instead of malloc for initializing to zeros */
    TargetDensityMult=(float*)malloc(sizeof(float) * cell_count);
    if(TargetComposition==NULL){return -600;} /* Cannot allocate memory */
    if(TargetDensityMult==NULL){return -601;} /* Cannot allocate memory */
    TargetImplantedIons=(int*)calloc(cell_count,sizeof(int));
    TargetReplacingIons=(int*)calloc(cell_count,sizeof(int));
    TargetTotalDisplacements=(int*)calloc(cell_count,sizeof(int));
    TargetTotalInterstitials=(int*)calloc(cell_count, sizeof(int));
    TargetTotalVacancies=(int*)calloc(cell_count, sizeof(int));
    TargetTotalReplacements=(int*)calloc(cell_count, sizeof(int));
    if(TargetImplantedIons==NULL){return -602;} /* Cannot allocate memory */
    if(TargetReplacingIons==NULL){return -603;} /* Cannot allocate memory */
    if(TargetTotalDisplacements==NULL){return -604;} /* Cannot allocate memory */
    if(TargetTotalInterstitials==NULL){return -605;} /* Cannot allocate memory */
    if(TargetTotalVacancies==NULL){return -606;} /* Cannot allocate memory */
    if(TargetTotalReplacements==NULL){return -607;} /* Cannot allocate memory */
    if(detailed_sputtering==1){ /* For detailed calucation of sputtering, we need these */
      TargetVacuumNeighbors=(char*)calloc(cell_count,sizeof(char));
      TargetTotalSputtered=(int*)calloc(cell_count,sizeof(int));
      if(TargetVacuumNeighbors==NULL){return -608;} /* Cannot allocate */ 
      if(TargetTotalSputtered==NULL){return -609;} /* Cannot allocate */ 
    }
    if(store_energy_deposit==1){ /* For detailed storing of deposited energy, we need these */
      TargetEnergyPhonons=(double*)calloc(cell_count,sizeof(double));
      TargetEnergyElectrons=(double*)calloc(cell_count,sizeof(double));
      if(TargetEnergyPhonons==NULL){return -610;} /* Cannot allocate */ 
      if(TargetEnergyElectrons==NULL){return -611;} /* Cannot allocate */ 
    }
		
    /* Go through materials, init arrays for interstitials and vacancies */
    for(i=0;i<NumberOfMaterials;i++){
      /* Create arrays of pointers to arrays */
      ListOfMaterials[i].TargetImplantedRecoilsInt=(int**)malloc(sizeof(int*)*ListOfMaterials[i].ElementCount);
      ListOfMaterials[i].TargetImplantedRecoilsRepl=(int**)malloc(sizeof(int*)*ListOfMaterials[i].ElementCount);
      ListOfMaterials[i].TargetElementalVacancies=(int**)malloc(sizeof(int*)*ListOfMaterials[i].ElementCount);
      ListOfMaterials[i].TargetElementalDisp=(int**)malloc(sizeof(int*)*ListOfMaterials[i].ElementCount);
      if(ListOfMaterials[i].TargetImplantedRecoilsInt==NULL){return -612;}
      if(ListOfMaterials[i].TargetImplantedRecoilsRepl==NULL){return -613;}
      if(ListOfMaterials[i].TargetElementalVacancies==NULL){return -614;}
      if(ListOfMaterials[i].TargetElementalDisp==NULL){return -615;}
      if(detailed_sputtering==1){ /* For detailed calucation of sputtering, we need these */
	ListOfMaterials[i].TargetSputteredAtoms=(int**)malloc(sizeof(int*)*ListOfMaterials[i].ElementCount);
	if(ListOfMaterials[i].TargetSputteredAtoms==NULL){return -616;}
      }
      for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* Go through elements, make arrays for each one */
	ListOfMaterials[i].TargetImplantedRecoilsInt[j]=(int*)calloc(cell_count,sizeof(int));
	ListOfMaterials[i].TargetImplantedRecoilsRepl[j]=(int*)calloc(cell_count,sizeof(int));
	ListOfMaterials[i].TargetElementalVacancies[j]=(int*)calloc(cell_count,sizeof(int));
	/* printf("DEBUG target.c line %i, Mat %i, Elem %i, VacPointer %i\n",__LINE__,i,j,(int)ListOfMaterials[i].TargetElementalVacancies[j]); */
	ListOfMaterials[i].TargetElementalDisp[j]=(int*)calloc(cell_count,sizeof(int));
	if(ListOfMaterials[i].TargetImplantedRecoilsInt[j]==NULL){return -617;}
	if(ListOfMaterials[i].TargetImplantedRecoilsRepl[j]==NULL){return -618;}
	if(ListOfMaterials[i].TargetElementalVacancies[j]==NULL){return -619;}
	if(ListOfMaterials[i].TargetElementalDisp[j]==NULL){return -620;}
	if(detailed_sputtering==1){ /* For detailed calucation of sputtering, we need these */
	  ListOfMaterials[i].TargetSputteredAtoms[j]=(int*)calloc(cell_count,sizeof(int));
	  if(ListOfMaterials[i].TargetSputteredAtoms[j]==NULL){return -621;}
	}
      }
    }

    message(1,"Memory for target arrays has been allocated.\n\n");
		
    /* We can now read in the target composition file */
    /* Target composition file for non-dyanmic version based on materials */
    result=ReadIntFileIntoArray(TargetCompositionFileName,TargetComposition, cell_count, TargetCompositionFileType);
    if(result!=0){
      message_error(result,"cannot read target composition file %s (error %i).\n",TargetCompositionFileName,result);
      return result;
    }
    message(0,"Target composition read from %s.\n",TargetCompositionFileName);
#ifdef DEBUG_MODE4
    printf("%s, l: %i: DEBUG. cc: %i\n",__FILE__,__LINE__,cell_count);fflush(stdout);
    printf("%s, l: %i: DEBUG. No. of Mat: %i\n",__FILE__,__LINE__,NumberOfMaterials);fflush(stdout);
#endif
    /* For each material, count how many cells consist of this. */
    for(i=0;i<cell_count;i++){
#ifdef DEBUG_MODE4
      //printf("%s, l: %i: DEBUG. Cell %i has mat: %i\n",__FILE__,__LINE__,i,TargetComposition[i]);fflush(stdout);
#endif
      ListOfMaterials[TargetComposition[i]].CellCount+=1;
    }
    for(i=0;i<NumberOfMaterials;i++){
      message(1,"Material %i found in %i cells.\n",i,ListOfMaterials[i].CellCount);
    }
    if(detailed_sputtering==1){ /* For detailed calucation of sputtering,
				   we need to know which neighbors of each cell are vacuum */
      DetermineVacuumNeighbors();
    }
    if(UseDensityMult==1){  /* If applicable, read the density multiplicators */
      result=ReadFloatFileIntoArray(TargetDensityMultFileName,TargetDensityMult, cell_count, TargetCompositionFileType);
      if(result!=0){
	message_error(result,"cannot read density multiplicator file %s (error %i).\n",TargetDensityMultFileName,result);
	return result;
      }
      message(0,"Density multiplicators read from %s.n",TargetDensityMultFileName);
    }
		
  } else { /* Only calc memory usage */
    lui_temp+=7*sizeof(int)*cell_count;
    lui_temp+=sizeof(float)*cell_count;
    if(detailed_sputtering==1){ /* For detailed calucation of sputtering, we need these */
      lui_temp+=sizeof(int)*cell_count;
      lui_temp+=sizeof(char)*cell_count;
    }
    if(store_energy_deposit==1){ /* For detailed storing of deposited energy */
      lui_temp+=2*sizeof(double)*cell_count;
    }
    for(i=0;i<NumberOfMaterials;i++){
      for(j=0;j<ListOfMaterials[i].ElementCount;j++){ /* Go through elements, make arrays for each one */
	lui_temp+=4*cell_count*sizeof(int);
	if(detailed_sputtering==1){ /* For detailed calucation of sputtering, we need these */
	  lui_temp+=cell_count*sizeof(int);
	}
      }
    }
    mem_usage+=lui_temp;
    if(mem_usage_details==1){printf("MEMORY Target historgrams:          %li bytes\n",lui_temp);}
  }


  return 0;
}

int TargetIndex(int x, int y, int z){
  /* Calcualtes index for accessing one-dimensional target
     arrays knowing the integer cell coordinates */
  return x*layer_count_yz+y*cell_count_z+z;
}


/* TO DO: Idea: use reciprocal values for 1/cell_count and so on, to speed up calculations */
int GetTargetXYZ(int index, int* x, int* y, int* z){
  /* Calculates the x,y,z coordinates (cell numbers) for a given index ... so the inverse function of TargetIndex() */ 
  *z=(index%layer_count_yz)%cell_count_z;
  *y=((index%layer_count_yz)-(*z))/cell_count_z;
  /* *y=((index%layer_count_yz)-(*z))/cell_count_z;  */
  *x=((index-(*z)-cell_count_z*(*y))/layer_count_yz); 
  return 0;
}

int GetCellIndex(float x, float y, float z){
  /* Calculates index for accessing one-dimensional target
     arrays from a float-value point in space */
  int ix,iy,iz;
  ix=x/cell_size_x;
  iy=y/cell_size_y;
  iz=z/cell_size_z;

  // SAFETY TEST:
  /*  if(ix>=cell_count_x){printf("Index Warning! ix: %i \n",ix);}
      if(iy>=cell_count_y){printf("Index Warning! iy: %i \n",iy);}
      if(iz>=cell_count_z){printf("Index Warning! iz: %i \n",iz);} */
  return TargetIndex(ix,iy,iz);
}


int PositionInCell(float x, float y, float z, float* rx, float* ry, float* rz){
  /* Calculates the relative position (rx, ry, rz) within a cell from an absolute position (x,y,z) */
  /* returns 0 on success */
  *rx= x-( ((int)(x/cell_size_x))*cell_size_x ); 
  *ry= y-( ((int)(y/cell_size_y))*cell_size_y ); 
  *rz= x-( ((int)(z/cell_size_z))*cell_size_z );
  return 0;
}

int DetermineVacuumNeighbors(){
  /* Creates an array where for each cell a char is stored that
     describes, which neighbors are vacuum.
     This is done bitwise:
     bit 0:  next cell in +x direction is vaccum
     bit 1:  -x
     2: +y
     3: -y
     4: +z
     5: -z
     Cells that are at the front of the target (x=0), have vacuum in -x direction!
  */
  int x,y,z;
  int cellindex;
  int isvac;
  for(x=0;x<cell_count_x;x++){ /* Loop through all cells */
    for(y=0;y<cell_count_y;y++){  /* if a cell is vacuum, tell all its neighbors, that they */
      for(z=0;z<cell_count_z;z++){  /* have vacuum in their neighborhood */
	cellindex=TargetIndex(x,y,z);
	isvac=ListOfMaterials[TargetComposition[cellindex]].Is_Vacuum;
	if(isvac==1){ /* cell is vacuum */
	  if((x<(cell_count_x-1))&&(cell_count_x>1)){ /* indicate vac neighbor for next cell in +x direc */
	    TargetVacuumNeighbors[TargetIndex(x+1,y,z)]+=2; /* Bit 1 indicates vacuum in -x direction */
	  } else {
	    if(boundary_x==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(0,y,z)]+=2;
	    }
	  }
	  if(x>0){ /* indicate vac neighbor for next cell in -x direc */
	    TargetVacuumNeighbors[TargetIndex(x-1,y,z)]+=1; /* Bit 0 indicates vacuum in +x direction */
	  } else {
	    if(boundary_x==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(cell_count_x-1,y,z)]+=1;
	    }
	  }

	  if((y<(cell_count_y-1))&&(cell_count_y>1)){ /* indicate vac neighbor for next cell in +y direc */
	    TargetVacuumNeighbors[TargetIndex(x,y+1,z)]+=8; /* Bit 3 indicates vacuum in -y direction */
	  } else {
	    if(boundary_y==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(x,0,z)]+=8;
	    }
	  }
	  if(y>0){ /* indicate vac neighbor for next cell in -y direc */
	    TargetVacuumNeighbors[TargetIndex(x,y-1,z)]+=4; /* Bit 2 indicates vacuum in +y direction */
	  } else {
	    if(boundary_y==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(x,cell_count_y-1,z)]+=4;
	    }
	  }

	  if((z<(cell_count_z-1))&&(cell_count_z>1)){ /* indicate vac neighbor for next cell in +z direc */
	    TargetVacuumNeighbors[TargetIndex(x,y,z+1)]+=32; /* Bit 5 indicates vacuum in -Z direction */
	  } else {
	    if(boundary_z==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(x,y,0)]+=32;
	    }
	  }
	  if(z>0){ /* indicate vac neighbor for next cell in -z direc */
	    TargetVacuumNeighbors[TargetIndex(x,y,z-1)]+=16; /* Bit 4 indicates vacuum in +z direction */
	  } else {
	    if(boundary_z==1){ /* PBC */
	      TargetVacuumNeighbors[TargetIndex(x,y,cell_count_z-1)]+=16;
	    }
	  }
	}
	if(x==0){ /* Cells that are at the front of the target (x=0), have vacuum is -x direction! */
	  TargetVacuumNeighbors[cellindex]|=2;
	}
      } /* EOF z */
    } /* EOF y */
  } /* EOF x */
  return 0;
}

int CheckSurfaceAtom(int cellindex, float x, float y, float z, float spacing){
  /* Returns 0 if the given position is in bulk, 1 or higher if the given position corresponds
     to a surface atom: the returned number indicates in which direction the neighbouring vacuum is.
     spacing is the maximum distance from atom to vacuum to be
     considered a surface atom.
     cellindex needs to be provided to speed up calculations... its known anyway. */

  float rx,ry,rz;   /* relative position inside cell */

  PositionInCell(x,y,z,&rx,&ry,&rz);   /* Determine relative position inside cell: */

  /* Check for each direction if this is close to vacuum: */

  if(ry<=spacing){ /* close in -y direction */
    if((TargetVacuumNeighbors[cellindex]&(char)8)==(char)8){return 8;}
  }
  if(ry>=(cell_size_y-spacing)){ /* close in +y direction */
    if((TargetVacuumNeighbors[cellindex]&(char)4)==(char)4){return 4;}
  }
  if(rx>=(cell_size_x-spacing)){ /* close in +x direction */
    if((TargetVacuumNeighbors[cellindex]&(char)1)==(char)1){return 1;}
  }
  if(rx<=spacing){ /* close in -x direction */
    if((TargetVacuumNeighbors[cellindex]&(char)2)==(char)2){return 2;}
  }
  if(rz>=(cell_size_z-spacing)){ /* close in +z direction */
    if((TargetVacuumNeighbors[cellindex]&(char)16)==(char)16){return 16;}
  }
  if(rz<=spacing){ /* close in -z direction */
    if((TargetVacuumNeighbors[cellindex]&(char)32)==(char)32){return 32;}
  }
  /* Not close to vacuum in any direction, return 0 */
  return 0;
}



/* ---- now come some helping functions concerning file reading ----- */


int PrepareScatteringMatrix(struct scattering_matrix * ScatMatrix, float ProjM, int ProjZ, float TarM, int TarZ){
  /* Creates the matrix with scattering results and stores it to the structure pointed to by ScatMatrix */

  /* Adapted from the corteo code */

  ScatMatrix->screening_length     = d2f(0.1f * SCREENCONST/ ( pow(ProjZ,0.23) + pow(TarZ,0.23) )); /* This is in nm! */
  ScatMatrix->inv_screening_length = d2f(1.0f/ScatMatrix->screening_length);
  ScatMatrix->mass_ratio           = ProjM/TarM;
  ScatMatrix->sqrt_mass_ratio      = sqrtdf(ProjM/TarM);
  ScatMatrix->kfactor_m            = 4.0f*ScatMatrix->mass_ratio / ( (1.0f+ScatMatrix->mass_ratio) * (1.0f+ScatMatrix->mass_ratio)  ); /* k factor without the angle part */
  ScatMatrix->red_E_conv           = 10.0f * ScatMatrix->screening_length / (  (1.0f+ScatMatrix->mass_ratio) * ProjZ * TarZ * E2); /* We need screening length in Angstroms since E2 is given in eVA */

  /* The maximum reduced impact parameter depends on density and flight length. This
     means that we cannot calculate it as a function of projectile and target element only.
     Thus we need to calculate it somewhere else (as a material property of target). This
     is different from corteo, where the scattering matrices are stored for each target
     layer with known density etc. */

  /* ScatMatrix->max_red_im_par   = 1.0f / ( ScatMatrix->screening_length * sqrtdf( PI * a )  ) */
  /* ion->smax[ilayer][ielem] = 1.0f/(ion->a[ilayer][ielem]*sqrtdf(PI*layerDensity[ilayer])*sqrtMeanFreePath[ilayer]);  */
  /* float max_red_im_par;       Maximum reduced impact parameter */

  ScatMatrix->CosScat          = (float*)malloc(DIME*DIMS*sizeof(float));
  if( (ScatMatrix->CosScat)==NULL ){ return -1; } /* Cannot allocate memory */
  ScatMatrix->SinScat          = (float*)malloc(DIME*DIMS*sizeof(float));
  if( (ScatMatrix->SinScat)==NULL ){ return -1; } /* Cannot allocate memory */

  fillCosSinTable(ScatMatrix->CosScat, ScatMatrix->SinScat, ScatMatrix->mass_ratio);

  return 0;

}

int TargetStructure_DataBlockReader(char* BlockName){
  /* Needs to be called from the ini file reader
     while the traget structure input file is read */
  /* We won't need this at the moment... */
  return 0;
}

int TargetStructure_DataReader(char* ParName, char* ParValue){
  /* Needs to be called from the ini file reader
     while the target structure input file is read */

  if(strcmp(ParName,"cell_count_x")==0){ /* ... */
    sscanf(ParValue,"%i",&cell_count_x);
  }
  if(strcmp(ParName,"cell_count_y")==0){ /* ... */
    sscanf(ParValue,"%i",&cell_count_y);
  }
  if(strcmp(ParName,"cell_count_z")==0){ /* ... */
    sscanf(ParValue,"%i",&cell_count_z);
  }

  if(strcmp(ParName,"cell_size_x")==0){ /* ... */
    sscanf(ParValue,"%f",&cell_size_x);
  }
  if(strcmp(ParName,"cell_size_y")==0){ /* ... */
    sscanf(ParValue,"%f",&cell_size_y);
  }
  if(strcmp(ParName,"cell_size_z")==0){ /* ... */
    sscanf(ParValue,"%f",&cell_size_z);
  }

  if(strcmp(ParName,"periodic_boundary_x")==0){ /* ... */
    sscanf(ParValue,"%i",&boundary_x);
  }
  if(strcmp(ParName,"periodic_boundary_y")==0){ /* ... */
    sscanf(ParValue,"%i",&boundary_y);
  }
  if(strcmp(ParName,"periodic_boundary_z")==0){ /* ... */
    sscanf(ParValue,"%i",&boundary_z);
  }

#ifdef INCLUDE_SPECIAL_GEOMETRY
  if(strcmp(ParName,"special_geometry")==0){ /* ... */
    sscanf(ParValue,"%i",&special_geometry);
    message(0,"Special geometry in use: %i\n",special_geometry);
  }
#endif
	
  if(strcmp(ParName,"CompositionFileType")==0){ /* ... */
    sscanf(ParValue,"%i",&TargetCompositionFileType);
  }

  if(strcmp(ParName,"CompositionFileName")==0){ /* ... */
    /* TargetCompositionFileName=malloc(sizeof(char)*(strlen(ParValue)+1));
       if(TargetCompositionFileName==NULL){return -1;} */
    if(single_input_file!=1){
      strcpy(TargetCompositionFileName,ParValue);
    }
  }

  if(strcmp(ParName,"UseDensityMultiplicator")==0){ /* ... */
    sscanf(ParValue,"%i",&UseDensityMult);
  }
  if(strcmp(ParName,"DensityMultiplicatorFileName")==0){ /* ... */
    TargetDensityMultFileName=malloc(sizeof(char)*(strlen(ParValue)+1));
    if(TargetDensityMultFileName==NULL){return -625;}
    strcpy(TargetDensityMultFileName,ParValue);
  }
  return 0;
}

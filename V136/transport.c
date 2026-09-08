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


/***************************************************************************
	This module contains the functions to simulate the ion transport, etc.   
	There are currently two transport functions:
	
	FullProjectileTransport() simulates full transport accurately.
	FastProjectileTransport() works similar to the corteo transport function.
	It is faster, but makes a few more approximations. It is accurate enough
	for ion and coarse damage distribution but not for sputtering!
	
	Modifications by Jean-Paul Crocombette, CEA Saclay for modified
	Kinchin-Pease calculation of damage and free flight path approximation,
	to allow comparison with SRIM.
	They  are indicated by "CROC" comments
	
	Error numbers supposed to be used in transport module: 1000-1999      
	
****************************************************************************/

//#define MONITOR_ION 1000000000000

#include "transport.h"

/* DJGPP will not link the program correctly under windows, because it doesn't know sqrtf.
   So we'll define it by using standard sqrt for doubles*/
#ifdef DJGPP
float sqrtf(float x){return (float)sqrt(x);}
#endif

int IrradiateTarget(){
  /* Let ions impinge at the entry side of the target */
  /* The entry positions can be random or more defined */

  /*CROC x=0 initialization to allow from inside bulk irradiation */
  float x,y,z,energy;
  int   i;    /* count ions */
  x=0.0;   y=0.0;   z=0.0;

  sim_start_time=time(NULL);

  /* init counters: */
	disp_c=0;
	miss_c=0;
	coll_c=0;
	sputter_c=0;
	ion_c=0;
	repl_c=0;
	leaving_ions_c=0;
	leaving_recoils_c=0;
	vac_c=0;
	int_c=0;
	escape_solid_c=0;
	replacer_escaped=0;
	recoil_counter=0;

	/*
	subreplace_c=0;
	superthres_leave_c=0;
	subthres_created_c=0;
	subthres_reach_egde_c=0;
	superthres_created_counter_c=0;
	sub_become_super_c=0;
	*/

  if(store_ion_paths==1){ /* if requested, store the exact paths of the ions */
    ion_paths_fp=OpenFileContinuous(OutputFileBaseName,".ionpaths");
    if(ion_paths_fp==NULL){
      message_error(-1,"Cannot open file for storing ion paths!\n");
      return -1;
    }
  }
  if(store_recoil_cascades==1){ /* if requested, store the exact paths of the ions */
    recoil_cascades_fp=OpenFileContinuous(OutputFileBaseName,".cascades");
    if(recoil_cascades_fp==NULL){
      message_error(-1,"Cannot open file for storing recoil cascades!\n");
      return -1;
    }
  }
  /*CROC if store_range3d=2 storage of interstitials and vacancy positions*/
  if(store_range3d>=1){ /* if requested, store the final positions of ions, similar to TRIMs range_3d.txt (only ions stopped inside target!) */
    store_range3d_fp=OpenFileContinuous(OutputFileBaseName,".ions.range3d");
    if(store_range3d_fp==NULL){
      message_error(-1002,"Cannot open range3d-file!\n");
      return -1002;
    }
  }
  if(store_PKA>=1){ /* if requested, store each PKA with type, positions, direction, energy, ...!) */
    store_PKA_fp=OpenFileContinuous(OutputFileBaseName,".PKAs");
    if(store_PKA_fp==NULL){
      message_error(-1005,"Cannot open PKA file!\n");
      return -1005;
    } else {
      /* header */
      fprintf(store_PKA_fp,"#x\ty\tz\tvx\tvy\tvz\tmat_index\telem_index\tZ\tm\tenergy\n");
      fprintf(store_PKA_fp,"#nm\tnm\tnm\t1\t1\t1\ti\ti\t1\tAMU\teV\n");
    }
  }
  if(store_range3d==2){ /* CROC: if requested, store the final positions of interstitials and vacancies*/
    store_range3dV_fp=OpenFileContinuous(OutputFileBaseName,".ions.range3dV");
    if(store_range3dV_fp==NULL){
      message_error(-1003,"Cannot open range3dV-file!\n");
      return -1003;
    }
    store_range3dI_fp=OpenFileContinuous(OutputFileBaseName,".ions.range3dI");
    if(store_range3dI_fp==NULL){
      message_error(-1004,"Cannot open range3dI-file!\n");
      return -1004;
    }
  }


  single_ion_sputter_counter=0;

  /*  fEnergy=OpenFileContinuous(OutputFileBaseName,".SurfEnergy"); */
  /* CROC call to building of KP tables if simulation_type==5*/
  if(simulation_type==5){ 
    prepare_KP_tables2 ();
  }

  if (simulation_type==7){
    x=target_size_x/2.0;
    y=target_size_y/2.0;
    z=target_size_z/2.0;
    minE=10;
    int N = 100;
    
    float r = powf(ionInitialEnergy / minE, 1.0f / (N - 1));
    
    energy=minE;
    
    while(energy<=ionInitialEnergy){
      totlP=0;
      totlE=0;
      for(i=0;i<max_no_ions;i++){ /* Create the ions and let them impinge on the target */
	FastProjectileTransport(ionZ,ionM,energy,x,y,z,ion_vx,ion_vy,ion_vz,1,0,0,0,0);
      }
      printf("%f %f   \n",energy,totlE/totlP);
      fprintf(ball_stopping_fp,"%g\t 0\t 0\t %g\n",energy,totlE/totlP);
      energy*=r;
      }
    }else{		  
      
      for(i=0;i<max_no_ions;i++){ /* Create the ions and let them impinge on the target */
     

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"\n");
#endif
    /* Write status file */
    if(create_status_file==1){ /* Status file, which can be read by other programs */
      if(i%status_update_interval==0){
	write_status_file("sim", i);
      }
    }

    /* Display progress */
    if(i%display_interval==0){message(-1,"Ions completed: %i\n",i);fflush(stdout);}
    ion_c++;

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"\n");
#endif

    /* Store ion paths */
    if(i==store_path_limit){ /* Stop storing paths after so many ions */
      /* Close files if necessary */
      if(store_ion_paths==1){
	fflush(ion_paths_fp);
	fclose(ion_paths_fp);
      }
      store_ion_paths=0; 	/* Stop further logging */
      message(1,"Ion paths are no longer stored after %i ions.\n",i);fflush(stdout);
    }
		if(i==store_path_limit_recoils){ /* Stop storing recoil paths after so many ions */
      if(store_recoil_cascades==1){
	fflush(recoil_cascades_fp);
	fclose(recoil_cascades_fp);
      }
      store_recoil_cascades=0; /* Stop further logging */
      message(1,"Recoil cascades are no longer stored after %i ions.\n",i);fflush(stdout);
    }

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"\n");
#endif

    if(store_ion_paths==1){ /* Make sure data get into file */
      fflush(ion_paths_fp);
    }
    if(store_recoil_cascades==1){
      fflush(recoil_cascades_fp);
    }
    if(store_range3d==1){fflush(store_range3d_fp);}
    if(store_range3d==2){fflush(store_range3dV_fp);}
    if(store_range3d==2){fflush(store_range3dI_fp);}
    if(store_PKA>=1)    {fflush(store_PKA_fp);}

    if(single_ion_sputter_yields==1){ /* store single ion sputter yields */
      if(single_ion_sputter_counter<=MAX_SPUTTERED){
	sputter_yield_histogram[single_ion_sputter_counter]++;
	single_ion_sputter_counter=0;
				
      }
    }
		
    if(((i%storage_interval)==0)&&(i>0)){ /* Intermediate storing */

      /*      message(-1,"Storing results after ion no. %i, %i displacements, %i replacements, %i missed collision...",i,disp_c,repl_c,miss_c);fflush(stdout);*/
      message(-1,"Storing results after ion no. %i ...\n",i);fflush(stdout);
      calculate_normalization_factor(i);
			
      store_results(OutputFileBaseName,i);
      message(-1,"done.\n");

    }

    if(store_ion_paths==1){ /* store empty line in ion paths file to separate ions*/
      fprintf(ion_paths_fp,"\n");
    }

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"\n");
#endif

    /**************************************** Define ion entry position and let ion start ******/
    /* Define ion entry positions: */
    x=0.0;
    switch(ion_distribution){
    case 0: /* Random */
      y=randomx()*target_size_y;
      z=randomx()*target_size_z;
      /*printf("Entry pos: y: %g, z: %g\n",y,z);*/
      break;
    case 1:  /* centered */
      y=target_size_y/2.0;
      z=target_size_z/2.0;
      break;
    case 2:  /* defined position */
      x=enter_x;
      y=enter_y;
      z=enter_z;
      break;
    case 3:  /* random square around predefined position */
      y=enter_y+(0.5-randomx())*beam_spread;
      z=enter_z+(0.5-randomx())*beam_spread;
      break;
    case 4:  /* CROC : BULK centered initial position */
      x=target_size_x/2.0;
      y=target_size_y/2.0;
      z=target_size_z/2.0;
      break;
    } /* Plan: gaussian beam */


#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"\n");
#endif

    /*CROC indicates in the defect file the beginning of a new cascade*/
    if(store_range3d==2){
      fprintf(store_range3dV_fp,"C %i\n",i);
      fprintf(store_range3dI_fp,"C %i\n",i);
    } 

    /* Call function to simulate ion */
    switch(transport_type){
    case 0: /* accurate */
      FullProjectileTransport(ionZ,ionM,ionInitialEnergy,x,y,z,ion_vx,ion_vy,ion_vz,1,0,0,0,0,0);
      break;
    default: /* fast, like corteo */
      FastProjectileTransport(ionZ,ionM,ionInitialEnergy,x,y,z,ion_vx,ion_vy,ion_vz,1,0,0,0,0);
      break;
    }

  }
  }

#ifdef DEBUG_MODE
  printf("DEBUG %s, l %i\n",__FILE__,__LINE__);fflush(stdout);
#endif


  /* Print the counters */
  /* Note that this sputtering type is extremely simplistic and should not be used for the sputtering yield. */
  message(3,"Counters:\n\tIons:\t\t%i\n\tDisp:\t\t%i\n\tMissed:\t\t%i\n\tColl:\t\t%i\n\tRepl:\t\t%i\n",i,disp_c,miss_c,coll_c,repl_c);

  message(3,"\tVacs:\t\t%i\n\tInts:\t\t%i\n\tIons left:\t%i\n\tRecoils left:\t%i\n\tEscaped solid:\t%i\n\tRepl.escaped:\t%i\n",vac_c,int_c,leaving_ions_c,leaving_recoils_c,escape_solid_c,replacer_escaped);
  /*  message(3,"\tSubCrea:\t%i\n\tSubReach:\t%i\n\tSuperCrea:\t%i\n\tSubSuper:\t%i\n\tSuperLeave:\t%i\n\tSubReplace:\t%i\n",subthres_created_c,subthres_reach_egde_c,superthres_created_counter_c,sub_become_super_c,superthres_leave_c,subreplace_c); */
  message(3,"\tAtoms lost via vacuum (simplistic sputtering):\t%i\n",sputter_c);

  /* Close files that have been opened for continuous output: */
  if(store_ion_paths==1){fclose(ion_paths_fp);}
  if(store_recoil_cascades==1){fclose(recoil_cascades_fp);}
  if(store_range3d>=1){fclose(store_range3d_fp);}
  if(store_range3d==2){fclose(store_range3dV_fp);}
  if(store_range3d==2){fclose(store_range3dI_fp);}
  if(store_PKA>=1)    {fclose(store_PKA_fp);}
  return 0;
}

int FastProjectileTransport(int ProjZ, float ProjM, double ProjE, float Proj_x, float Proj_y, float Proj_z, float Proj_vx, float Proj_vy, float Proj_vz, int is_ion, int OrgMaterial, int OrgElement, int OrgCell, int RD){
  /* Calculates the path of a projectile with given properties (proton number, mass,
     energy, coordiantes, velocity unit vector (directional cosines) through the target material. 
     Can be called recursively to follow recoils. The "is_ion" parameter must be 1 when an initial
     ion is handled or 0 when recoils are simulated, so the function will know where to store
     results and which scattering matrix is needed (shouldn't be necessary for the physics though).
     If the projectile is a recoil from within the target, the material and element and cell that
     it originated from must be known in order to store correctly the distribution of implanted
     recoils for each material. */
  /* RD: recursion depth. 0: ion, 1: PKA, 2: secondary... */

  /* What this function does:
     It performs a loop with:
     - Find out in what material we currently are
     - Obtain a free flight path
     - Calculate electronic stopping along the path
     - calculate new position (obey boundary conditions)...
     - simulate a collision
     - if recoil gets enough energy:
     -   store a vacancy
     -   simulate the recoil recursively
     - if projectile has low energy:
     -   stop it and store it as ion or interstitial
     - else: advance to next collision site
  */
	
  int cell_i=0;               /* Current cell index */

  int current_material_index=0;
  struct material* current_material; /* Pointer to material we are using */
  //  struct material* old_material;     /* Index of material into which the projectile moves after collision */
  int i=0;
  int ion_left_target=0;
  int leaving_direction=0;    /* Direction in which a sputtered atom leaves the target */
  float flightpath=0;         /* length of free flight until next collision [nm] */
  float x,y,z;                /* Current point in space where projectile is [nm]*/
  float vx,vy,vz;             /* current velocity unit vector */
  float old_vx,old_vy,old_vz; /* previous velocity unit vector */
  double energy;              /* projectile energy [eV] */
  float stopping;             /* electronic energy loss */
  float straggling;           /* electronic energy loss straggling */

  float random;               /* temporary random number */
  float conc_sum;             /* Summed up relative elemental concentrations */

  float target_mass=0;        /* for collision partner */
  int target_Z=1; 
  int target_index=0;         /* index of collision partner (within current material) */
  float impact_par=0.0;       /* impact parameter for collision */
  float red_impact_par;       /* reduced impact parameter for collision */
  float recoil_energy;        /* Energy transferred to recoil */
  float recoil_vx,recoil_vy,recoil_vz;   /* velocity unit vector of recoil */
  int recoil_is_at_surface;              /* is >0 if recoil is a surface atom. Notes also the direction of neighbouring vacuum cell */
  long int recoil_number;                /* Identifying unique number of this recoil */
  float E_frac;                          /* Fraction of energy directed into vacuum */
  float E_compare;                       /* the displacement energy that has to be used.
					    (In some cases it is the surface energy instead of the displacement energies) */
  int surface_sputtered;                 /* is 1, if atom was sputtered from surface */
  float e_disp,e_latt,e_surf=0,e_repl;   /* Displacement, lattice, surface and replacement energy of target atoms */
  int replaced=0;                        /* is 1 if a replacement occured */
  int proj_eq_target=0;                  /* if projectile and target are the same, this is 1 */
  struct scattering_matrix * ScatMatrix; /* Points to currently needed scattering matrix */
  int smIndex;                           /* Index inside scattering matrix */
  float sin2thetaby2;                    /* sin^2(theta/2) of collision */ 

  float temp1,temp2;                     /* some temporary values */
  int itemp1;
  float vel, svel;                       /* velocity and 1/sqrt(vel). Needed to correct
					    inaccuracies in rotation fct. */
  /*  int multi_col_counter; */          /* To allow multiple collisions as suggested in Eckstein, p. 92 */
  int isvac;                             /* is 1 if vacuum */
  float epsil;                           /* CROC SRIMepsillon*/
  float randomSRIM;                      /* CROC  */
  float   e_d, g_ed, E_v,  E_div,k_dr;   /* CROC : SRIM variables for free flight path*/
  float xsi ,pmax=.1,bmax ;              /* 7.40-7.41 page 7-16 */
  /* pmax requires to be initialized in case the starting cell contains no material */

  /* init values */
  recoil_number=recoil_counter++; /* set unique ID of this recoil */
  x=Proj_x;
  y=Proj_y;
  z=Proj_z;
  vx=Proj_vx;
  vy=Proj_vy;
  vz=Proj_vz;
  energy=ProjE;
  //randomSRIM=d2f(randomx());           /* is done later if required */
	
  if(RD==1){ /* PKA! */
    if(store_PKA>=1){ /* if requested, store each PKA with type, positions, direction, energy, ...!) */
      fprintf(store_PKA_fp,"%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%f\t%g\n",x,y,z,vx,vy,vz,OrgMaterial,OrgElement,ProjZ,ProjM,ProjE);
    }
  }
	
  recoil_vx=0;
  recoil_vy=0;
  recoil_vz=0;

  /* if(is_ion==1){printf("DEBUG transport.c: line %i, ION Energy %g\n",__LINE__,energy);}
     else { 
     printf("DEBUG transport.c: line %i,  Recoil Energy %g\n",__LINE__,energy);
     printf("DEBUG info: xyz: %f,%f,%f energy: %g \n",x,y,z,energy);
     printf("DEBUG info: vx vy vz: %f,%f,%f energy: %g \n",vx,vy,vz,energy);
     } */
	
	
  /*CROC: if from inside bulk cascade the initial site is a vacancy! */
  if (ion_distribution==4){
    if(is_ion==1){
      if(store_range3d==2){fprintf(store_range3dV_fp,"%g\t%g\t%g\n",x,y,z);}             
    }
  }
	
  /* Determine initial cell and material */
  cell_i=GetCellIndex(x,y,z);
	
#ifdef INCLUDE_SPECIAL_GEOMETRY
  if(special_geometry==0){ /* use standard grid to determine material */
    current_material_index=TargetComposition[cell_i];
    current_material=&(ListOfMaterials[current_material_index]);
  } else { /* for special geometries, this might be different */
    current_material_index=GetMaterialFromPosition(cell_i,x,y,z);
  }
#else
  current_material_index=TargetComposition[cell_i];
#endif
  current_material=&(ListOfMaterials[current_material_index]);

  while(energy>0){ /* The projectile still has energy and can move on */
    replaced=0;    /* Assume at first, that no replacement occurs */
    if(is_ion==1){
      if(store_ion_paths==1){ /* if requested, store position of ion */
	fprintf(ion_paths_fp,"%g\t%g\t%g\t%g\n",x,y,z,energy);
      }
    } else { /* Recoil */
      if(store_recoil_cascades==1){ /* if requested, store position of recoil */
	fprintf(recoil_cascades_fp,"%g\t%g\t%g\t%g\t%li\n",x,y,z,energy,recoil_number);
      }
    }
		
    isvac=current_material->Is_Vacuum;
    if(isvac==0){ /* if real material, not vacuum */
      /* Calculate free flight path (note that this is strongly correlated with the calculation of the impact parameter (see below) ): */
      switch(flight_length_type){
      case 0: /* Poisson distributed flight length and impact pars */
	/* Flight length between collisions: Poisson-distributed value with interatomic spacing as mean free path */
	/* flightpath = -current_material->AtomicDistance * log(randomx()); */
	temp1= sqrtloglist[iranloglist]*current_material->SqrtAtomicDistance; /* Use of random list faster (as in corteo)*/
	flightpath=temp1*temp1;
	break;
      case 1: /* atomic spacing */
	flightpath=current_material->AtomicDistance;
	break;
      case 2: /* constant flightpath in nm*/
	flightpath=flight_length_constant;
	break;
      case 3: /*CROC SRIM like free flight path*/
	epsil=energy*current_material->MeanF;
	/*	if(is_ion==1) {	printf(" %f  %f \n ", energy,epsil);}*/
	/* 	printf(" %f  %f \n ", energy,epsil);*/
	xsi=pow( epsil*current_material->MeanMinRedTransfer,0.5) ;
	bmax=1/(xsi+pow(xsi,0.5)+0.125*pow(xsi,0.1));
	pmax=bmax*current_material->MeanA;
	/*printf("pmax=%f\n",pmax);*/
	flightpath=1.0f / (PI*current_material->DensityNM*pmax*pmax) ;
	break;
      }
      if (simulation_type==7  ){
	totlP+=flightpath;
      }
      /* TO DO: Perhaps impose a limit on unrealstically long flight paths? --> this has only little influence */
      /* printf("DEBUG transport.c: line %i  flight length %f, ad %f\n",__LINE__,flightpath,current_material->AtomicDistance); */
			
      /* Calculate electronic stopping(+-straggling): */
      switch (stopping_model) {
        case 0:
	  stopping   = (flightpath * 10.0) * ElectronicStopping(ProjZ,ProjM,energy,current_material_index); /* factor 10, because flightpath is in nm */
	  break;
      case 1:
	stopping   = (flightpath * 1.0) * ElectronicStoppingdpass(ProjZ,ProjM,energy,current_material_index); 
	break;
      case 2:
	stopping   = (flightpath ) * ElectronicStoppingfile(ProjZ,energy,current_material_index); 
	break;
      default:
	printf("wong electronic stopping_model\n");
	exit(1);
	  }

      /* if (stopping_model==0){ */
      /* 	stopping   = (flightpath * 10.0) * ElectronicStopping(ProjZ,ProjM,energy,current_material_index); /\* factor 10, because flightpath is in nm *\/ */
      /* } else { */
      /* 	stopping   = (flightpath * 1.0) * ElectronicStoppingdpass(ProjZ,ProjM,energy,current_material_index); /\* factor 10, because flightpath is in nm *\/ */
      /* } */
      straggling = 3.16227766016838 * sqrt(flightpath) *  ElectronicStraggling(ProjZ,ProjM,energy,current_material_index) * inverse_erf_list[erflist_pointer++] ;      /* factor is sqrt(10) */
      /* Due to gaussian distribution, the straggling can in some cases get so large that the projectile gains energy or suddenly looses a huge amount of energy. Is this realistic? This might actually happen. However, in the simulation, ions may have higher energy than the initial energy.
	 We will therefore limit the |straggling| to |stopping|.  Furthermore, with hydrogen, the straggling is often so big, that ions gain huge amount of energy, and the phononic system would actually gain energy. */
      if(fabs(straggling)>stopping){
	if(straggling<0){
	  straggling=-stopping;
	} else {
	  straggling=stopping;
	}
      }
      /* The stopping tables have no values below 16 eV. Therefore, we do simple linear downscaling of electronic stopping below 16 eV. */
      if(energy<16){ /* Experimental scaling down */
	stopping*=(energy*0.0625);
	straggling*=(energy*0.0625);
      }
      if(erflist_pointer>=MAXERFLIST){erflist_pointer=0;} /* Check and adjust boundary for ERF index*/
      /* Now subtract electronic stopping from projectile energy: */
      energy-= (stopping+straggling);
      if(store_energy_deposit==1){
	TargetEnergyElectrons[cell_i]+=(double)(stopping+straggling);
      }
    } else { /* We are in vacuum: no stopping happens here, fly along some way */

      /* flightpath=0.3+randomx()*0.1; */
      flightpath=0.3;
      /* Calcualting the end point where the projectile exits the current cell would be nicer,
	 however this would also take some cpu time... so we just let is fly in small steps.
	 With 0.3 nm we shouldn't miss target material */
      /* Note further: A constant path is a little dangerous here; it can lead to artifacts
	 when calculating sputtering. Anyway, for the fast transport process this is not so important.
	 Nevertheless you might consider using the upper line with the random value */
    }

    /* Now let projectile fly and calculate new position: */
    x+=flightpath*vx;
    y+=flightpath*vy;
    z+=flightpath*vz;

    /* Check if position has exceeded boundaries, correct if necessary */
    ion_left_target=0;
    ion_left_target+=CheckAndCorrectBoundary(&x,target_size_x,&target_max_x,boundary_x);
    ion_left_target+=CheckAndCorrectBoundary(&y,target_size_y,&target_max_y,boundary_y);
    ion_left_target+=CheckAndCorrectBoundary(&z,target_size_z,&target_max_z,boundary_z);
		
    if(ion_left_target!=0){ /* projectile left sample */
      /*  if it is the ion that has left the target, store transmitted ion: */
      if(is_ion==1){
	if(store_transmitted_ions==1){
	  transmit_list[transmission_pointer].x=x;
	  transmit_list[transmission_pointer].y=y;
	  transmit_list[transmission_pointer].z=z;
	  transmit_list[transmission_pointer].vx=vx;
	  transmit_list[transmission_pointer].vy=vy;
	  transmit_list[transmission_pointer].vz=vz;
	  transmit_list[transmission_pointer].energy=energy;
	  transmission_pointer+=1;
	  if(detailed_sputtering==1){ /* if detailed sputtering is on, we also take a look on the ions leaving the target */
	    if(x>target_size_x){leaving_direction=0;}
	    if(x<=0){leaving_direction=1;}
	    if(y>target_size_y){leaving_direction=2;}
	    if(y<=0){leaving_direction=3;}
	    if(z>target_size_z){leaving_direction=4;}
	    if(z<=0){leaving_direction=5;}
	    leaving_ions[leaving_direction]+=1;
	  }
	}
      } else { /* it's a recoil leaving the sample */
	if(store_exiting_recoils==1){ 
	  itemp1 = (ListOfMaterials[OrgMaterial]).leaving_recoils_pointer[OrgElement];
	  if(itemp1<store_exiting_limit){ /* not too many recoils yet */
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).x = x;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).y = y;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).z = z;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).vx = vx;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).vy = vy;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).vz = vz;
	    (((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]).energy = energy;
	    (ListOfMaterials[OrgMaterial]).leaving_recoils_pointer[OrgElement]+=1;
	  }
	}
	/* Simple sputter counting: a sputtered atom is one that went through vacuum before it left the sample: */
	isvac=current_material->Is_Vacuum;
	if(isvac==1){
	  sputter_c+=1;
	}
	/* More complex sputter counting (only takes place, if detailed sputtering is switched on:) */
	/* for each target material and element, store the number of atom ejected out of the simulation volume into each of the possible 6 directions. */
	if(detailed_sputtering==1){	
	  /* Determine direction of leaving out of target: */
	  if(x>target_size_x){leaving_direction=0;}
	  if(x<=0){leaving_direction=1;}
	  if(y>target_size_y){leaving_direction=2;}
	  if(y<=0){leaving_direction=3;}
	  if(z>target_size_z){leaving_direction=4;}
	  if(z<=0){leaving_direction=5;}
	  ListOfMaterials[OrgMaterial].SputterCounter[(OrgElement*6)+leaving_direction]++;
	  (ListOfMaterials[OrgMaterial].TargetSputteredAtoms[OrgElement])[OrgCell]++; /* With this, we can later see where ions are sputtered from */
	}

	if(store_recoil_cascades==1){ /* if cascades are to be stored, store empty line into file to separate record */
	  fprintf(recoil_cascades_fp,"\n");
	}
      }
      energy=-0.001; /* in order to exit loop */

      /* ########################################################################## */
    } else { /* Projectile is still in the sample, might collide */
      /* Determine current cell and material at new position */
      cell_i=GetCellIndex(x,y,z);
#ifdef INCLUDE_SPECIAL_GEOMETRY
      if(special_geometry==0){ /* use standard grid to determine material */
	current_material_index=TargetComposition[cell_i];
      } else { /* for special geometries, this might be different */
	current_material_index=GetMaterialFromPosition(cell_i,x,y,z);
      }
#else
      current_material_index=TargetComposition[cell_i];
#endif
      current_material=&(ListOfMaterials[current_material_index]);
      isvac=current_material->Is_Vacuum;
			
      /* If material has changed, we should correct stopping. However, since the
	 flightlengths are significantly smaller than the cell dimensions, this
	 can only cause very small errors and only in cases where the material
	 changes and also the stopping powers are very different! */
			
      /* Now, let's see if we make a collision: */
      if(isvac==0){ /* if real material, not vacuum */
	coll_c+=1;    /* increase collision counter */
	proj_eq_target=0; /* Assume at first, that projectile is not the same as target */

	/* Randomly select collision target according to concentration: */
	if(iranlist>=MAXRANLIST-2)iranlist = 0;
	random   = randomlist[(iranlist)++]; /* A list of precalculated random values is used like in corteo */
	conc_sum = 0;

	for(i=0;i<current_material->ElementCount;i++){    /* go through elements in material... */
	  conc_sum+=(current_material->ElementsConc)[i];  /* ... and sum up concentrations */
	  if(random<conc_sum){ /* pick this element as target*/
	    target_index = i;
	    target_Z     = (current_material->ElementsZ)[target_index];
	    target_mass  = (current_material->ElementsM)[target_index];
	    if(target_Z==ProjZ){ /* check if projectile and target are the same: */
	      if(target_mass==ProjM){proj_eq_target=1;} /* Note, it can be dangerous to compare floats for
							   equality, however, the masses are never calculated but
							   always loaded, thus they should be exactly the same in memory
							   for the same target atoms */
	    }
	    break;
	  }
	}
	/* Now that we have selected the target atom, we know which of the scattering_matrices to use */
	if(is_ion==1){ /* It's the initial ion */
	  ScatMatrix=&(ion_scattering_matrix[target_Z]);
	} else { /* it's some other projectile */
	  ScatMatrix=&(scattering_matrices[ProjZ][target_Z]);
	}
				
	/* Select impact parameter (this depends also on the flight lengths used): */
	switch(flight_length_type){
	case 0: /* Poisson distributed flight length and impact pars */
	  impact_par = sqrtrandomlist[(iranlist)++] * current_material->MeanImpactPar * sqrtloglist1[iranloglist];
	  /* Note: the irandomlist pointer is already checked and set back to 0 at the selection of collision partner */
	  if(++iranloglist>=MAXLOGLIST) iranloglist=0;
	  break;
	case 1: /* atomic spacing */
	  impact_par = sqrtrandomlist[(iranlist)++] * current_material->MeanImpactPar; /* without the log list! */
	  break;
	case 2: /* constant */
	  impact_par = sqrtrandomlist[(iranlist)++] * current_material->SqrtRecFlDensity; /* without the log list! */
	  break;
	case 3: /*CROC SRIMlike */
	  randomSRIM=d2f(randomx());  /* Attention! Using random tables does not seem to be "sufficiently" random for large ion numbers in der KP case! --> randomx() is better here to get better statistics! */
	  impact_par = sqrt(randomSRIM) *pmax ;/* pmax decided before  */
	  /*	  impact_par = sqrtrandomlist[(iranlist)++] *pmax ; */ /* pmax decided before  */
	  /* Note: in case the projectile just came from vacuum into material: pmax is not properly defined.
	   * Solution: Use last used pmax. However, the very first pmax needs to be defined, if projectile starts in vacuum!
	   * */
	  break;
	}
				
	/* calculate reduced impact parameter s=p/a */
	red_impact_par=impact_par * ScatMatrix->inv_screening_length;
	/* Note that in contrast to corteo, the screening length AND impact par are in units of nm here,
	   so it doesn't matter for the reduced impact parameter */
				
	/* Calculate scattering angle */
	temp1=2.0*current_material->LayerDistance;
	if(impact_par >= temp1&& simulation_type!=7) {  /* impact parameter larger than interatomic distance: assume collision has missed */
	  /* Is this valid? Yes: tests have shown, that increasing the limit to 10 times the atomic distance
	     yields practically the same distribution of implanted ions and damage! */
	  if(flight_length_type <3) { /* missing should not occur for KP calculations! */
	    sin2thetaby2 = 0.0f;
	    recoil_energy = 0.0f;
	    miss_c+=1;
	  }
	} else { /* Collision takes place */
	  /* Get matrix index (using corteo's indexing functions): */
	  smIndex      = Eindex(energy * ScatMatrix->red_E_conv) * DIMS + Sindex(red_impact_par);
	  sin2thetaby2 = Matrix(smIndex); /* Obtain angle from general collision matrix using
					     reduced values of energy and impact par */
	  /* Energy transfer to recoil: */
	  recoil_energy = energy * ScatMatrix->kfactor_m * sin2thetaby2;
	  if (simulation_type==7  ){
	    totlE+=recoil_energy;
	    return;
	  }
		
	  energy       -= recoil_energy;
	  /*if one follows only the ion this energy loss should be attributed to phonons, i.e. ballistic losses*/
	  if(simulation_type==3){	  
	    if(store_energy_deposit==1){TargetEnergyPhonons[cell_i]+=(double)recoil_energy;} /* Energy goes to ballistic losses */
	  }

	  /* Store old flying direction */
	  old_vx = vx;
	  old_vy = vy;
	  old_vz = vz;
					
	  /* Calculate new ion direction by using cos and sin of scattering angle from precalculated matrices */
	  /* This is done by using the rotation routine from the corteo code. */
	  /* In the original version of corteo, there is a rare case of theta=180, which leads to sinTheta
	     being a "nan". This happens at extremely small impact parameters. --> has been corrected here */
	  rotate(&vx,&vy,&vz,&iazimAngle,(ScatMatrix->CosScat)[smIndex],(ScatMatrix->SinScat)[smIndex]);

	  /* after many collisions, the velocity starts to deviate from 1, because of the limited accuracy of floats. We might therefore correct it */
#ifdef RENORM_VELOCITIES
	  vel=vx*vx+vy*vy+vz*vz;
	  if(fabs(vel-1.9)>0.01){
	    svel=1.0/sqrtf(vel);
	    vx*=svel;vy*=svel;vz*=svel;
	  }
#endif
					
	  /* Errors should not occur often, but if they do, we can check if the velocity became "not a number" */
#ifdef CHECK_NAN_VECTORS  /* we need to check if vx becomes NaN */
	  if(isnan(vx)) {
	    if(is_ion==1){
	      message_error(-100,"error occured (caused by ion no. %i).\n",ion_c);
	    } else {
	      message_error(-100,"Rotation error occured (caused by a recoil in cascade from ion no. %i).\n",ion_c);
	    }
	    /* Print some detailed information: */
	    message_error(-100,"Info:\n old vx: %f, new vx: %f, energy: %f, reduced ip: %f\n,cos(theta)=%f, sin(theta)=%f, matrix index:%i\n",old_vx,vx,energy,red_impact_par,ScatMatrix->CosScat[smIndex],ScatMatrix->SinScat[smIndex],smIndex);
	    message_error(-100,"The projectile is discarded.\n");
	    return -100;
	  }
#endif

	  if(simulation_type<3){ /* We need to consider the recoil */
	    /* Obtain relevant energy barriers */
	    e_disp=(current_material->ElementsDispEnergy)[target_index];
	    e_latt=(current_material->ElementsLattEnergy)[target_index];
	    e_repl=(current_material->ElementsReplEnergy)[target_index];
	    /* e_surf=(current_material->ElementsSurfEnergy)[target_index]; */

	    /* If all of the following conditions are fulfilled, we need to check the surface binding energy instead of the displacement energy:
	       - surface sputtering is considered
	       - the recoil is a surface atom
	       - the flying direction of the recoil points into vacuum
	       - the directional fraction of the recoil energy pointing perpendicular into the vacuum is larger than the surface binding energy
	    */
	    recoil_is_at_surface=0;
	    E_compare=e_disp;
	    surface_sputtered=0;        /* default: not sputtered from surface */
	    if(detailed_sputtering==1){ /* if sputtering is switched on, we need to check if the
					   recoil is a surface atom, and in which direction the vacuum is */
#ifdef INCLUDE_SPECIAL_GEOMETRY
	      if(special_geometry==1){
		recoil_is_at_surface=special_CheckSurfaceAtom(target_index,x,y,z,current_material->LayerDistance);
	      } else {
		recoil_is_at_surface=CheckSurfaceAtom(target_index,x,y,z,current_material->LayerDistance);
	      }
#else
	      recoil_is_at_surface=CheckSurfaceAtom(cell_i,x,y,z,current_material->LayerDistance);
#endif
	      if(recoil_is_at_surface>0){
		/* printf("RIAS: x:\t%f   y:\t%f   z:\t%f \t %i\n",x,y,z,recoil_is_at_surface);
		   printf("TVN %i \n",(int)TargetVacuumNeighbors[cell]); */

		/* Calculate recoil direction: */
		temp1=ScatMatrix->sqrt_mass_ratio*sqrt(  (energy+recoil_energy)/recoil_energy );
		temp2=ScatMatrix->sqrt_mass_ratio*sqrt(  energy/recoil_energy );
		recoil_vx = temp1 * old_vx - temp2 * vx;
		recoil_vy = temp1 * old_vy - temp2 * vy;
		recoil_vz = temp1 * old_vz - temp2 * vz;
		/* (regarding the comment in corteo code: A mon avis le rapport est bon) */
								
#ifdef INCLUDE_SPECIAL_GEOMETRY
		if(special_geometry==0){
#endif
		  /* Check in which direction the vacuum is and if the velocity points into the vacuum: */
		  if( ((recoil_is_at_surface==1)&&(  (E_frac=recoil_vx) > 0)) ||
		      ((recoil_is_at_surface==2)&&(  (E_frac=recoil_vx) < 0)) ||
		      ((recoil_is_at_surface==4)&&(  (E_frac=recoil_vy) > 0)) ||
		      ((recoil_is_at_surface==8)&&(  (E_frac=recoil_vy) < 0)) ||
		      ((recoil_is_at_surface==16)&&( (E_frac=recoil_vz) > 0)) ||
		      ((recoil_is_at_surface==32)&&( (E_frac=recoil_vz) < 0)) ){ /* flying direction points into vacuum */
		    E_frac = E_frac*E_frac*recoil_energy; /* fraction of energy perpendicular to surface */
		    e_surf = (current_material->ElementsSurfEnergy)[target_index]; /* Get surface energy */
		    if(E_frac>e_surf){ /* if relevant energy is larger than surface binding energy, then it is sputtered (simpler sputter counting) */
		      surface_sputtered=1;
		      E_compare=E_frac;
		    }
		  }
#ifdef INCLUDE_SPECIAL_GEOMETRY
		} else { /* We must consider the special geometry */
		  E_frac = special_CalcDirectionalFractionSqr(x,y,z,recoil_vx,recoil_vy,recoil_vz); /* This depends on the geometry.
												       E_frac is already squared! */
		  E_frac = E_frac * recoil_energy; /* fraction of energy perpendicular to surface */
		  e_surf = (current_material->ElementsSurfEnergy)[target_index]; /* Get surface energy */
		  if(E_frac>e_surf){ /* if relevant energy is larger than surface binding energy, then it is sputtered. */
		    surface_sputtered=1;
		    E_compare=E_frac;
		  }
		}
#endif
								
	      } /* else: recoil is not at surface */
	    } /* End of: if detailed sputtering */
						
	    if(recoil_energy>=E_compare){ /* The recoil is displaced from its lattice site */
	      disp_c+=1;

	      /* Increment displacement counter for that element and that cell */
	      ((current_material->TargetElementalDisp)[target_index])[cell_i]+=1;
							
	      if(surface_sputtered==1){ /* has been sputtered from surface */

		/* Since the projectile should not suffer more energy loss etc., we should directly advance it a little bit
		   and move it into vacuum! (by about the spacing!) However this might move it outside target which could be a problem!*/
		/* Start recoil as new projectile: */
		recoil_energy-=E_compare; /* substract surface binding energy */
		FastProjectileTransport(target_Z, target_mass, recoil_energy,
					x+recoil_vx*current_material->LayerDistance,
					y+recoil_vy*current_material->LayerDistance,
					z+recoil_vz*current_material->LayerDistance,
					recoil_vx, recoil_vy, recoil_vz, 0,current_material_index,target_index,cell_i,RD+1); 
								
	      } else { /* displaced from bulk */
		/* we may not know the recoil direction yet, so calculate: */
		/* Calculate flying direction of recoil: */
		/* The formula can be obtained by considering conservation of momentum, and using
		   known energies from projectile before, after and from recoil after collision: */
		temp1=ScatMatrix->sqrt_mass_ratio*sqrt(  (energy+recoil_energy)/recoil_energy );
		temp2=ScatMatrix->sqrt_mass_ratio*sqrt(  energy/recoil_energy );
		recoil_vx=temp1 * old_vx - temp2 * vx;
		recoil_vy=temp1 * old_vy - temp2 * vy;
		recoil_vz=temp1 * old_vz - temp2 * vz;
		/* (regarding the comment in corteo code: A mon avis le rapport est bon) */
								
		/* Energy of recoil has to be reduced by lattice binding energy: */
		recoil_energy-=e_latt;
		if(store_energy_deposit==1){TargetEnergyPhonons[cell_i]+=(double)e_latt;} /* Energy goes to phonons */

		/* Follow recoil as new projectile: */
		FastProjectileTransport(target_Z, target_mass, recoil_energy, x, y, z,
					recoil_vx, recoil_vy, recoil_vz, 0,current_material_index,target_index,cell_i,RD+1);
								
	      } /* end of IF surface_sputtered */

	      replaced=0; /* Assume first, that no replacement occurs */

	      /* Now check whether the projectile might replace the recoil */
	      if(proj_eq_target==1){ /* this only happens if they are the same, otherwise projectile can only become an interstitial if stopped */
		// E_compare=e_disp;    OLD VERSION
		//	if(surface_sputtered==1){
		//	E_compare=e_surf;
		//}
		E_compare=e_repl; /* Flexible replacement energy introduced in version 1.1.1 of iradina. */
		if(energy<E_compare){ /* the projectile cannot leave the site and replaces the recoil, no vacancy! */
		  repl_c+=1;
		  replaced=1;
		  if(is_ion==0){   /* Old recoil (atom from target) replaces new recoil */
		    (ListOfMaterials[OrgMaterial].TargetImplantedRecoilsRepl[OrgElement])[cell_i]+=1;
		  } else { /* The ion replaces the recoil */
		    TargetReplacingIons[cell_i]+=1;
		  }
		  if(store_energy_deposit==1){TargetEnergyPhonons[cell_i]+=(double)energy;} /* remaining projectile energy is released to phonons */
		  energy=-.001; /* set the energy slightly negative to stop the ion. */
		} else { /* The projectile has enough energy to leave site, a vacancy is created */
		  ((current_material->TargetElementalVacancies)[target_index])[cell_i]+=1;
		  /* CROC  store positions of vacancies*/
		  if(store_range3d==2){fprintf(store_range3dV_fp,"%g\t%g\t%g\n",x,y,z);} 
		  /*		  ((current_material->TargetElementalVacancies)[target_index])[cell_i]+=1;*/
		}
	      } else { /* projectile and target not the same */
		/* A vacancy is created in the target */
		((current_material->TargetElementalVacancies)[target_index])[cell_i]+=1;
		/* CROC  store positions of vacancy: */
		if(store_range3d==2){fprintf(store_range3dV_fp,"%g\t%g\t%g\n",x,y,z);} 
	      }
	    } else { /* Recoil cannot be displaced, because e<e_comp, its energy goes into phonons (unless sputtered) */

// TEST TEST
	      
	    // if(proj_eq_target==1){ /* this only happens if they are the same, otherwise projectile can only become an interstitial if stopped */
		
			// if(energy<e_repl){
				// if (recoil_energy<e_disp) {
					// if (energy<recoil_energy) {
		      
		  // /* the projectile cannot leave the site and replaces the recoil, no vacancy! */
		      		      // // repl_c+=1;
		      // // replaced=1;
		      // disp_c+=1;
		      // ((current_material->TargetElementalDisp)[target_index])[cell_i]+=1;
		      // if(is_ion==0){   /* Old recoil (atom from target) replaces new recoil */
			// (ListOfMaterials[OrgMaterial].TargetImplantedRecoilsRepl[OrgElement])[cell_i]+=1;
		      // } else { /* The ion replaces the recoil */
			// TargetReplacingIons[cell_i]+=1;
		      // }

		  
					// }
				// }
			// }
	    // }
     // TESTEST	      
	      
	      if(store_energy_deposit==1){TargetEnergyPhonons[cell_i]+=(double)recoil_energy;}
	    }
	  } else { /* no recoils considered as cascade, but perhaps as KP: */
	    if(simulation_type==5){    /*CROC : KP INCLUDED HERE.  based on average material <> SRIM but FOLLOWS page 7-28 of SRIM book by ZBZ*/
	      E_div=2.5*current_material->MeanEd;
	      e_d = current_material->ed_oE * recoil_energy;
	      g_ed = 3.4008 * pow (e_d, 1.0 / 6.0) + 0.40244 * pow (e_d, 0.75) + e_d;
	      E_v = recoil_energy/ (1.0 + current_material->k_d * g_ed);
	      if(store_energy_deposit==1){
		cell_i=GetCellIndex(x,y,z);
		TargetEnergyElectrons[cell_i]= TargetEnergyElectrons[cell_i]+recoil_energy-E_v;
		TargetEnergyPhonons[cell_i]+= E_v;
	      }
	      if (E_v < current_material->MeanEd) {
		((current_material->TargetElementalVacancies)[0])[cell_i]+=0;
	      }
	      else if (E_v >= current_material->MeanEd && E_v < E_div) {
		((current_material->TargetElementalVacancies)[0])[cell_i]+=1;
	      }
	      else if (E_v >= E_div) {
		((current_material->TargetElementalVacancies)[0])[cell_i]+=floor(E_v / E_div) ;
	      }
	    }
	    if(simulation_type==4){  /*CROC : KP INCLUDED HERE. pure elemental solid after collisions apparently corresponding to  SRIM (page 7-28 book by ZBZ*/
	      e_disp=(current_material->ElementsDispEnergy)[target_index];
	      if (recoil_energy>e_disp) {
	      E_div=2.5*e_disp;
	      k_dr=0.1334 * pow ( target_Z, 2.0 / 3.0) / pow ( target_mass, 0.5);
	      e_d = 0.01014 * pow (target_Z , -7.0 / 3.0)  * recoil_energy;
	      g_ed = 3.4008 * pow (e_d, 1.0 / 6.0) + 0.40244 * pow (e_d, 0.75) + e_d;
	      E_v = recoil_energy/ (1.0 + k_dr * g_ed);
	      /* CROC outputs the ballistic energy and the electronic energy of the cascade as estimated by the above formulas */	      
	      if(store_energy_deposit==1){
		cell_i=GetCellIndex(x,y,z);
		//		printf(" 1 %f %f \n", TargetEnergyElectrons[cell_i], TargetEnergyPhonons[cell_i] );
		TargetEnergyElectrons[cell_i]= TargetEnergyElectrons[cell_i]+recoil_energy-E_v;
		TargetEnergyPhonons[cell_i]+= E_v;
		//		printf(" 2 %f %f \n", TargetEnergyElectrons[cell_i], TargetEnergyPhonons[cell_i]) ;
		//		printf(" %f  %f %i \n ", recoil_energy,E_v,cell_i);
	      }
	      if (E_v < e_disp) {
		((current_material->TargetElementalVacancies)[target_index])[cell_i]+=0;
	      }
	      else if (E_v >= e_disp && E_v < E_div) {
		((current_material->TargetElementalVacancies)[target_index])[cell_i]+=1;
	      }
	      else if (E_v >= E_div) {
		((current_material->TargetElementalVacancies)[target_index])[cell_i]+=floor(E_v / E_div) ;
	      }
	    }
	    }
	  } /* end of: sim type >3, perhaps KP included. */
	} /* Collision took place */
				
	/* Check what happens to the projectile after possible collision: */
	if(energy<min_energy){ /* projectile has to stop. Store as implanted ion or recoil */
	  /*CROC Interstitials and vacancies  are stored in store RANGE 3D V/I  and not only the ions */
	  if(is_ion==1){ /* the ion comes to rest */
	    if(replaced==0){
	      if(store_range3d==2){fprintf(store_range3dI_fp,"%g\t%g\t%g\n",x,y,z);} 
	    }
	    if(store_range3d>=1){fprintf(store_range3d_fp,"%g\t%g\t%g\n",x,y,z);}
	    TargetImplantedIons[cell_i]+=1;
	  } else { /* a target atom stopped */
	    if(replaced==0){ /* if it wasn't a replacement, it becomes interstitial, store as such */
	      if(store_range3d==2){fprintf(store_range3dI_fp,"%g\t%g\t%g\n",x,y,z);} 
	      ((ListOfMaterials[OrgMaterial]).TargetImplantedRecoilsInt[OrgElement])[cell_i]+=1;
	    } /* else: it was a replacement, so we do not need to store interstitial */
	    if(store_recoil_cascades==1){ /* if cascades are to be stored, store empty line in file to separate record */
	      fprintf(recoil_cascades_fp,"\n");
	    }
	  }
	  if(store_energy_deposit==1){TargetEnergyPhonons[cell_i]+=(double)energy;} /* Assume that remaining energy is lost to phonons */
	  energy=-0.001; /* Take away remainig energy to make projectile stop completely. */
	} /* else: Enough energy to advance to next collision site */
      } /* else: There are no collisions in vacuum */
    } /* end of (IF projectile is inside sample) */
  } /* end of energy>0 loop */
  return 0;
}



/**************************************** Alternate completely rewritten transport function               *******/
/**************************************** This is more accurate but slower than the corteo computing flow *******/

int FullProjectileTransport(int ProjZ, float ProjM, double ProjE, float Proj_x, float Proj_y, float Proj_z, float Proj_vx, float Proj_vy, float Proj_vz, int is_ion, int OrgMaterial, int OrgElement, int OrgCell, char ProjState, int RD){
  /* What this function does:
     It performs a loop with:
     - Fly some free flight path 
     - Calculate electronic stopping along the path
     - Check if we have crossed boundary into vacuum or vice versa
     - calculate new position (obey boundary conditions)...
     - simulate a collision
     - if recoil gets enough energy:
     -   store a vacancy
     -   simulate the recoil recursively
     - if projectile has low energy:
     -   stop it and store it as ion or interstitial
     - else: advance to next collision site
  */

  /* projectile properties: */
  double energy;                        /* energy */
  float x,y,z;                          /* position */
  float old_x,old_y,old_z;              /* last position */
  float vx,vy,vz;                       /* velocity unit vector (directional cosines) */
  float old_vx,old_vy,old_vz;           /* last vel unit vector */
  float flightpath=0.0;                 /* free flight path length */
  float flightpath_sqr=0.0;             /* squareroot of flight path */
  float stopping;                       /* electronic stopping */
  float straggling;                     /* electronic energy loss straggling */
  int ion_left_target;                  /* indicates if projectile has left sim volume */
  int isvac;                            /* to store if a material is vacuum */
  int isvac2;                           /* for new material */

  /* materials and positioning: */
  int   cell_i;                         /* index of current cell */
  int   new_cell=0;                     /* cell index at new pos */
  int   new_material_index=0;           /* index of material at new pos */
  int   current_material_index;         /* index of current material */
  struct material* current_material;    /* pointer to current mat */
  struct material* new_material;        /* pointer to that */
  float nx,ny,nz;                       /* surface normal for cell border crossings */

  /* target nucleus properties */
  float targetM=0.0;                    /* mass */
  int   targetZ=0;                      /* protons */
  int   targetIndex=0;                  /* elemental index within material */
  float target_x,target_y,target_z;     /* Position of target nucleus before collision */
  int   target_cell;                    /* index of cell at target pos */
  int   target_material_index=0;        /* index of material at target pos */
  struct material* target_material;     /* pointer to target material */
  int   target_outside_simulation;      /* is 0, when target atom is within simulation volume */
  float recoil_vx,recoil_vy,recoil_vz;  /* recoil velocity unit vector */
  float recoil_energy;                  /* energy transferred to recoil */
  char  next_ProjState=0;               /* Status flags for newly generated recoil */

  /* collisional stuff: */
  float max_impact_par=0.0;             /* maximum impact parameter (IP) */
  float impact_par=0.0;                 /* actual IP */
  //  float red_energy;                 /* reduced projectile energy */
  float red_impact_par;                 /* reduced IP */
  struct scattering_matrix* ScatMatrix; /* pointer to matrix with scattering data */
  int   matrix_index=0;                 /* index of ScatBase matrix for red_energy and red_IP */
  float sin2thetaby2;                   /* sin^2(theta/2), where theta is the scattering angle */
  int   proj_eq_target;                 /* denotes whether projectile and target nucleas are of the same type */
  float theta,sintheta,costheta;        /* scattering angle and stuff */

  /* Other stuff: */
  float temp1,temp2;                    /* temporary floats */
  int   itemp1,i;                       /* temporary int */
  int   replaced=0;                     /* 1, if projectile replaced recoil */
  int   replace_cell=0;                 /* cell index, where replacement occured */
  float vel,svel;                       /* for velocity normalization */
  int   coll_counter;                   /* if multiple collisions are allowed for one flightpath, they need to be counted */
  long int recoil_number;               /* Identifying unique number of this recoil */
  float E_surf, E_disp, E_latt,E_repl;  /* Energy barriers: surface binding, perpendicular energy component, displacement threshold, lattice binding energy, and replacement threshold */
  int   leaving_direction=0;            /* denotes into which direction a particle left the sim volume */
  float random;                         /* a random number */
  float conc_sum;                       /* sum of concentration of elements */
  struct transmitted_ion* leave_p;      /* Pointer to a projectile leaving the simulation volume */

#ifdef PROJ_HANGUP_SAFETY
  int proj_step_counter;
#endif

  /* init some values: */
  recoil_number=recoil_counter++; /* set unique ID of this recoil */
  new_material=ListOfMaterials;
  energy=ProjE;
  x=Proj_x;
  y=Proj_y;
  z=Proj_z;
  vx=Proj_vx;
  vy=Proj_vy;
  vz=Proj_vz;
  cell_i=GetCellIndex(x,y,z);
#ifdef INCLUDE_SPECIAL_GEOMETRY
  if(special_geometry==0){
#endif
    current_material_index=TargetComposition[cell_i];
#ifdef INCLUDE_SPECIAL_GEOMETRY
  } else {
    current_material_index=GetMaterialFromPosition(cell_i,x,y,z);
  }
#endif

  if(RD==1){ /* PKA! */
    if(store_PKA>=1){ /* if requested, store each PKA with type, positions, direction, energy, ...!) */
      fprintf(store_PKA_fp,"%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%f\t%g\n",x,y,z,vx,vy,vz,OrgMaterial,OrgElement,ProjZ,ProjM,ProjE);
    }
  }

  current_material=&(ListOfMaterials[current_material_index]);

#ifdef PROJ_HANGUP_SAFETY
  proj_step_counter=0;
#endif

  while(energy>0){ /* let projectile proceed until all its energy is lost */

    /* if((ion_c>MONITOR_ION)&&(1)){
       printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);
       printf(" Collisions:      %i\n",proj_step_counter);
       printf(" Ion?\t %i\n",is_ion);
       printf(" Velocity vector: %g\t%g\t%g\n",vx,vy,vz);
       printf(" Position:        %g\t%g\t%g\n",x,y,z);
       printf(" Energy:          %g\n",energy);
       printf(" Ion number:      %i\n",ion_c);
       }*/
		
		
#ifdef PROJ_HANGUP_SAFETY
    if((proj_step_counter++)>PROJ_HANGUP_SAFETY){
      //if(is_ion==1){ /* To do: take out this line! */
      printf("Attention: projectile discarded after %i steps (possible hang-up situation).\n",proj_step_counter);
      printf(" Ion?\t %i\n",is_ion);
      printf(" Velocity vector: %g\t%g\t%g\n",vx,vy,vz);
      printf(" Position:        %g\t%g\t%g\n",x,y,z);
      printf(" Ion number:      %i\n",ion_c);
      //}
      return 0;
    }
#endif
		
    if(is_ion==1){
      if(store_ion_paths==1){ /* if requested, store position of ion */
	fprintf(ion_paths_fp,"%g\t%g\t%g\t%g\n",x,y,z,energy);
      }
    } else { /* Recoil */
      if(store_recoil_cascades==1){ /* if requested, store position of recoil */
	fprintf(recoil_cascades_fp,"%g\t%g\t%g\t%g\t%li\t%i\n",x,y,z,energy,recoil_number,RD);
      }
    }

    isvac=current_material->Is_Vacuum;

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g\n",x,y,z);
#endif

    /* Calculate flight path (depends on material we're in) */
    if(isvac==1){ /* vacuum */
      /* The free flight path in vacuum should not be constant, because this can cause artifacts in sputter yield calculations.*/
      /* So, some random variations must be used (this is new in version 1.0.1 of iradina!). */
      flightpath=0.15+randomlist[(iranlist)++]*0.14; 
      flightpath_sqr=sqrtf(flightpath);
      if(iranlist>=MAXRANLIST-2)iranlist = 0; /* reset list pointer if needed */
      stopping=0.0;                           /* no stopping or straggling in vacuum */
      straggling=0.0;
    } else { /* no vacuum */
      /* Calculate the flight path, depending on material: */
      switch(flight_length_type){
      case 0: /* Poisson distributed flight length and impact pars */
	/* Flight length between collisions: Poisson-distributed value with interatomic spacing as mean free path */
	/* flightpath = -current_material->AtomicDistance * log(randomx()); */
	flightpath_sqr = sqrtloglist[iranloglist]*current_material->SqrtAtomicDistance; /* Use of random list faster (as in corteo)*/
	flightpath = flightpath_sqr*flightpath_sqr;
	break;
      case 1: /* atomic spacing */
	flightpath=current_material->AtomicDistance;
	flightpath_sqr=sqrtf(flightpath);
	break;
      case 2: /* constant flightpath in nm*/
	flightpath=flight_length_constant;
	flightpath_sqr=sqrtf(flightpath);
	break;
      }

#ifdef DEBUG_MODE
      message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g\n",x,y,z);
#endif

      /* calculate electronic energy loss (+- straggling): */
      stopping   = (flightpath * 10.0) * ElectronicStopping(ProjZ,ProjM,energy,current_material_index);
      /* factor 10, because flightpath is in nm, stopping in A */
      straggling = 3.16227766016838 * flightpath_sqr *  ElectronicStraggling(ProjZ,ProjM,energy,current_material_index) * inverse_erf_list[erflist_pointer++];
      /* Due to Gaussian distribution, the straggling can in some cases get so large that the projectile gains energy or suddenly looses a huge amount of energy. Is this realistic? This might actually happen. However, in the simulation, ions may have higher energy than the initial energy. We will therefore limit the |straggling| to |stopping|.  Furthermore, with hydrogen, the straggling is often so big, that ions gain huge amount of energy, and the phononic system would actually gain energy. */
      if(fabs(straggling)>stopping){
	if(straggling<0){
	  straggling=-stopping;
	} else {
	  straggling=stopping;
	}
      }

      /* The stopping tables have no values below 16 eV. Therefore, we do simple linear downscaling of electronic stopping below 16 eV. */
      if(energy<16){ /* scaling down */
	stopping*=(energy*0.0625);
	straggling*=(energy*0.0625);
      }
      if(erflist_pointer>=MAXERFLIST){erflist_pointer=0;} /* Check and adjust boundary */
    }

#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"Stop: %g, Stragg: %g\n",stopping,straggling);
    message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g\n",x,y,z);
    message_debug(__FILE__,__LINE__,"proj vx vy vz: %g %g %g\n",vx,vy,vz);
    message_debug(__FILE__,__LINE__,"flightpath: %g\n",flightpath);
#endif

    /* Calculate end position of flight path: */
    old_x=x; old_y=y; old_z=z;
    x += flightpath*vx;
    y += flightpath*vy;
    z += flightpath*vz;
#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g\n",x,y,z);
#endif

    /* Now check if projectile is still within the simulation volume */
    /* If not, check if periodic boundary conditions must be applied */
    ion_left_target=0;
    ion_left_target+=CheckAndCorrectBoundary(&x,target_size_x,&target_max_x,boundary_x);
    ion_left_target+=CheckAndCorrectBoundary(&y,target_size_y,&target_max_y,boundary_y);
    ion_left_target+=CheckAndCorrectBoundary(&z,target_size_z,&target_max_z,boundary_z);
#ifdef DEBUG_MODE
    message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g, ion_left_target: %i\n",x,y,z,ion_left_target);
#endif

    if(ion_left_target!=0){ /* projectile left sample! */

      /* Outside the simulation there might be vacuum...! If the last position was material,
	 then we should make sure that the projectile has sufficient energy to overcome
	 the surface binding energy, otherwise it should be refleced. This is not implemented.
	 Instead, the user should make sure that there is some vacuum around the structure,
	 if this is required. */
      /*  if it is the ion that has left the target, store transmitted ion: */
      if(is_ion==1){
	if(store_transmitted_ions==1){
	  transmit_list[transmission_pointer].x=x;
	  transmit_list[transmission_pointer].y=y;
	  transmit_list[transmission_pointer].z=z;
	  transmit_list[transmission_pointer].vx=vx;
	  transmit_list[transmission_pointer].vy=vy;
	  transmit_list[transmission_pointer].vz=vz;
	  transmit_list[transmission_pointer].energy=energy;
	  transmission_pointer++;
	}
	if(detailed_sputtering==1){ /* if detailed sputtering is on, we also take a look on the ions leaving the target */
	  if(x>target_size_x){leaving_direction=0;}
	  if(x<=0){leaving_direction=1;}
	  if(y>target_size_y){leaving_direction=2;}
	  if(y<=0){leaving_direction=3;}
	  if(z>target_size_z){leaving_direction=4;}
	  if(z<=0){leaving_direction=5;}
	  leaving_ions[leaving_direction]+=1;
	}
	leaving_ions_c++;
      } else { /* its a recoil leaving the simulation volume! */
	/* If the recoil was subthreshold, we will reset it to its original site. That means it neither creates a vac, nor an int, nor is it leaving the sim volume */
	if( (ProjState&1)==0){ /* Let particle leave only if it was a real (super threshold) recoil */
	  if(store_exiting_recoils==1){ 
	    itemp1  = (ListOfMaterials[OrgMaterial]).leaving_recoils_pointer[OrgElement];
	    leave_p = &(((ListOfMaterials[OrgMaterial]).ElementalLeavingRecoils[OrgElement])[itemp1]);
	    if(itemp1<store_exiting_limit){ /* not too many recoils yet */
	      leave_p->x = x;
	      leave_p->y = y;
	      leave_p->z = z;
	      leave_p->vx = vx;
	      leave_p->vy = vy;
	      leave_p->vz = vz;
	      leave_p->energy = energy;
	      (ListOfMaterials[OrgMaterial]).leaving_recoils_pointer[OrgElement]+=1;
	    }
	  }
					
	  isvac=current_material->Is_Vacuum;

#ifdef DEBUG_MODE
	  message_debug(__FILE__,__LINE__,"\n");
#endif
					
	  /* Simple sputter counting: a sputtered atom is one that went through vacuum before it left the sample */
	  if(isvac==1){
	    sputter_c+=1;
	  }
	  leaving_recoils_c++;

	  /* More complex sputter counting (only takes place, if detailed sputtering is switched on:) */
	  /* for each target material and element, store the number of atoms ejected out of the simulation
	     volume into each of the possible 6 directions. */
	  if(detailed_sputtering==1){
	    single_ion_sputter_counter++;
	    /* Determine direction of leaving out of target: */
	    leaving_direction=0;  /* This is necessary. In some extremely rare cases, leaving direction is undefined otherwise --> crash! */
	    /* if(x>target_size_x){leaving_direction=0;} */
	    if(x<=0){leaving_direction=1;}
	    if(y>target_size_y){leaving_direction=2;}
	    if(y<=0){leaving_direction=3;}
	    if(z>target_size_z){leaving_direction=4;}
	    if(z<=0){leaving_direction=5;}
	    ListOfMaterials[OrgMaterial].SputterCounter[(OrgElement*6)+leaving_direction]++;
	    (ListOfMaterials[OrgMaterial].TargetSputteredAtoms[OrgElement])[OrgCell]++; /* With this, we can later see where particles are sputtered from */
	  }
	  if(store_recoil_cascades==1){ /* if cascades are to be stored, store empty line into file to separate record */
	    fprintf(recoil_cascades_fp,"\n");
	  }
	} else { /* if( (ProjState&1)==1), meaning it was a sub threshold recoil */
	  /* This recoil will not be allowed to leave the target */
	  /* The remaining energy of this recoil must be deposited at its original cell, because it is assumed to jump back */
	  if((ProjState&2)==0){ /* it was not a replacer. Its energy goes to its origin: */
	    if(store_energy_deposit==1){TargetEnergyPhonons[OrgCell]+=(double)energy;} /* remaining projectile energy is released to phonons */
	  } else { /* It was a replacer. Its energy goes to the replacement cell */
	    if(store_energy_deposit==1){TargetEnergyPhonons[replace_cell]+=(double)energy;} /* remaining projectile energy is released to phonons */	    
	  }
	} /* end of: ProjState LSB was 0 or 1 */
      } /* end of: it's a recoil leaving the sim volume */

      energy=-0.001; /* in order to exit the loop and stop the particle's movement.*/
			
    } else { /* Projectile is still within sample: */
#ifdef DEBUG_MODE
      message_debug(__FILE__,__LINE__,"proj xzy: %g %g %g\n",x,y,z);
#endif

      /* Determine material at new position */
      new_cell=GetCellIndex(x,y,z);
#ifdef DEBUG_MODE
      message_debug(__FILE__,__LINE__,"new cell index: %i\n",new_cell);
#endif

#ifdef INCLUDE_SPECIAL_GEOMETRY
      if(special_geometry==0){
#endif
				
	new_material_index=TargetComposition[new_cell];
#ifdef INCLUDE_SPECIAL_GEOMETRY
      } else {
	new_material_index=GetMaterialFromPosition(new_cell,x,y,z);
      }
#endif
      new_material=&(ListOfMaterials[new_material_index]);
      //if((ion_c>MONITOR_ION)&&(1)){  printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);}

#ifdef DEBUG_MODE
      message_debug(__FILE__,__LINE__,"new mat index: %i\n",new_material_index);
#endif

      /* Boundary crossing checks; only necessary if detailed sputtering is of interest! */
      if(detailed_sputtering==1){
	isvac  = current_material->Is_Vacuum;
	isvac2 = new_material->Is_Vacuum;
	//  printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);

	if(isvac==1){ /* old material is vacuum */
	  if(isvac2==0){   /* new material is not vacuum */
#ifdef DEBUG_MODE
	    if(ion_c>MONITOR_ION){  printf("DEBUG %s, line %i. Proj tries to cross in to material.\n",__FILE__,__LINE__);}
#endif
	    /* --> Projectile crosses into material */
	    /* determine surface binding energy (SBE) */
	    if(is_ion==1){ /* its the ion. Use ion SBE of new material */
	      E_surf=new_material->IonSurfEnergy;
	    } else {  /* its not the ion. Use SBE of that element (though it may differ
			 in another material...)  */
	      E_surf=(ListOfMaterials[OrgMaterial]).ElementsSurfEnergy[OrgElement];
	    }
#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__, "E_Surf: %g\n",E_surf);
#endif
						
	    /* crossing from vacuum into material: gain SBE and refract direction: */
#ifdef INCLUDE_SPECIAL_GEOMETRY
	    if(special_geometry==0){
#endif
	      /* Attention: when leaving the material, the role of new and old cell must be reverted, because the classical CalcSurfaceNormal returns a vector towards the flying direction */
	      CalcSurfaceNormal(new_cell,cell_i,&nx,&ny,&nz);
	      /*CalcSurfaceNormal(cell_i,new_cell,&nx,&ny,&nz); */ /* would results into energy loss instead of gain */
#ifdef DEBUG_MODE
	      message_debug(__FILE__,__LINE__, "normal vec: %g %g %g.\n",nx,ny,nz);
#endif
#ifdef INCLUDE_SPECIAL_GEOMETRY
	    } else {
	      /* The special version always returns the correct vector out of the material */
	      special_CalcSurfaceNormal(cell_i,new_cell,x,y,z,&nx,&ny,&nz,1.0);
	    }
#endif
						
#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__, "\n");
#endif
						
#ifdef DEBUG_MODE3
	    if(ion_c>MONITOR_ION){
	      message_debug(__FILE__,__LINE__, "into mat: Normal vec: nx=%f, ny=%f, nz=%f\n",nx,ny,nz);
	      message_debug(__FILE__,__LINE__, "into mat: Before refr: E=%f, vx=%f, vy=%f, vz=%f\n",energy,vx,vy,vz);
	      message_debug(__FILE__,__LINE__,"vx: %f\n",vx);
	    }
#endif
						
	    RefractProjectile2(&vx,&vy,&vz,nx,ny,nz,&energy,E_surf);
						
#ifdef DEBUG_MODE3
	    if(ion_c>MONITOR_ION){
	      message_debug(__FILE__,__LINE__, "DEBUG %s, line %i.\n",__FILE__,__LINE__);
	      message_debug(__FILE__,__LINE__, "into mat: After  refr: E=%f, vx=%f, vy=%f, vz=%f\n\n",energy,vx,vy,vz);
	    }
#endif
	    /* subtract fractional electronic stopping in new material */
	    /* we should calculate here the fractional paths, but this takes some time.
	       Instead, we assume so little energy is lost, we ignore it, because
	       flight paths are short (not as in TRIM)! */
	  } else { /* New material is also vacuum, do nothing */
	  }
#ifdef DEBUG_MODE
	  message_debug(__FILE__,__LINE__, "\n");
#endif
					
	} else { /* old material is not vacuum */
	  if(isvac2==0){   /* new material is also not vacuum */
#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__,"\n");
#endif
	    /* projectile stays within materials, no refraction required */
	    /* subtract electronic stopping: */
	    energy -= (stopping+straggling);
	    if(store_energy_deposit==1){
	      TargetEnergyElectrons[new_cell]+=(double)(stopping+straggling);
	    }
	    if(energy<min_energy){ /* energy below cut-off, stop projectile */
	      /* will be done below */
	      if(store_energy_deposit==1){TargetEnergyPhonons[new_cell]+=(double)energy;} /* remaining projectile energy is released to phonons */
	      energy=-.00001; /* set the energy slightly negative to stop the projectile. */
	    }
#ifdef DEBUG_MODE
	    //printf("DEBUG %s, l %i\n",__FILE__,__LINE__);fflush(stdout);
	    message_debug(__FILE__,__LINE__, "After energy check\n");
#endif

	  } else { /* New material is vacuum. That means projectile tries to leave the solid! */
	    /* subtract fractional electronic stopping in old material */
	    /* we should calculate here the fractional paths, but this takes some time.
	       Instead, we assume that full stopping occurs in old cell. This is ok, 
	       because our flight paths are short (not as in TRIM)! */
	    energy -= (stopping+straggling);
	    if(store_energy_deposit==1){
	      TargetEnergyElectrons[cell_i]+=(double)(stopping+straggling);
	    }

	    /* To DO: Hier schon: wenn die Energie jetzt schon kleiner als E_min ist, kann man sofort aufhoeren. Dadurch wird es ggf. schneller. */

	    /* determine surface binding energy: */
	    if(is_ion==1){ /* its the ion. Use ion SBE of new material */
	      E_surf=current_material->IonSurfEnergy;
	    } else {  /* its not the ion. Use SBE of that element (though it may differ in another material...)  */
	      E_surf=(ListOfMaterials[OrgMaterial]).ElementsSurfEnergy[OrgElement];
	    }
#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__,"\n");
#endif
	    /* calculate energy fraction of velocity perp. to surface: */
#ifdef INCLUDE_SPECIAL_GEOMETRY
	    if(special_geometry==0){
#endif
	      CalcSurfaceNormal(cell_i,new_cell,&nx,&ny,&nz);
#ifdef INCLUDE_SPECIAL_GEOMETRY
	    } else {
	      special_CalcSurfaceNormal(cell_i,new_cell,x,y,z,&nx,&ny,&nz,+1.0);
	    }
#endif

#ifdef DEBUG_MODE2
	    if(ion_c>MONITOR_ION){
	      printf("Try escape: Before refr: E=%g, vx=%f, vy=%f, vz=%f\n",energy,vx,vy,vz);
	      printf("Try escape: new pos:     x=%f,  y=%f,  z=%f\n",x,y,z);
	      printf("Try escape: Normal vec: nx=%f, ny=%f, nz=%f\n",nx,ny,nz);
	    }
#endif
	    itemp1 = RefractProjectile2(&vx,&vy,&vz,nx,ny,nz,&energy,E_surf);
#ifdef DEBUG_MODE2
	    if(ion_c>MONITOR_ION){
	      printf("from mat: After  refr: E=%f, vx=%f, vy=%f, vz=%f\n\n",energy,vx,vy,vz);
	    } 
#endif
						
	    if(itemp1==1){ /* projectile escaped */
	      escape_solid_c++;
	      if( (ProjState&1)==1){ 
		/* projectile was sub-threshold! But since it has reached the surface and was able to escape from the solid,
		   it is not assumed to jump back! Therefore, we have to create a vacancy at its original position! */
		if((ProjState&2)==0){ /* the projectile was not a replacer */
		  /* create the original vacancy and disp, which was not created upon displacement. */
		  ((ListOfMaterials[OrgMaterial]).TargetElementalVacancies[OrgElement])[OrgCell]+=1;
		  ((ListOfMaterials[OrgMaterial]).TargetElementalDisp[OrgElement])[OrgCell]+=1;
		  vac_c++;disp_c++;
		} else { /* the sub-threshold-particle was about to replace another atom, this has now become obsolete. Therefore we need to
			    create a vacancy instead: */
		  ((ListOfMaterials[OrgMaterial]).TargetElementalVacancies[OrgElement])[replace_cell]+=1;
		  replacer_escaped++;
		  ProjState^=2; /* make second last bit zero to relabel projectile as a non-replacer.  */
		  replaced=0;
		}
		ProjState^=1; /* make last bit zero to relabel projectile as super-threshold in retrospective */
	      }
	    } else { /* projectile was reflected because normal energy component smaller than SBE. Set back position */
							
	      /* Approximation: */
	      x=old_x; y=old_y; z=old_z;
	      new_cell=cell_i;
	      new_material_index=current_material_index;
	      new_material=current_material;
	      if(itemp1!=2){
		if(special_geometry==0){
		  //printf("Error! line %i. This should not be printed (or only really not so often)! ision: %i,itemp1=%i\n",__LINE__,is_ion,itemp1);
		  /* although this case should not happen it can happen in some rare cases: when a projectile leaves the 
		     solid at an edge or corner of a cell and the velocity is almost perpendicular to one cartesian direction. */
		}
	      }
	    } /* end of: proj escaped a solid or not */
	  } /* end of: proj tries to leave solid into vacuum*/
	} /* end of: old material was not vac */
      } else { /* no detailed sputtering */
	/* simply subtract stopping */
	energy -= (stopping+straggling);
	if(store_energy_deposit==1){
	  TargetEnergyElectrons[new_cell]+=(double)(stopping+straggling);
	}
      }

#ifdef DEBUG_MODE
      message_debug(__FILE__,__LINE__,"\n");
#endif

      /* if projectile has some energy left, it may now collide */
      if(energy>min_energy){ /* projectile has not stopped */
	for(coll_counter=0;coll_counter<=max_annular_coll_volumes;coll_counter++){  /* there may be multiple collisions (see W. Eckstein, p. 92ff) */
	  isvac2=new_material->Is_Vacuum;

	  if(isvac2==0){ /* not in vacuum */
	    /* Select maximum impact parameter - this depends also on the flight length distribution model */
	    switch(flight_length_type){
	    case 0: /* Poisson distributed flight length and impact pars */

	      //	      max_impact_par = current_material->MeanImpactPar * sqrtloglist1[iranloglist];
	      /* There was a bug here in ealier version: current_material was used instead of new_material. This sometimes caused crashes, when projectiles moved from vacuum into material */
	      max_impact_par = new_material->MeanImpactPar * sqrtloglist1[iranloglist];
#ifdef DEBUG_MODE
							
	      message_debug(__FILE__,__LINE__,"mean impact par: %g\n",new_material->MeanImpactPar);
#endif
	      if(++iranloglist>=MAXLOGLIST) iranloglist=0;
	      break;
	    case 1: /* atomic spacing */
	      max_impact_par = new_material->MeanImpactPar; /* without the log list! */
	      break;
	    case 2: /* constant */
	      max_impact_par = new_material->SqrtRecFlDensity;  /* without the log list! */
	      break;
	    }
	  } else { /* In vacuum: what sort of impact parameter should we chose? hm. Difficult, this depends on the material nearby */
	    /* let's first try 0.3 nm as maximum IP. At this distance we might have interaction with surface atoms. */
	    max_impact_par = 0.3;
	  }

	  /* select randomly between 0 and max (or higher, as disucssed in Eckstein's book) */
	  impact_par = (sqrtf(coll_counter) + sqrtrandomlist[(iranlist)++] * sqrtf(coll_counter+1)) * max_impact_par;
	  if(iranlist>=MAXRANLIST-2) iranlist = 0;  /* Make sure index pointer is ok */

#ifdef DEBUG_MODE
	  message_debug(__FILE__,__LINE__,"ip: %g, max ip: %g\n",impact_par,max_impact_par);
#endif

	  /* Calculate position of target nucleus (by a sort of rotation with theta=270 degree): */
	  /* Determine direction first: */
	  CalculateRelativeTargetAtomPosition(vx,vy,vz,&target_x,&target_y,&target_z,iazimAngle);
	  /* Calculate actual position: */
	  target_x = x+(target_x*impact_par);
	  target_y = y+(target_y*impact_par);
	  target_z = z+(target_z*impact_par);
#ifdef DEBUG_MODE
	  message_debug(__FILE__,__LINE__,"t_vec: %g %g %g.\n",target_x,target_y,target_z);
#endif
					
	  /* Check and correct boundary for target position! */
	  target_outside_simulation=0;
	  target_outside_simulation+=CheckAndCorrectBoundary(&target_x,target_size_x,&target_max_x,boundary_x);
	  target_outside_simulation+=CheckAndCorrectBoundary(&target_y,target_size_y,&target_max_y,boundary_y);
	  target_outside_simulation+=CheckAndCorrectBoundary(&target_z,target_size_z,&target_max_z,boundary_z);

	  /* SAFETY */
	  /* if(target_x>=target_size_x){printf(" WARNING! x position at edge! x: %g, xedge: %g\n",target_x,target_size_x);}
	     if(target_y>=target_size_y){printf(" WARNING! y position at edge! y: %g, yedge: %g\n",target_y,target_size_y);}
	     if(target_z>=target_size_z){printf(" WARNING! z position at edge! z: %g, zedge: %g\n",target_z,target_size_z);} */
					
#ifdef DEBUG_MODE
	  message_debug(__FILE__,__LINE__,"\n");
#endif
	  if(target_outside_simulation==0){ /* target nucleus is within simulation volume */
	    //	    if(ion_c>=DEBUGIONLIMIT){printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);fflush(stdout);}
	    target_cell=GetCellIndex(target_x,target_y,target_z);

#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__,"\n");
#endif

#ifdef INCLUDE_SPECIAL_GEOMETRY
	    if(special_geometry==0){
#endif
#ifdef DEBUG_MODE
	      message_debug(__FILE__,__LINE__,"TC: %i, t_vec: %g %g %g.\n",target_cell,target_x,target_y,target_z);
#endif
	      target_material_index=TargetComposition[target_cell];
#ifdef INCLUDE_SPECIAL_GEOMETRY
	    } else {
#ifdef DEBUG_MODE
	      message_debug(__FILE__,__LINE__,"TC: %i, t_vec: %g %g %g.\n",target_cell,target_x,target_y,target_z);
#endif
	      target_material_index=GetMaterialFromPosition(target_cell,target_x,target_y,target_z);
	    }
#endif
	    target_material=&(ListOfMaterials[target_material_index]);
	    isvac=target_material->Is_Vacuum;

#ifdef DEBUG_MODE
	    message_debug(__FILE__,__LINE__,"\n");
#endif

	    if(isvac==0){ /* target position is not in vacuum, collision occurs! */

	      /**************** COLLISION ********************************************************/
	      coll_c++;
	      proj_eq_target = 0; /* Assume at first, that projectile is not the same as target */
	      //replaced       = 0; /* Assume at first, that no replacement occured */
	      /* Randomly select collision target atom type according to concentration: */
	      random = randomlist[(iranlist)++]; /* A list of precalculated random values is used like in corteo */
	      conc_sum = 0;

	      for(i=0;i<target_material->ElementCount;i++){    /* go through elements in material... */
		conc_sum+=(target_material->ElementsConc)[i];  /* ... and sum up concentrations */
		if(random<conc_sum){ /* pick this element as target*/
		  targetIndex = i;
		  targetZ  = (target_material->ElementsZ)[targetIndex];
		  targetM  = (target_material->ElementsM)[targetIndex];
		  if(targetZ==ProjZ){ /* check if projectile and target are the same: */
		    if(targetM==ProjM){proj_eq_target=1;} /* Note, it can be dangerous to compare floats for equality, however, the
							     masses are never calculated but always loaded, thus they should be
							     exactly the same in memory for the same target atoms */
		  }
		  break;
		}
	      }

	      /* Now that we have selected the target nucleus, we know which of the scattering_matrices to use */
	      if(is_ion==1){ /* It's the initial ion */
		ScatMatrix=&(ion_scattering_matrix[targetZ]);
	      } else { /* it's some other projectile */
		ScatMatrix=&(scattering_matrices[ProjZ][targetZ]);
	      }

	      /* calculate reduced impact parameter s=p/a */
	      red_impact_par=impact_par * ScatMatrix->inv_screening_length;
	      /* (Note that in contrast to corteo, the screening length and impact par are in units of nm here,
		 so it doesn't matter for the reduced impact parameter) */

	      /* Now calculate the scattering angle */
	      /* This can be done by the corteo database method or by MAGIC. */  
	      /* Already the original ScatBase approach by Yuan is reported to be approx. 20 times faster than MAGIC
		 and 200 times faster the guass-mehler quadrature of the integral (NIMB83, 413 (1993)) 
		 But: MAGIC has advantages also: the scattering matrices need a lot of memory, as every atoms/atom combination
		 of the target needs to be calculated! Furthermore MAGIC offers the advantage, that target isotopes different
		 from the most abundant could be used */


#ifdef DEBUG_MODE
	      message_debug(__FILE__,__LINE__,"\n");
#endif
	      /* if((ion_c>MONITOR_ION)&&(1)){  printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);}
		 if(ion_c>=DEBUGIONLIMIT){printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);fflush(stdout);}
		 if(ion_c>=MONITOR_ION){
		 printf("DEBUG %s, line %i.\n",__FILE__,__LINE__);
		 printf("x: %g, y: %g, z: %g\n",x,y,z);
		 printf("target_x: %g, target_y: %g, target_z: %g\n",target_x,target_y,target_z);
		 printf("vx: %g, vy: %g, vz: %g\n",vx,vy,vz);
		 printf("Target: Z: %i, M: %g\n",targetZ,targetM);
		 printf("ScatMatrix: SL %g, ISL %g,MR %g, MR %g, SMR %g, k %g,%p,%p\n",ScatMatrix->screening_length,ScatMatrix->inv_screening_length,ScatMatrix->mass_ratio,ScatMatrix->sqrt_mass_ratio,ScatMatrix->kfactor_m,ScatMatrix->red_E_conv,ScatMatrix->CosScat,ScatMatrix->SinScat);
		 fflush(stdout);
		 }*/
							
							
	      temp1=10.0*target_material->LayerDistance;
							
	      if(impact_par >= temp1) {  /* impact parameter much larger than interatomic distance: assume collision has missed */
		sin2thetaby2 = 0.0f;
		miss_c+=1;
	      } else { /* Collision actually takes place */
		if(scattering_calculation==0){ /* Use corteo's database method to determine scattering angle. */
		  /* Get matrix index (using corteo's indexing functions and reduced values for E and impact par): */
		  matrix_index = Eindex(energy * ScatMatrix->red_E_conv) * DIMS + Sindex(red_impact_par);
		  sintheta=(ScatMatrix->SinScat)[matrix_index];
		  costheta=(ScatMatrix->CosScat)[matrix_index];
		  sin2thetaby2 = Matrix(matrix_index);
		} else { /* MAGIC */
		  /* Using MAGIC (all angles in CMS): */
		  temp1 = MAGIC(red_impact_par,energy*ScatMatrix->red_E_conv);
		  if(fabs(temp1)>=1){temp1=0.99999;}
		  /*if(abs(temp1)<=0.00001){temp1=0.00001;}*/
		  sin2thetaby2 = 1.0 - (temp1*temp1);
		  costheta=2*temp1*temp1-1;
		  sintheta=sqrt(1.0-costheta*costheta);
		  /* now conversion to lab frame of reference: */
		  theta=atan( sintheta/(costheta+ProjM/targetM));
		  /*if(isnan(theta)){theta=0;}*/
		  sintheta=sin(theta);
		  costheta=cos(theta);
		}
		/* Energy transfer to recoil: */
		recoil_energy = energy * ScatMatrix->kfactor_m * sin2thetaby2;
		energy       -= recoil_energy;

		/* Store old flying direction */
		old_vx = vx;
		old_vy = vy;
		old_vz = vz;
								
		/* Calculate new projectile direction */
		/* This is done by using the rotation routine from the corteo code (equivalent to TRIM's DIRCOS) */
		/* In the original version of corteo, there is a rare case of theta=180, which leads to sinTheta
		   being a "not a number". This happens at extremely small impact parameters. --> has been corrected here */
#ifdef DEBUG_MODE
		message_debug(__FILE__,__LINE__,"proj vx vy vz: %g %g %g\n",vx,vy,vz);
#endif
		rotate(&vx,&vy,&vz,&iazimAngle,costheta,sintheta);
#ifdef DEBUG_MODE
		message_debug(__FILE__,__LINE__,"proj vx vy vz: %g %g %g\n",vx,vy,vz);
#endif

		/* Possible safety factor: */
#ifdef RENORM_VELOCITIES
		vel=vx*vx+vy*vy+vz*vz;
		if(fabs(vel)>1.5){
		  //printf("Strange velocity warning: %f\n",vel);
		  svel=1.0/sqrtf(vel);
		  vx*=svel;vy*=svel;vz*=svel;
		}
#endif
								
#ifdef CHECK_NAN_VECTORS  /* we need to check if vx becomes NaN */
		if(isnan(vx)) {
		  if(is_ion==1){
		    message_error(-100,"Rotation error occured (caused by ion no. %i).\n",ion_c);
		  } else {
		    message_error(-100,"Rotation error occured (caused by a recoil in cascade from ion no. %i).\n",ion_c);
		  }
		  /* Print some detailed information: */
		  message_error(-100,"Info:\n old vx: %f, new vx: %f, energy: %f, reduced ip: %f\n,cos(theta)=%f, sin(theta)=%f, matrix index:%i\n",old_vx,vx,energy,red_impact_par,ScatMatrix->CosScat[matrix_index],ScatMatrix->SinScat[matrix_index],matrix_index);
		  message_error(-100,"The projectile is discarded.\n");
		  return -100;
		}
#endif
								
		/* Check what happens to target nucleus: */

		/* Determine lattice binding energy and displacement energy: */
		E_latt = target_material->ElementsLattEnergy[targetIndex];
		E_disp = target_material->ElementsDispEnergy[targetIndex];
		//		if(ion_c>=DEBUGIONLIMIT){printf("DE TR, l %i i%i\n",__LINE__,ion_c);fflush(stdout);}

		if((recoil_energy>min_energy)&&(recoil_energy>E_latt)){ /* recoil gains at least the minimum energy, it might be displaced */
		  next_ProjState=0;
		  if(recoil_energy<E_disp){ /* Mark recoil as sub-threshold, because E_Rec < E_Disp */
		    next_ProjState=1;
		  }
		  /* Recoil starts as new projectile */
		  /* Calc recoil direction: */
		  temp1=ScatMatrix->sqrt_mass_ratio*sqrt(  (energy+recoil_energy)/recoil_energy );
		  temp2=ScatMatrix->sqrt_mass_ratio*sqrt(  energy/recoil_energy );
		  recoil_vx=temp1 * old_vx - temp2 * vx;
		  recoil_vy=temp1 * old_vy - temp2 * vy;
		  recoil_vz=temp1 * old_vz - temp2 * vz;
		  /* (regarding the comment in corteo code: A mon avis le rapport est bon) */
		  //		  if(ion_c>=DEBUGIONLIMIT){printf("DE TR, l %i i%i\n",__LINE__,ion_c);fflush(stdout);}
									
		  /* Subtract lattice binding energy */
		  recoil_energy-=E_latt;

		  /* Recursive call to start recoil as new projectile: */
		  FullProjectileTransport(targetZ, targetM, recoil_energy,
					  target_x,target_y,target_z, recoil_vx, recoil_vy, recoil_vz,
					  0,target_material_index,targetIndex,target_cell,next_ProjState,RD+1); 
									
		  /* Update counters: */
		  if(next_ProjState==0){ /* really displaced! (Energy was larger than displacement energy) */
		    /* create displacement and vacancy: */
		    ((target_material->TargetElementalDisp[targetIndex])[target_cell])++;
		    ((target_material->TargetElementalVacancies[targetIndex])[target_cell])++;
		    disp_c++;vac_c++;
		    /* check if projectile might replace the recoil. Conditions for replacement collision:
		       1. equal particles, 2. proj has less than lattice energy, 3. proj was not marked sub-threshold, 4. recoil is not sub-threshold */
		    if( (proj_eq_target==1)){
		      E_repl = target_material->ElementsReplEnergy[targetIndex];
		      if(energy<E_repl){  /* changed in version 1.1.1 of iradina: compare to e_repl instead of e_disp */
			if((ProjState&1)==0){
			  /* all conditions fulfilled, replace! */
			  if(is_ion==0){ /* replacing recoil: */
			    replaced=1;
			    ((ListOfMaterials[OrgMaterial].TargetImplantedRecoilsRepl[OrgElement])[target_cell])++;
			    repl_c++;
			    /* Remove vacancy: */
			    ((target_material->TargetElementalVacancies[targetIndex])[target_cell])--;
			    vac_c--;
			  } else { /* replacing ion! */
			    replaced=1;
			    TargetReplacingIons[target_cell]++;
			    repl_c++;
			    /* Remove vacancy: */
			    ((target_material->TargetElementalVacancies[targetIndex])[target_cell])--;
			    vac_c--;
			  }
			  replace_cell=target_cell;  /* store cell index, where replacement occured */
			  ProjState|=3; /* mark replacing projectile as sub-threshold (LSB) and as replacer (next bit) */
			  /* Transfer remaining energy of replacing particle to the phononic system:  */
			  /* no, the replacer is allowed to fly onward. It could cause some sputtering, or leave the surface itself, before it comes to rest at its replacing site  (new since version 1.0.7). */
			  /* if(store_energy_deposit==1){TargetEnergyPhonons[target_cell]+=(double)energy;}
			     energy=-0.000001; */
			}
		      }
		    } /*end of: if proj_eq_target */
		  } /* end of: if projstate=0 */
		} else { /* recoil not displaced */
		  /* energy goes to phonons: */
		  if(store_energy_deposit==1){TargetEnergyPhonons[target_cell]+=(double)recoil_energy;}
		}
	      } /* End of: Collision actually takes place */
	    } else { /* Target position is in vacuum -> no collision occurs */
	      /* Nothing */
	    }
	  } else { /* position of target nucleus would be out of sim volume. Therefore, no collision will occur */
	    /* Nothing */
	  }
	} /* end of for(coll_counter<3) */

      } else {  /* energy is smaller than min_energy, no collision happens, projectile might be stopped */
	isvac2=new_material->Is_Vacuum;
	if(isvac2==0){   /* no vacuum, stop, because of too little energy */
	  if(is_ion==1){                  /* its an ion that has to stop */
	    TargetImplantedIons[new_cell]+=1;
	    if(store_range3d>=1){fprintf(store_range3d_fp,"%g\t%g\t%g\n",x,y,z);}
	  } else {                        /* its a recoil that has to stop */
	    if( (ProjState&1)==0){        /* projectile was not sub-threshold, so it was really displaced. It's allowed to become an interstitial */
	      if(replaced==0){            /* replacing recoils are not counted as interstitials */
		((ListOfMaterials[OrgMaterial]).TargetImplantedRecoilsInt[OrgElement])[new_cell]+=1;
		int_c++;
	      }
	    } else { /* It was sub-threshold: it jumps back to its original place, no interstitial created */
	      /*if(replaced==1){
		printf("Warning! subthres recoil wants to replace atom.");
		}*/
	    }
	    if(store_recoil_cascades==1){ /* if cascades are to be stored, store empty line in file to separate record */
	      fprintf(recoil_cascades_fp,"\n");
	    }
	  }
	  if(store_energy_deposit==1){TargetEnergyPhonons[new_cell]+=(double)energy;} /* remaining projectile energy is released to phonons */
	  energy=-.00001; /* set the energy slightly negative to stop the projectile. */
	} else { /* in vacuum: may proceed */
	}
      }
    }     /* end of projectile has left simulation volume or not */

    cell_i=new_cell;
    current_material_index=new_material_index;
    current_material=new_material;
  }       /* end of while(energy>0) */
  return 0;
}

int CalcSurfaceNormal(int old_cell, int new_cell, float* nx, float* ny, float *nz){
  /* if a projectile moves from one cell to another, it is necessary to know the surface
     normal between the two cells. This function calculates the nomalized surface vector
     for two cell indices. 
     For rectangular cells the surface normal should consist of integers (1s or 0s),
     but for general geometries it might differ, so we will allow float values */
  /* function returns values via the normal vector nx, ny, nz - which will be normalized */
  int ox,oy,oz;  /* integer cell coordinates of the old cell*/
  int x,y,z;     /* integer coords of new cell */
  float temp;   /* */
  /* Get the coords from cell index: */
  GetTargetXYZ(old_cell,&ox,&oy,&oz);
  GetTargetXYZ(new_cell,&x,&y,&z);
  *nx = x-ox;
  *ny = y-oy;
  *nz = z-oz;
  if(boundary_x==1){
    if(*nx>1){ /* assume proj flew from lowest cell 0 to uppermost cell */
      *nx=-1;
    }else{
      if(*nx<-1){ /* assume proj flew from uppermost cell to lowest cell */
	*nx=+1;
      }
    }
  }
  if(boundary_y==1){
    if(*ny>1){ /* assume proj flew from lowest cell 0 to uppermost cell */
      *ny=-1;
    }else{
      if(*ny<-1){ /* assume proj flew from uppermost cell to lowest cell */
	*ny=+1;
      }
    }
  }
  if(boundary_z==1){
    if(*nz>1){ /* assume proj flew from lowest cell 0 to uppermost cell */
      *nz=-1;
    }else{
      if(*nz<-1){ /* assume proj flew from uppermost cell to lowest cell */
	*nz=+1;
      }
    }
  }
  /* normalize n-vector */
  temp=1.0/sqrtf(*nx * *nx + *ny * *ny + *nz * *nz);
  //  if((temp>1.001)||(temp<0.999)){
  //    printf("DEBUG %s, l %i, leave csn. Nvec: %g   %g   %g\n",__FILE__,__LINE__,*nx,*ny,*nz);
  //  }
  *nx *= temp;
  *ny *= temp;
  *nz *= temp;
  //  printf("DEBUG %s, l %i, leave csn. Nvec: %g   %g   %g\n",__FILE__,__LINE__,*nx,*ny,*nz);
  return 0;
}

int RefractProjectile2(float *vx, float* vy, float* vz, float nx, float ny, float nz, double* energy, float E_surf){
  /* performs surface refraction of a projectile.
     The direction and energy (passed as refs) are adjusted.
     n is the surface normal vector and E_surf the surface binding energy.
     The surface vector must point to vacuum!
     returns 0 if projectile went into solid.
     returns 1 if projectile tried to leave solid and succeded.
     returns 2 if projectile tried to leave solid but did not succeed. */
  /* Idea of this routine: decompose velocity unit vector into normal and inplane component,
     with respect to the local surface. Use normal component to add/subtract surface binding
     energy, recalculate normal component of velocity, compose new velocity vector. */ 
  float vnx, vny, vnz;  /* velocity component normal to surface */
  float vlx, vly, vlz;  /* velocity component in surface */
  float vlx2,vly2,vlz2; /* velocity component in surface after refraction */
  float vnx2,vny2,vnz2; /* velocity component normal to surface after refraction */
  float vnabs2, vlabs2; /* squared lengths of vn and vl */
  float nv;             /* scalar product of vn and n */
  double Energy2;       /* projectile energy after refraction */
  double En1, En2;      /* Energy normal to surface after and before refraction */
  double El;            /* Energy parallel to surface */
  float temp;
  int   result=10;      /* To return status */

  /* Decompose velocity into normal and in-plane component: */
  nv=(*vx * nx + *vy * ny + *vz * nz); /* inner product */
  vnx=nx*nv; vny=ny*nv; vnz=nz*nv;
  vlx = *vx - vnx; vly = *vy - vny; vlz = *vz - vnz;

#ifdef DEBUG_MODE
  if(ion_c>MONITOR_ION){  printf("DEBUG %s, line %i, vuv: %g.\n",__FILE__,__LINE__,nv);}
#endif

  /* Calculate |v|^2 for both components */
  vnabs2=(vnx*vnx + vny*vny + vnz*vnz);
  vlabs2=(vlx*vlx + vly*vly + vlz*vlz);

  if(nv>0.99999){ /* velocity is practically along surface vector, particle tries to move out */
    if(*energy>E_surf){ /* can get out */
			/* velcoity unit vector not changed, energy reduced */
      *energy -= (double)E_surf;
      return 1;
    } else { /* energy too small */
      /* reverse particle */
      *vx = -*vx; *vy = -*vy; *vz = -*vz;
      return 2;
    }
  } else {
    if(nv<-0.99999){ /* velocity is along surface vector, particle moves in */
      /* velocity unit vector doesn't change, energy increases */
      *energy += E_surf;
      return 0;
    }
    if(fabs(nv)<0.000001){ /* velocity is practically inplane */
      /* vnabs2 must not be zero. So let it be small */
      vnabs2=0.000001;
    }
  }

  /* Calculate energies: */
  En1 = *energy * (double)vnabs2;  /* normal energy */
  El  = *energy * (double)vlabs2;  /* in-plane energy */

  if(nv<0){ /* projectile moves into solid */
    result=0;
    En2     = En1 + (double)E_surf; /* gains energy */
    Energy2 = En2 + El;     /* total energy */
    temp = sqrtf(En2/(Energy2*vnabs2));
    vnx2 = vnx * temp;  vny2 = vny * temp;  vnz2 = vnz * temp;  /* direction of vn stays the same, length changes */
    temp = sqrtf(El/(Energy2*vlabs2));
    vlx2 = vlx * temp;  vly2 = vly * temp;  vlz2 = vlz * temp;  /* direction of vl stays the same, length changes */
  } else { /* projectile wants to escape, n and v in same direction */ 
    if(En1>E_surf){ /* enough energy to escape */
      result=1;
      En2     = En1 - (double)E_surf; /* looses energy */
      Energy2 = En2 + El;     /* total energy */
      temp = sqrtf(En2/(Energy2*vnabs2));
      vnx2 = vnx * temp;  vny2 = vny * temp;  vnz2 = vnz * temp;  /* direction of vn stays the same, length changes */
      temp = sqrtf(El/(Energy2*vlabs2));
      vlx2 = vlx * temp;  vly2 = vly * temp;  vlz2 = vlz * temp;  /* direction of vl stays the same, length changes */
    }else{          /* Cannot escape */
      result  = 2;
      Energy2 = *energy;
      /* Reflect on surface: */
      vnx2 = -vnx; vny2 = -vny; vnz2 = -vnz;
      vlx2 =  vlx; vly2 =  vly; vlz2 =  vlz; 
    }
  }
  /* Compose new velocity vector */
  *vx = vnx2 + vlx2;
  *vy = vny2 + vly2;
  *vz = vnz2 + vlz2;
  *energy = Energy2;   /* update energy */
  return result;

  /* Think about how to make this faster */
}


float ElectronicStopping(int ionZ, float ionM, float ionE, int material){
  /* Calculates the electronic dedx for ion with Z, M, energy Z in given material [eV/A] */
  /* material version */
  float amuRatio;
  amuRatio=MostAbundantIsotope[ionZ]/ionM; /* energy-velocity correction for isotopes which are not the most abundant. */
  /* TODO precalculate this, because division needs time! */
  /* Or change to inverse mass */
  return (ListOfMaterials[material].StoppingZE[ionZ])[Dindex(ionE*amuRatio)];  /* Energy * amuRatio. This ratio is Z/M */
}
float ElectronicStoppingdpass(int ionZ, float ionM ,float E, int material){
  /* Calculates the electronic dedx for ion with Z, M, energy Z in given material [eV/A] */
  /* material version */
 
  double Epu;
  double E0 = Enerdpass[0];
  double Y0 = ListOfMaterials[material].StoppingZEdpass[ionZ][0];;
  double EN = Enerdpass[MAX_STOPPING_ENTRIES_dpass-1];
  double YN = ListOfMaterials[material].StoppingZEdpass[ionZ][MAX_STOPPING_ENTRIES_dpass-1];

  Epu=E/ionM;
  /* printf("energies %f %f %f %f\n",E,Epu ,E0,EN); */
  /* printf("lossmin %f %f\n",E0,Y0); */
  /*   printf("lossmax %f %f\n ",EN,YN); */
  
  //    if (E <= Ener[0]) return dEdX[0];
  if (Epu <= Enerdpass[0]) {
    if (Epu <= 0.0) return 0.0;
    
    
    return Y0 * sqrt(Epu / E0);
  }
  if (Epu >= Enerdpass[MAX_STOPPING_ENTRIES_dpass-1]) {
    
    return YN * (EN / Epu);
  }
  
  int low = 0;
  int high = MAX_STOPPING_ENTRIES_dpass - 1;
  
  while (high - low > 1) {
    int mid = (low + high) / 2;
    if (Epu < Enerdpass[mid])
      high = mid;
    else
      low = mid;
  }

  double E1 = Enerdpass[low];
  double E22 = Enerdpass[high];
  double Y1 ;
  double Y2 ;
  Y1 = ListOfMaterials[material].StoppingZEdpass[ionZ][low];
  Y2 = ListOfMaterials[material].StoppingZEdpass[ionZ][high];
  
  if (Y1 <= 0 || Y2 <= 0)
    return Y1;
  
  double t = (log(Epu) - log(E1)) / (log(E22) - log(E1));
  
  return exp(log(Y1) + t * (log(Y2) - log(Y1)));
}
 
float ElectronicStoppingfile(int ionZ, float E, int material){
  /* Calculates the electronic dedx for ion with Z, M, energy Z in given material [eV/A] */
  /* material version */
 
  double Epu;
  double E0 = ListOfMaterials[material].Efile[ionZ][0];
  double Y0 = ListOfMaterials[material].StoppingZEfile[ionZ][0];;
  int MSE;
  MSE=	ListOfMaterials[material].StoppingNLfile[ionZ];
    
  double EN = ListOfMaterials[material].Efile[ionZ][MSE-1];
  double YN = ListOfMaterials[material].StoppingZEfile[ionZ][MSE-1];

  Epu=E;
  /* printf("energies %f %f %f %f\n",E,Epu ,E0,EN); */
  /* printf("lossmin %f %f\n",E0,Y0); */
  /*   printf("lossmax %f %f\n ",EN,YN); */
  
  //    if (E <= Ener[0]) return dEdX[0];
  if (Epu <= E0) {
    if (Epu <= 0.0) return 0.0;
    
    
    return Y0 * sqrt(Epu / E0);
  }
  if (Epu >= EN) {
    
    return YN * (EN / Epu);
  }
  
  int low = 0;
  int high = MSE - 1;
  
  while (high - low > 1) {
    int mid = (low + high) / 2;
    if (Epu < ListOfMaterials[material].Efile[ionZ][high])
      high = mid;
    else
      low = mid;
  }

  double E1 = ListOfMaterials[material].Efile[ionZ][low];
  double E22 = ListOfMaterials[material].Efile[ionZ][high];
  double Y1 ;
  double Y2 ;
  Y1 = ListOfMaterials[material].StoppingZEfile[ionZ][low];
  Y2 = ListOfMaterials[material].StoppingZEfile[ionZ][high];
  
  if (Y1 <= 0 || Y2 <= 0)
    return Y1;
  
  double t = (log(Epu) - log(E1)) / (log(E22) - log(E1));
  
  return exp(log(Y1) + t * (log(Y2) - log(Y1)));
}
 
  


float ElectronicStraggling(int ionZ, float ionM, float ionE, int material){
  /* Calculates the electronic energy loss straggling for ion with Z, M, energy Z in given material  [eV/A] */
  /* material version */
  float amuRatio;
  amuRatio=MostAbundantIsotope[ionZ]/ionM; /* energy-velocity correction for isotopes which are not the most abundant. */
  return (ListOfMaterials[material].StragglingZE[ionZ])[Dindex(ionE*amuRatio)];  /* Energy * amuRatio. This ratio is Z/M */
}
int CheckAndCorrectBoundary(float* dim, float target_size, float* target_size_max, int boundary){
  /* Check wheter a position is inside the target for dimension x,y or z. If
     not, the position is corrected according to the boundary conditions. 0 is
     returned if new position is inside target, 1 is returned if position is
     still outside the target (projectile has left the target) */
  int result=0;
  if((*dim)>=target_size){
    if(boundary==1){ /* PBC */
      (*dim)-=((int)((*dim)/target_size))*target_size;
    } else { /* no PBC, leave */
      result=1;
    }
  }
  if((*dim)<0){
    if(boundary==1){ /* PBC */
      (*dim)-=(((int)((*dim)/target_size))-1)*target_size;  /* this looks more complicated than dim=targetsize-dim.
							       But it is safe even if dim is < -targetsize. */
      /* Note: the previous command may sometimes results in dim being exactly equal to the target size due to
	 limited floating point precision. We therefore include a safety branch and set the value of dim to
	 slighther smaller than the target_size */

      if( (*dim)>(*target_size_max) ){
	/*	printf("dim: %f, targt: %f, maxt: %f,\n dim: %x, targ %x, maxt: %x\n",(*dim),target_size,(*target_size_max),*((int*)dim),*((int*)(&(target_size))),*((int*)target_size_max)); */
	(*dim)=(*target_size_max);
	/* 	printf("Warning! Float precision insufficient for exact boundary\ncrossing in this case, adjusting position.\n");fflush(stdout);	printf("safepos: %f, hex: %x\n",*target_size_max,*((int*)(target_size_max))); */
      }
    } else { /* no PBC, leave */
      result=1;
    }
  }
  return result;
}

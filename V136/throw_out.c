int WriteNDoubleArraysToFile(char* Filename,int ArrayCount, int ElemCount, double ** SourceArrays, char * ColumnHeads[], int FileType); /* Write multiple double arrays into one file */


int WriteNDoubleArraysToFile(char* Filename,int ArrayCount, int ElemCount, double ** SourceArrays, int** SourceArrays, char * ColumnHeads[], int FileType){      /* Writes designated array into file */
	/* Writes the designated arrays (n=ArrayCount) of ElemCount elements into a file. 
	File can have n value columns; then set FileType to 1.
	If the file should contains 3 coordinate columns (x, y, z) and n value columns, then set it to 0.
	The caller needs to make sure, that the arrays exist and are large enough. */

	FILE* fp;
	int   i,j; /* i indicates Elemens, j indicates arrays */
	int   x,y,z;

	fp=fopen(Filename,"wt");
	if(fp==NULL){return -5006;}
	i=0;

	switch(FileType){
	case 1: /* just data value column */
		fprintf(fp,"#");
		for(j=0;j<ArrayCount;j++){ /* Column heads : */
			fprintf(fp,"%s\t",ColumnHeads[j]);
		}
		fprintf(fp,"\n");
		for(i=0;i<ElemCount;i++){
			for(j=0;j<ArrayCount;j++){
				fprintf(fp,"%lg\t",((SourceArrays[j])[i])*unit_conversion_factor);
			}
			fprintf(fp,"\n");
		}
		break;
	default: /* coordinate columns x,y,z and data values. */
		fprintf(fp,"#x\ty\tz\t");
		for(j=0;j<ArrayCount;j++){ /* Column heads : */
			fprintf(fp,"%s\t",ColumnHeads[j]);
		}
		fprintf(fp,"\n");
		for(i=0;i<ElemCount;i++){
			GetTargetXYZ(i,&x,&y,&z); /* convert linear index to 3 coords */
			fprintf(fp,"%i\t%i\t%i\t",x,y,z);
			for(j=0;j<ArrayCount;j++){
				fprintf(fp,"%lg\t",((SourceArrays[j])[i])*unit_conversion_factor);
			}
			fprintf(fp,"\n");
		}
		break;
	}
	fclose(fp);
	return 0;
}


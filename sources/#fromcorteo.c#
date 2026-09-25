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


/*********************************************************************/
/* This file contains functions and declarations, taken or adapted with
small changes from Corteo (version Corteo20090527). Corteo was
written by Francois Schiettekatte.
Corteo was released under the GNU General Public License as
published by the Free Software Foundation, version 3.

You may obtain the original corteo source code from:
http://www.lps.umontreal.ca/~schiette/index.php?n=Recherche.Corteo

Modifications of the functions were done by C. Borschel.
Significant modifications are marked as such.

May 2014: some changes have been made in order to adapt from 32 to 64 bit
(some long ints have been replaced by ints like in corteo20130715).
*/
/*********************************************************************/


#include "fromcorteo.h"

/*************** Functions adapted from corteo.c ***********************/

/* fill the table of cos and sin of the scattering angle in the lab frame given a mass ratio */
void fillCosSinTable(float * cosTable, float * sinTable, float mr) {

	double sin2thetaby2, costheta, costhetaLab, sinthetaLab;
	unsigned long i;

	/* compute scattering angle components */
	for(i=0; i<DIME*DIMS; i++) {
		sin2thetaby2 = Matrix(i);
		costheta = 1.-2.*sin2thetaby2;

		if(costheta==-1.0 && mr==1.0) {
			costhetaLab = 0.0;  /* peculiar case of head-on collision of identical masses */
			sinthetaLab = 1.0;
		} else {
			costhetaLab = (costheta+mr)/sqrt(1.+2.*mr*costheta+mr*mr);
			sinthetaLab = sqrt(1.-costhetaLab*costhetaLab);
		}
		cosTable[i] = d2f(costhetaLab);
		sinTable[i] = d2f(sinthetaLab);

		/* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
		/* In some rare cases, when cos=1, then sin becomes "Not a Number". To prevent this, I will set the sine to 0 in those cases. */
		if(! ((sinTable[i]>=-1)&&(sinTable[i]<=1)) ) {
			sinTable[i]=0;
		}
	}
}

/*  Compute new directional cosines from cos & sin components of scattering angle theta
	and a random azimutal angle (no trigonometric functions called!)
	This function is critical in terms of computation efficiency */
void rotate(float *l, float *m, float *n, unsigned int * iazimAngle, float costheta, float sintheta) {
	float k, kk, kinv;
	float ll=*l, mm=*m, nn=*n;
	float k2 = 1.0f-nn*nn;
#ifdef DEBUG_MODE
	message_debug(__FILE__,__LINE__,"k2: %g\n",k2);
#endif
	/* random azimutal rotation angle components */
	float cosomega = cosAzimAngle[*iazimAngle];
	float sinomega = sinAzimAngle[*iazimAngle];
	if(++*iazimAngle==MAXAZILIST) *iazimAngle = 0;

	if(k2<=0) {	// originally: if(k2==0).  Note, Borschel 2019: in very rare cases, when v
			// is not normalized and vz is larger than 1, k2 will become slightly negative.
			// the sqrt will then produce nan, and the program will crash in transport.c
			// during determination of the cell index. Using k2<=0 avoids this problem. 
		/* extremely rare case */
		*n = costheta;
		*m = sintheta*cosomega;
		*l = sintheta*sinomega;
		return;
	} 

	/* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
	/* Errors occur in rare cases. Using the "real" sqrt is safer */
#ifndef SAFE_SQR_ROTATION
	kinv = myInvSqrt(k2);  /* 1/sqrt() (approximate) */
	k = k2*kinv;           /* so we get k and 1/k by a table lookup and a multiplication */
#else
	/* the last two lines can be replaced by the following two lines but... */
	k  = sqrtf(k2);   /*  ...using a sqrt() here makes the program 25% slower! */
	kinv = 1.0f/k;
#endif

#ifdef DEBUG_MODE
	message_debug(__FILE__,__LINE__,"k: %g, kinv: %g\n",k,kinv);
#endif

	kk = sintheta*kinv;
	*l = ll * costheta+kk*(ll*nn*cosomega+mm*sinomega);
	*m = mm * costheta+kk*(mm*nn*cosomega-ll*sinomega);
	*n = nn * costheta-k*sintheta*cosomega;

	/* MODIFICATION FROM CORTEO FOR IRADINA, C. Borschel 2011: */
	/* makes iradina slower, but safer */
#ifdef SAFE_ROTATION
	if(*l>1){*l=1};
	if(*l<-1){*l=-1};
	if(*m>1){*m=1};
	if(*m<-1){*m=-1};
	if(*n>1){*n=1};
	if(*n<-1){*n=-1};
#endif
}


/* generate and randomize lists randomlist,sqrtloglist, sinAzimAngle & cosAzimAngle   */
void computelists() {
	unsigned int irandomlist, iloglist, iazimAngle;
	/* generate a uniformly spaced list of values between .5/MAXRANLIST and 1-.5/MAXRANLIST */
	for(irandomlist=0; irandomlist<MAXRANLIST; irandomlist++)
	randomlist[irandomlist] = (irandomlist+0.5f)/MAXRANLIST;
	randomizelist(randomlist, MAXRANLIST); /* put the list in random order */
	for(irandomlist=0; irandomlist<MAXRANLIST; irandomlist++)
	sqrtrandomlist[irandomlist] = sqrtdf(randomlist[irandomlist]); /* also compute sqrt of these values */

	/* produce a list of MAXLOGLIST sqrt(-log(x)) values for x uniformly distributed */
	/* ...between x=.5/MAXLOGLIST and x=1-.5/MAXLOGLIST */
	for(iloglist=0; iloglist<MAXLOGLIST; iloglist++)
	sqrtloglist[iloglist] = sqrtdf(-log((iloglist+.5)/MAXLOGLIST));
	randomizelist(sqrtloglist, MAXLOGLIST); /* put the list in random order */
	/* precompute 1/sqrtloglist[] */
	for(iloglist=0; iloglist<MAXLOGLIST; iloglist++)
	sqrtloglist1[iloglist] = 1.0f/sqrtloglist[iloglist];

	/* produce a list uniformly distributed but randomly ordered azimutal angles  */
	for(iazimAngle=0; iazimAngle<MAXAZILIST; iazimAngle++)
	/* cosAzimAngle temporarly contains angles */
	cosAzimAngle[iazimAngle]= 2.0f*(float)PI*(float)(iazimAngle)/(float)(MAXAZILIST);
	randomizelist(cosAzimAngle, MAXAZILIST); /* put the list in random order */
	for(iazimAngle=0; iazimAngle<MAXAZILIST; iazimAngle++) {
		/* compute the cos and sine of these angles */
		sinAzimAngle[iazimAngle] = d2f(sin(cosAzimAngle[iazimAngle]));
		cosAzimAngle[iazimAngle] = d2f(cos(cosAzimAngle[iazimAngle]));
	}
}


/*************** Adapted from constants.c ***********************/

/* Extracted from NIST data in May 2007 (physics.nist.gov/PhysRefData/Compositions/index.html)
Developers and Contributors:
J. S. Coursey, D. J. Schwab, and R. A. Dragoset 
NIST, Physics Laboratory, Office of Electronic Commerce in Scientific and Engineering Data
(There are 100 data but only 92 are used with SRIM 2006) */

float MostAbundantIsotope[]={ 92,
	1.0078250321f,    4.0026032497f,       7.0160040f,       9.0121821f,      11.0093055f,
	12.0000000f,   14.0030740052f,   15.9949146221f,     18.99840320f,   19.9924401759f,
	22.98976967f,     23.98504190f,     26.98153844f,   27.9769265327f,     30.97376151f,
	31.97207069f,     34.96885271f,    39.962383123f,      38.9637069f,      39.9625912f,
	44.9559102f,      47.9479471f,      50.9439637f,      51.9405119f,      54.9380496f,
	55.9349421f,      58.9332002f,      57.9353479f,      62.9296011f,      63.9291466f,
	68.925581f,      73.9211782f,      74.9215964f,      79.9165218f,      78.9183376f,
	83.911507f,      84.9117893f,      87.9056143f,      88.9058479f,      89.9047037f,
	92.9063775f,      97.9054078f,       97.907216f,     101.9043495f,      102.905504f,
	105.903483f,      106.905093f,     113.9033581f,      114.903878f,     119.9021966f,
	120.9038180f,     129.9062228f,      126.904468f,     131.9041545f,      132.905447f,
	137.905241f,      138.906348f,      139.905434f,      140.907648f,      141.907719f,
	144.912744f,      151.919728f,      152.921226f,      157.924101f,      158.925343f,
	163.929171f,      164.930319f,      165.930290f,      168.934211f,     173.9388581f,
	174.9407679f,     179.9465488f,      180.947996f,     183.9509326f,     186.9557508f,
	191.961479f,      192.962924f,      194.964774f,      196.966552f,      201.970626f,
	204.974412f,      207.976636f,      208.980383f,      208.982416f,      209.987131f,
	222.0175705f,     223.0197307f,     226.0254026f,     227.0277470f,     232.0380504f,
	231.0358789f,     238.0507826f,     237.0481673f,      244.064198f,     243.0613727f,
	247.070347f,      247.070299f,      251.079580f,      252.082970f,     257.0950990f  };

float AtomicMass[]={ 92,
	1.00794f,        4.002602f,           6.941f,        9.012182f,          10.811f,
	12.0107f,         14.0067f,         15.9994f,      18.9984032f,         20.1797f,
	22.989770f,         24.3050f,       26.981538f,         28.0855f,       30.973761f,
	32.065f,          35.453f,          39.948f,         39.0983f,          40.078f,
	44.955910f,          47.867f,         50.9415f,         51.9961f,       54.938049f,
	55.845f,       58.933200f,         58.6934f,          63.546f,          65.409f,
	69.723f,           72.64f,        74.92160f,           78.96f,          79.904f,
	83.798f,         85.4678f,           87.62f,        88.90585f,          91.224f,
	92.90638f,           95.94f,             98.f,          101.07f,       102.90550f,
	106.42f,        107.8682f,         112.411f,         114.818f,         118.710f,
	121.760f,          127.60f,       126.90447f,         131.293f,       132.90545f,
	137.327f,        138.9055f,         140.116f,       140.90765f,          144.24f,
	145.f,          150.36f,         151.964f,          157.25f,       158.92534f,
	162.500f,       164.93032f,         167.259f,       168.93421f,          173.04f,
	174.967f,          178.49f,        180.9479f,          183.84f,         186.207f,
	190.23f,         192.217f,         195.078f,       196.96655f,          200.59f,
	204.3833f,           207.2f,       208.98038f,            209.f,            210.f,
	222.f,            223.f,            226.f,            227.f,        232.0381f,
	231.03588f,       238.02891f,            237.f,            244.f,            243.f,
	247.f,            247.f,            251.f,            252.f,            257.f  };


char  AtomicNames[][3]= {"XX",
"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", 
"K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
"Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", 
"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", 
"Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
"Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf"};



/*************** Adapted from corteoindex.c ***********************/

/* Xindex: compute an integer index using the binary representation of a floating point value val/MIN
(X is E for reduced energy, S for reduced impact parameter or D for energy of stopping).
These indexes are used to access sin/cos tables of the scattering angle and stopping power tables.
The division by MIN is actually carried out by subtracting the appropriate BIAS from the resulting index

Input: val, with val >= MIN, (MIN = MINE, MINS, or MIND)

Returns an integer equal to 16 times the exponent base 2 of val 
+ 1 for each 1/16 interval between each power of 2
(Index of input value = MIN is 0)

examples
val/MIN    index
1           0
1.0625      1
1.0626      1
1.1249      1
1.1250      2
1.2         3
1.25        4
2           16
5           40

BEWARE: architecture dependent!
Assuming 32-bit 'unsigned long' integer and float
Assuming IEEE Standard 754 representation for float:
bit 31: sign
bit 30-23: exponent (base 2), biased by 127 (i.e. value of exponent for 2^0 is 127)
bit 22-0: mantissa
*/

/* version that compute the index of a reduced energy */
unsigned int Eindex(float val) {
	unsigned int ll;
	ll = *(unsigned int *)&val;
	ll = (ll >> SHIFTE)-BIASE;

#ifdef INDEX_BOUND_CHECKING
	if(val<MINE) {
		EminErr++;
		return 0;
	}
	if(ll>=DIME) {
		EmaxErr++;
		return DIME-1;
	}
#endif

	return ll;
}

/* version that compute the index of a reduced impact parameter */
unsigned int Sindex(float val) {
	unsigned int ll;
	ll= *(unsigned int *)&val;
	ll = (ll >> SHIFTS)-BIASS;

#ifdef INDEX_BOUND_CHECKING
	if(val<MINS) {
		SminErr++;
		return 0;
	}
	if(ll>=DIMS) {
		SmaxErr++;
		return DIMS-1;
	}
#endif

	return ll;
}

/* version that compute the index of the energy of a stopping power and enegy straggling */
unsigned int Dindex(float val) {
	unsigned int ll;
	ll = *(unsigned int *)&val;
	ll = (ll >> SHIFTD)-BIASD;

#ifdef INDEX_BOUND_CHECKING
	if(val<MIND) {
		DminErr++;
		return 0;
	}
	if(ll>=DIMD) {
		DmaxErr++;
		return DIMD-1;
	}
#endif

	return ll;
}

/* return the float value corresponding to index
(Actually returns the average value between the value of this index and the next
because when an index is computed, values are truncated to the lower bound of the interval.)
*/
float Eval(unsigned int index) {
	float temp1, temp2;
	unsigned int ll;

	ll = (index+BIASE)<<SHIFTE;
	temp1 = (*(float *)&ll);
	ll = (index+1+BIASE)<<SHIFTE;
	temp2 = (*(float *)&ll);

	return (temp1+temp2)*0.5f;
}

float Sval(unsigned int index) {
	float temp1, temp2;
	unsigned int ll;

	ll = (index+BIASS)<<SHIFTS;
	temp1 = (*(float *)&ll);
	ll = (index+1+BIASS)<<SHIFTS;
	temp2 = (*(float *)&ll);

	return (temp1+temp2)*0.5f;
}

float Dval(unsigned int index) {
	float temp1, temp2;
	unsigned int ll;

	ll = (index+BIASD)<<SHIFTD;
	temp1 = (*(float *)&ll);
	ll = (index+1+BIASD)<<SHIFTD;
	temp2 = (*(float *)&ll);

	return (temp1+temp2)*0.5f;
}



/*************** Adapted from corteomatrix.c ***********************/

/* scattering matrix: sin^2(theta_CM/2) */
float matrix[DIME*DIMS];

/* matrix elements calulation parameters */
#define NSUM  1000   /* number of terms in Gauss-Mehler quadrature sum when computing the matrix */
#define NSUM2 100    /* number of terms in Gauss-Mehler quadrature sum when evaluation the screened cross-section */
#define THETAERR -1000.0

/* parameters used to compute matrix, stored in matrix file */
#define HEADERSIZE 7
float headerRef[HEADERSIZE] = {MINE, DIME, MAXE, MINS, DIMS, MAXS, NSUM}; 


// pointer to screening function (I like C!)
double (*PHI)(double x);

// Universal screening function 
double PHIUniv(double x) {
	return 0.1818*exp(-3.*x)+0.5099*exp(-0.9423*x)+0.2802*exp(-0.4028*x)+0.02817*exp(-0.2016*x);
}

// Lenz-Jensen screening function
double PHILJ(double x) {
	double y = 3.108*sqrt(x);
	return exp(-y)*(1+y+0.3344*y*y+0.0485*y*y*y+2.647e-3*y*y*y*y);
}

// no screening (screening = 1)
double PHInone(double x) {
	return 1.0;
}


// Following functions used to compute scattering angle following Gauss-Mehler quadrature
double H(double u, double x0, double epsilon, double s) {
	double x0u = x0/u;
	double eps1 = 1./epsilon;
	double ss = s*s;
	return sqrt((1-u*u)/(1-PHI(x0u)/x0u*eps1-ss/x0u/x0u));
}

// the root of this function (vs x0) is the minimal 
// approach distance at energy epsilon and impact parameter s
double funcX0(double x0, double epsilon, double s) {
	double sx0 = s/x0;
	return 1.-PHI(x0)/(x0*epsilon)-sx0*sx0;
}

// use bisection (Newton) method to find the minimal approach distance x0
// (return THETAERR on error)
double findX0(double epsilon, double s) {
	double funcX0val1, funcX0val2, funcX0valm;

	double temp = 1.0/(2.0*epsilon);
	double x2 = temp+sqrt(temp+s*s); // inital guesses: Mendenhall & Weller NIMB 58(1991)11, eq. 15 

	double x1 = x2/10.;
	funcX0val1 = funcX0(x1, epsilon, s);  // should be always negative
	funcX0val2 = funcX0(x2, epsilon, s);  // should be always positive (starting ~1.0)
	funcX0valm = funcX0((x1+x2)/2., epsilon, s);
	if(funcX0val1>=0.) {
		// initial guess for x1 too optimistic, start with a safe value (much longer)
		x1 = MINS;
		funcX0val1 = funcX0(x1, epsilon, s);  // should be always negative
	}

	if(funcX0val1>=0. || funcX0val2<=0.)
	return THETAERR;    // error, values should be on each side of 0

	do {
		if(funcX0valm<0.) {
			funcX0val1 = funcX0valm;
			x1 = (x1+x2)/2.;
		} else {
			funcX0val2 = funcX0valm;
			x2 = (x1+x2)/2.;
		}
		funcX0valm = funcX0((x1+x2)/2., epsilon, s);
	} while (fabs(funcX0valm)>1e-14);
	// 1.-XXX wont be more precise than 2^(-52)=2e-16
	// but using 1e-10 saves time and is enough precise for float conversion later on

	return (x1+x2)/2.;
}

// use bisection method to find the reduced 
// impact parameter s given epsilon and thetaCM 
// (return THETAERR on error)
double finds(double epsilon, double thetaCM) {
	double funcX0val1, funcX0val2, funcX0valm;

	// inital guesses: Mendenhall & Weller NIMB 58 (1991) 11, eqs. 23-25
	double gamma = (PI-thetaCM)/PI;
	double x0 = findX0((1.0-gamma*gamma)*epsilon, MINS);
	double x1 = 0.7*gamma*x0;
	double x2 = 1.0/(2.0*epsilon*tan(thetaCM/2.0)); 

	funcX0val1 = thetaCM-THETA(epsilon, x1, NSUM2);  // should be always negative
	funcX0val2 = thetaCM-THETA(epsilon, x2, NSUM2);  // should be always positive (starting ~1.0)
	funcX0valm = thetaCM-THETA(epsilon, (x1+x2)/2., NSUM2);

	if(funcX0val1>=0. || funcX0val2<=0.)
	return THETAERR;    // error, values should be on each side of 0

	do {
		if(funcX0valm<0.) {
			funcX0val1 = funcX0valm;
			x1 = (x1+x2)/2.;
		} else {
			funcX0val2 = funcX0valm;
			x2 = (x1+x2)/2.;
		}
		funcX0valm = thetaCM-THETA(epsilon, (x1+x2)/2., NSUM2);
	} while (fabs(funcX0valm)>1e-5);

	return (x1+x2)/2.;
}

// scattering angle calculated using Gauss-Mehler quadrature
// (return THETAERR on error)
double THETA(double epsilon, double s, unsigned int nsum) {
	double temp, sum = 0.;//, sum1 = 0;
	double x0 = findX0(epsilon, s);
	unsigned int j;
	double theta;

	if(x0==THETAERR) {
		return THETAERR;
	}

	for(j=0; j<nsum/2; j++) {
		temp = H(cos((2.*j-1.)*PI/(2.*nsum)),x0, epsilon, s);
		sum += temp;
	}
	theta = 2.*PI*s/nsum/x0*sum;

	return PI-theta;

}

/* compute all the elements of matrix and write matrix to file 'corteo.mat'
user sets showProgress!=0 to display the progress of this (long) calculation to the console
return 1 if successful, 0 if not able to write file */
int calcMatrix(int showProgress, char* DirectoryData) {
	unsigned long i, j, nThetaErr = 0;
	double theta, sinThetaBy2;
	//double s, ds, dsdTheta, eps;
	//FILE *off;

	FILE *ofp;
	char *mode = "wb";
	char dataFile[MAX_FILENAME_LENGTH];

	sprintf(dataFile, "%s/corteo.mat", DirectoryData); /* dir + corteo.mat */
	fprintf(stdout, "create %s", dataFile);

	PHI = PHIUniv;  // matrix computed using Universal potential 
	// (TODO: let the user decide; WARNING then: screening ion->a to be computed with the right screening length)

	if(showProgress) fprintf(stdout, "Computing scattering matrix");
	// compute matrix for each reduced energy, reduced impact parameter pair

	for(j=0; j<DIME; j++) {
		if(showProgress && j%(DIME/10)==0) {
			fprintf(stdout, ".");
			fflush(stdout);
		}
		for(i=0; i<DIMS; i++) {
			// calculations (summation) made using double to decrease numerical noise
			theta = THETA(Eval(j), Sval(i), NSUM);
			if(theta == THETAERR) nThetaErr++;
			sinThetaBy2 = sin(theta/2.);

			// store in matrix sin^2(theta_CM/2) as float 
			setMatrix(i+j*DIMS,(float)(sinThetaBy2*sinThetaBy2));
		}
	}

	// open file corteo.mat and write matrix, including a header indicating the parameters used to compute matrix
	ofp = fopen(dataFile, mode); // ex "./data/corteo.mat"
	if (ofp == NULL) return 0;
	fwrite((void*)headerRef, HEADERSIZE, sizeof(float), ofp);
	fwrite((void*)matrix,    DIME*DIMS,  sizeof(float), ofp);
	fclose(ofp);

	if(showProgress) { 
		fprintf(stdout, "done,\n");
		fflush(stdout);
	}
	if(nThetaErr) message_error(0,"ERROR: %lu error(s) evaluating theta.\n", nThetaErr);

	return 1;
}

/* read scattering matrix from file corteo.mat
return 1 if file read correctly, 0 otherwise */
int loadMatrix(char* DirectoryData) {
	unsigned int i;
	FILE *ifp;
	char *mode = "rb";

	float header[HEADERSIZE];
	unsigned int tail;

	char dataFile[1000];

	sprintf(dataFile, "%s/corteo.mat", DirectoryData); /* dir + corteo.mat */
	message(1,"Open %s\n", dataFile);
	ifp = fopen(dataFile, mode); // "../data/corteo.mat"

	if (ifp == NULL) {
		message_error(0,"ERROR: Can't open '%s'.", dataFile);
		return 0;
	}

	// file exists, get header
	fread((void*)header, HEADERSIZE, sizeof(float), ifp);
	if(print_level>2){
		for(i=0; i<HEADERSIZE; i++) {
			message(3,"----- Header elements %i (header %.4e and %.4e expected)\n", i, header[i], headerRef[i]);
		}
	}
	for(i=0; i<HEADERSIZE; i++)
	if(header[i] != headerRef[i]) {
		/* "corteo.mat" header does not correspond to currently defined parameters */
		fclose(ifp);
		message_error(0,"ERROR: Header element %i does not correspond to current parameters in corteo.mat\n",i);
		return 0;
	}
	message(2,"----- Header elements correspond to current parameters in corteo.mat\n");

	// read file into matrix
	fread((void*)matrix, DIME*DIMS, sizeof(float), ifp);

	// read control data
	fread((void*)(&tail), 1, sizeof(int), ifp);
	fclose(ifp);

	return 1;
}


// return element i of matrix
float Matrix(unsigned long i) { 
	return matrix[i]; 
}

// set element i of matrix
void setMatrix(unsigned long i, float val) {
	matrix[i] = val; 
}


/* returns the cross section in the center of mass frame considering a screened potential */
double crossSectionScreenPot(double E, unsigned int Z1, unsigned int Z2, double massRatio, double thetaCM, unsigned int screeningType) {
	double ds, dsdTheta, s, screenLength, epsilon;

	// set screening lenth and screening function
	switch (screeningType) {
	case 0:
		screenLength = 0.5291772108; // Bohr radius
		PHI = PHInone; // unscreened potential (screening = 1)
		break;
	case 2:
		screenLength = SCREENCONST/(pow(Z1,0.23)+pow(Z2,  0.23));   // Universal screening length
		PHI = PHIUniv; // Universal screening function
		break;
	case 3:
		screenLength = SCREENCONST/sqrt(pow(Z1,2/3)+pow(Z2, 2/3));  // Lindhard screening length
		PHI = PHILJ; // Lenz-Jensen screening function
		break;
	default:
		message_error(0,"screening function of type %u unknown.", screeningType);
		exit(1);
	}

	// reduced energy according to screening length
	epsilon = E*screenLength/(Z1*Z2*E2);

	// find corresponfding reduced impact parameter 
	s = finds(epsilon, thetaCM);
	if(s==THETAERR) return 0;

	// ds/dTheta using five-point stencil 
	ds = s*0.01;
	dsdTheta = (12.0*ds)/(-THETA(epsilon,s+2.0*ds,NSUM2)+8.0*THETA(epsilon,s+ds,NSUM2)-8.0*THETA(epsilon,s-ds,NSUM2)+THETA(epsilon,s-2.0*ds,NSUM2));

	// return "non-reduced" center-of-mass cross section 
	return screenLength*screenLength*s/sin(thetaCM)*fabs(dsdTheta);

}


/*************** Adapted from corteoutil.c ***********************/
// randomly swap 3 times each elements of list so they come out in random order
void randomizelist(float *list, unsigned int maxlist) {
	unsigned int i, ilist, j;
	float x;

	for(i=0; i<maxlist*3; i++) {
		ilist = i%maxlist;
		j = (unsigned int)floor(randomx()*maxlist);
		x = list[ilist];
		list[ilist] = list[j];
		list[j] = x;
	}
}
// skip to end of line, including end of line character
void ignoreline(FILE *ifp) {
	fscanf(ifp, "%*[^\n]");
	fscanf(ifp, "%*1[\n]");
}

// convert double to float while preventing underflow
float d2f(double val) {
	if(val<0.0) return -val<FLT_MIN?0.0f:(float)val;
	return val<FLT_MIN?0.0f:(float)val;
}

// return sqrt in float format while preventing underflow
float sqrtdf(double val) {
	return d2f(sqrt(val));
}

// convert a string to float while preventing underflow
float a2f(char * s) {
	return d2f(atof(s));
}

// next 3 functions: table-based sqrt and 1/sqrt evaluation
// tables are computed for the sqrt and 1/sqrt of mantissa and exponent separately
// (one table, the shortest one related to the exponent, could be eliminated by simply 
//  bit-shifting the exponent by one bit... to be explored)

// fill tables
void mySqrtTableFill() {
	//  int safe;
	unsigned long i, j, n = 1<<16; // 1<<16: mantissa tables contain 65536 values
	float val;

	// mantissa from 0 to 2 in 65536 steps, giving off precision!
	for(i=0;i<n;i++) {
		// generate a float values 
		j = (i<<7) + ((1<<23)*127); 
		//    for(safe=0;safe<10;safe++);
		val = *(float *)&j;
		// store its sqrt and 1/sqrt in tables
		mySqrtTable[i]    = sqrtdf(val);  
		myInvSqrtTable[i] = 1.0f/mySqrtTable[i];
	}

	// exponent from 2^-128 to 2^128
	for(i=0;i<(1<<8);i++) {
		// generate a float values 
		j = i<<23;
		val = *(float *)&j;
		// store its sqrt and 1/sqrt in tables
		mySqrtTableExp[i]    = sqrtdf(val);
		myInvSqrtTableExp[i] = 1.0f/mySqrtTableExp[i];
	}
}

// sqrt of val is the product of the sqrt of its mantissa and sqrt of its exponent
// WARNING: approximate solution, precise to 0.0015%, use when precision is not critical
float mySqrt(float val) {
	if( *(unsigned long *)&val==0 ) return 0.0f;
	return mySqrtTable[((*(unsigned long *)&val)>>7)&0xFFFF]*mySqrtTableExp[((*(unsigned long *)&val)>>23)&0xFF];
}

// 1/sqrt of val is the product of the 1/sqrt of its mantissa and 1/sqrt of its exponent
// WARNING: approximate solution, precise to 0.0015%, use when precision is not critical
float myInvSqrt(float val) {
	if((*(unsigned long *)&val)==0) return 1.0f/val; // prevent division by 0
	return myInvSqrtTable[((*(unsigned long *)&val)>>7)&0xFFFF]*myInvSqrtTableExp[((*(unsigned long *)&val)>>23)&0xFF];
}


/*************** Adapted from randomx.c ***********************/

// random function from P. L'Ecuyer, Communications of the ACM 31 (1988) 742
double randomx() {
	int z,k;
	k = seed1 / 53668;
	seed1 = 40014 * (seed1 - k * 53668) - k * 12211;
	if (seed1 < 0) seed1 = seed1 + 2147483563;
	k = seed2 / 52774;
	seed2 = 40692 * (seed2 - k * 52774) - k * 3791;
	if (seed2 < 0) seed2 = seed2 + 2147483399;
	z = seed1 - seed2;
	if (z < 1) z = z + 2147483562;
	return z * 4.65661305956E-10;
}

/************** Adapted from corteo20130715 *********************/
/* Length of float and int must be 4 byte, otherwise indexing doesn't
work correctly. */
unsigned int check_type_representation() {
	float f = 3328.625f;
	if( sizeof(unsigned int)!=4 ) {
		message_error(0,"sizeof(unsigned int) = %i instead of 4.\n", (int)sizeof(unsigned int));
		return 1;
	}
	if(sizeof(float)!=4)  {
		message_error(0,"sizeof(float) = %i instead of 4.\n", (int)sizeof(float));
		return 2;
	}
	if(*(unsigned int *)&f != 0x45500a00)  {
		message_error(0,"looks like your machine does not follow IEEE-754 standard for binary representation of float.\n");
		return 3;
	}
	return 0;
}
